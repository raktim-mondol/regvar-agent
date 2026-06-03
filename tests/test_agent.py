"""Tests for the agent loop (agent.py).

Uses mocked OpenAI client with scripted responses so tests run without
API keys or network access.
"""

import json
from unittest.mock import MagicMock, patch

import pytest

from regvar.agent import (
    DEFAULT_MODEL,
    MAX_TURNS,
    SYSTEM_PROMPT,
    _build_system_prompt,
    build_task_message,
    run_agent_turn,
)
from regvar.variants import CandidateVariant


def _make_candidates():
    return [
        CandidateVariant("chr8", 127401060, "G", "T", region_id="enh1", note="test SNP"),
        CandidateVariant("chr10", 46046326, "A", "G", region_id="enh2"),
    ]


def test_build_task_message():
    candidates = _make_candidates()
    msg = build_task_message(candidates)

    assert "Candidate regulatory variants:" in msg
    assert "chr8:127401060:G>T" in msg
    assert "enh1" in msg
    assert "test SNP" in msg
    assert "chr10:46046326:A>G" in msg
    assert "enh2" in msg
    assert "Score them" in msg


def test_build_task_message_no_metadata():
    candidates = [CandidateVariant("chr1", 100, "A", "C")]
    msg = build_task_message(candidates)

    assert "chr1:100:A>C" in msg
    assert "()" not in msg


def test_build_system_prompt_no_hints():
    result = _build_system_prompt()
    assert result == SYSTEM_PROMPT


def test_build_system_prompt_with_assay_hint():
    result = _build_system_prompt(assay_hint=["ATAC-seq", "RNA-seq"])

    assert SYSTEM_PROMPT in result
    assert "Preferred assays: ATAC-seq, RNA-seq" in result


def test_build_system_prompt_with_tissue_hint():
    result = _build_system_prompt(tissue_hint=["UBERON:0002367"])

    assert "Preferred tissue ontology terms: UBERON:0002367" in result


def test_build_system_prompt_with_top_n_hint():
    result = _build_system_prompt(top_n_hint=15)

    assert "Return the top 15 effects" in result


def test_build_system_prompt_all_hints():
    result = _build_system_prompt(
        assay_hint=["RNA-seq"],
        tissue_hint=["CL:0000057"],
        top_n_hint=5,
    )

    assert "Preferred assays: RNA-seq" in result
    assert "CL:0000057" in result
    assert "top 5 effects" in result


class _FakeChoice:
    def __init__(self, message, finish_reason="stop"):
        self.message = message
        self.finish_reason = finish_reason


class _FakeResponse:
    def __init__(self, choices):
        self.choices = choices


def _make_final_message(content="Here is the analysis."):
    msg = MagicMock()
    msg.content = content
    msg.model_dump.return_value = {"role": "assistant", "content": content}
    msg.tool_calls = None
    return msg


def _make_tool_call_message(tool_name, tool_args, tool_call_id="call_1"):
    tc = MagicMock()
    tc.function.name = tool_name
    tc.function.arguments = json.dumps(tool_args)
    tc.id = tool_call_id

    msg = MagicMock()
    msg.content = None
    msg.model_dump.return_value = {
        "role": "assistant",
        "content": None,
        "tool_calls": [{"id": tool_call_id, "function": {"name": tool_name, "arguments": json.dumps(tool_args)}}],
    }
    msg.tool_calls = [tc]
    return msg


@patch("regvar.agent.run_tool", return_value=json.dumps({"assays": {}}))
def test_run_agent_single_turn(mock_run_tool):
    """Agent returns a final answer immediately — no tool calls."""
    fake_msg = _make_final_message("Analysis complete.")
    mock_client = MagicMock()
    mock_client.chat.completions.create.return_value = _FakeResponse([_FakeChoice(fake_msg)])

    messages = [
        {"role": "system", "content": "test"},
        {"role": "user", "content": "do it"},
    ]
    result = run_agent_turn(messages, mock_client, reasoning_effort=None)

    assert result == "Analysis complete."
    mock_run_tool.assert_not_called()
    assert mock_client.chat.completions.create.call_count == 1


@patch("regvar.agent.run_tool", return_value=json.dumps({"result": "ok"}))
def test_run_agent_multi_turn(mock_run_tool):
    """Agent makes one tool call, then produces a final answer."""
    tool_msg = _make_tool_call_message("list_supported_assays", {})
    final_msg = _make_final_message("Done after tool call.")

    mock_client = MagicMock()
    mock_client.chat.completions.create.side_effect = [
        _FakeResponse([_FakeChoice(tool_msg, finish_reason="tool_calls")]),
        _FakeResponse([_FakeChoice(final_msg)]),
    ]

    messages = [{"role": "system", "content": "test"}]
    result = run_agent_turn(messages, mock_client, reasoning_effort=None)

    assert result == "Done after tool call."
    assert mock_run_tool.call_count == 1
    assert mock_client.chat.completions.create.call_count == 2

    tool_message = messages[-2]
    assert tool_message["role"] == "tool"
    assert tool_message["content"] == json.dumps({"result": "ok"})


@patch("regvar.agent.run_tool", return_value=json.dumps({"result": "ok"}))
def test_run_agent_max_turns(mock_run_tool):
    """Agent always returns tool calls, hitting the MAX_TURNS limit."""
    tool_msg = _make_tool_call_message("list_supported_assays", {})
    mock_client = MagicMock()
    mock_client.chat.completions.create.return_value = _FakeResponse(
        [_FakeChoice(tool_msg, finish_reason="tool_calls")]
    )

    messages = [{"role": "system", "content": "test"}]
    result = run_agent_turn(messages, mock_client, reasoning_effort=None, max_turns=3)

    assert "MAX_TURNS" in result
    assert mock_client.chat.completions.create.call_count == 3


@patch("regvar.agent.run_tool", return_value=json.dumps({"ok": True}))
def test_run_agent_on_tool_call_callback(mock_run_tool):
    """on_tool_call callback is invoked with correct arguments."""
    tool_msg = _make_tool_call_message("list_supported_assays", {}, tool_call_id="call_42")
    final_msg = _make_final_message("Done.")

    mock_client = MagicMock()
    mock_client.chat.completions.create.side_effect = [
        _FakeResponse([_FakeChoice(tool_msg, finish_reason="tool_calls")]),
        _FakeResponse([_FakeChoice(final_msg)]),
    ]

    callback_calls = []

    def _cb(name, args, result_json):
        callback_calls.append((name, args, result_json))

    messages = [{"role": "system", "content": "test"}]
    run_agent_turn(messages, mock_client, reasoning_effort=None, on_tool_call=_cb)

    assert len(callback_calls) == 1
    name, args, result = callback_calls[0]
    assert name == "list_supported_assays"
    assert args == {}
    assert json.loads(result) == {"ok": True}


@patch("regvar.agent.run_tool", return_value=json.dumps({"error": "bad json"}))
def test_run_agent_tool_call_with_bad_json(mock_run_tool):
    """Tool call with malformed arguments still dispatches (empty args)."""
    tc = MagicMock()
    tc.function.name = "list_supported_assays"
    tc.function.arguments = "not json at all"
    tc.id = "call_bad"

    tool_msg = MagicMock()
    tool_msg.content = None
    tool_msg.model_dump.return_value = {"role": "assistant"}
    tool_msg.tool_calls = [tc]

    final_msg = _make_final_message("Recovered.")

    mock_client = MagicMock()
    mock_client.chat.completions.create.side_effect = [
        _FakeResponse([_FakeChoice(tool_msg, finish_reason="tool_calls")]),
        _FakeResponse([_FakeChoice(final_msg)]),
    ]

    messages = [{"role": "system", "content": "test"}]
    result = run_agent_turn(messages, mock_client, reasoning_effort=None)

    assert result == "Recovered."
    mock_run_tool.assert_called_once_with("list_supported_assays", {})
