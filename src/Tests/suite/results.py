"""Result store of the TNL-SPH unified test suite (pure data, stdlib only).

This module MUST NOT print anything (all console output lives in
:mod:`suite.reporting`) and MUST NOT import third-party packages.

Data contract (pinned; consumed by todo 27's session log and todo 28's
``--report`` reprint):

``session_log(entries, machine)`` returns exactly::

    {"schema": 1,
     "started": <iso8601 str>,
     "duration_s": <float>,
     "machine": {"gpu": <str>, "host": <str>},
     "entries": [entry.to_json(), ...]}

``ResultEntry`` carries per-case results: ``tag``, ``group``, ``tier``,
``phases`` (mapping of canonical phase names ``"init"``/``"build"``/``"run"``
to :class:`PhaseResult`), test counts, ``comp_time`` and the two failure
channels: ``error`` (hard failure) and ``skip_reason`` (explicit skip, e.g.
an unfrozen reference -- skips are NOT failures, see
:func:`result_skipped`/:func:`print_final_summary` in ``suite.reporting``).

The markdown writer renders a GitHub-flavored table WITHOUT pandas with the
pinned columns::

    | Case tag | Group | Init | Build | Run | Comp. time | Ref. time | Time diff % | Tests | Tests result | Error |

Cell rules are ported from ``src/tools/multirunner.py:156-177``
(``format_result_html``/``format_tests_html``/``format_tests_summary_html``)
and ``src/tools/testSolver.py:149-162`` (time-diff coloring) so the output
stays compatible with ``src/tools/printTestLog.py``'s span parsing
(``<span style="color:...">__text__</span>``).

Selftest (writes only into a temp dir)::

    PYTHONPATH=src/Tests python3 -m suite.results --selftest
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping

# Canonical phase keys and their markdown column headings (pinned order).
PHASE_KEYS = ("init", "build", "run")
PHASE_HEADINGS = ("Init", "Build", "Run")

_PHASE_STATUSES = ("PASS", "FAIL", "SKIP")

# Pinned markdown columns (todo 27/28 contract).
MD_COLUMNS = (
    "Case tag",
    "Group",
    "Init",
    "Build",
    "Run",
    "Comp. time",
    "Ref. time",
    "Time diff %",
    "Tests",
    "Tests result",
    "Error",
)

# Time-diff span thresholds ported from src/tools/testSolver.py:155-160.
_TIME_DIFF_RED = 5.0
_TIME_DIFF_GREEN = -5.0


@dataclass
class PhaseResult:
    """Outcome of a single runner phase (``status`` in PASS/FAIL/SKIP)."""

    status: str
    returncode: int = 0
    seconds: float = 0.0
    log: str = ""

    def __post_init__(self) -> None:
        if self.status not in _PHASE_STATUSES:
            raise ValueError(
                f"PhaseResult: invalid status {self.status!r} "
                f"(expected one of {_PHASE_STATUSES})"
            )

    def to_json(self) -> dict:
        return {
            "status": self.status,
            "returncode": self.returncode,
            "seconds": self.seconds,
            "log": self.log,
        }

    @classmethod
    def from_json(cls, data: Mapping) -> "PhaseResult":
        return cls(
            status=str(data.get("status", "")),
            returncode=int(data.get("returncode", 0)),
            seconds=float(data.get("seconds", 0.0)),
            log=str(data.get("log", "")),
        )


@dataclass
class ResultEntry:
    """Per-case result row; ``to_json``/``from_json`` round-trip ALL fields."""

    tag: str
    group: str
    tier: str
    phases: dict[str, PhaseResult] = field(default_factory=dict)
    tests_passed: int = 0
    tests_total: int = 0
    comp_time: float = 0.0
    error: str = ""
    skip_reason: str = ""

    def to_json(self) -> dict:
        return {
            "tag": self.tag,
            "group": self.group,
            "tier": self.tier,
            "phases": {name: phase.to_json() for name, phase in self.phases.items()},
            "tests_passed": self.tests_passed,
            "tests_total": self.tests_total,
            "comp_time": self.comp_time,
            "error": self.error,
            "skip_reason": self.skip_reason,
        }

    @classmethod
    def from_json(cls, data: Mapping) -> "ResultEntry":
        return cls(
            tag=str(data.get("tag", "")),
            group=str(data.get("group", "")),
            tier=str(data.get("tier", "")),
            phases={
                str(name): PhaseResult.from_json(phase)
                for name, phase in (data.get("phases") or {}).items()
            },
            tests_passed=int(data.get("tests_passed", 0)),
            tests_total=int(data.get("tests_total", 0)),
            comp_time=float(data.get("comp_time", 0.0)),
            error=str(data.get("error", "")),
            skip_reason=str(data.get("skip_reason", "")),
        )


def result_failed(entry: ResultEntry) -> bool:
    """True when the entry carries a hard failure: error text or any FAIL phase."""
    return bool(entry.error) or any(p.status == "FAIL" for p in entry.phases.values())


def result_skipped(entry: ResultEntry) -> bool:
    """True when the entry was explicitly skipped (``skip_reason`` set).

    Skips are kept strictly separate from failures by every reporter.
    """
    return bool(entry.skip_reason)


def _span(color: str, text: str) -> str:
    """printTestLog-compatible colored-bold marker."""
    return f'<span style="color:{color}">__{text}__</span>'


# Ported from src/tools/multirunner.py:156-159 (returncode -> phase/error rule).
def format_result_html(entry: ResultEntry) -> str:
    """Green ``__Success__`` when no FAIL phase and no error, else red ``__Failed__``."""
    if result_failed(entry):
        return _span("red", "Failed")
    return _span("green", "Success")


# Ported from src/tools/multirunner.py:161-168.
def format_tests_html(entry: ResultEntry) -> str:
    """Green ``__x/y__`` when passed==total>0, ``-`` when total==0, red otherwise."""
    passed = entry.tests_passed
    total = entry.tests_total
    if passed == total:
        if total > 0:
            return _span("green", f"{passed}/{total}")
        return "-"
    return _span("red", f"{passed}/{total}")


# Ported from src/tools/multirunner.py:170-177.
def format_tests_summary_html(entry: ResultEntry) -> str:
    """Green ``__Passed__`` when passed==total>0, ``-`` when total==0, red ``__Failed__``."""
    passed = entry.tests_passed
    total = entry.tests_total
    if passed == total:
        if total > 0:
            return _span("green", "Passed")
        return "-"
    return _span("red", "Failed")


def format_phase_html(phase: PhaseResult | None) -> str:
    """Phase cell: PASS green / FAIL red / SKIP yellow / absent ``-``."""
    if phase is None:
        return "-"
    colors = {"PASS": "green", "FAIL": "red", "SKIP": "yellow"}
    return _span(colors.get(phase.status, "red"), phase.status)


# Ported from src/tools/testSolver.py:149-162 (span rules and ``.1f %`` format).
def format_time_diff(comp_time: float, ref_time: float | None) -> str:
    """Red span when diff > 5 %, green when < -5 %, plain otherwise, ``-`` without a ref."""
    if not ref_time:
        return "-"
    diff = 100.0 * (comp_time - ref_time) / ref_time
    if diff > _TIME_DIFF_RED:
        return _span("red", f"{diff:.1f} %")
    if diff < _TIME_DIFF_GREEN:
        return _span("green", f"{diff:.1f} %")
    return f"{diff:.1f} %"


def session_log(
    entries: list[ResultEntry],
    machine: Mapping | None,
    started: float | None = None,
) -> dict:
    """Build the pinned session-log document (todo 27 contract).

    ``machine`` is a mapping with ``gpu``/``host`` keys (only these two are
    kept in the document). ``started`` is an optional epoch timestamp; when
    omitted, "now" is used and ``duration_s`` is ``0.0``.
    """
    now = time.time()
    start_ts = now if started is None else float(started)
    machine = machine or {}
    return {
        "schema": 1,
        "started": datetime.fromtimestamp(start_ts, tz=timezone.utc).astimezone().isoformat(),
        "duration_s": round(now - start_ts, 3),
        "machine": {
            "gpu": str(machine.get("gpu", "")),
            "host": str(machine.get("host", "")),
        },
        "entries": [entry.to_json() for entry in entries],
    }


def _require_writable_parent(path: Path, what: str) -> Path:
    """Return ``path`` or raise RuntimeError naming it; never creates directories."""
    path = Path(path)
    parent = path.parent
    if not parent.is_dir():
        raise RuntimeError(
            f"cannot write {what} to {path}: parent directory {parent} does not exist"
        )
    if not os.access(parent, os.W_OK):
        raise RuntimeError(
            f"cannot write {what} to {path}: parent directory {parent} is not writable"
        )
    return path


def write_json_log(session: dict, path: Path) -> Path:
    """Write the session log as 2-space-indented JSON; returns the path.

    Raises ``RuntimeError`` (naming the path) when the parent directory is
    missing or unwritable -- directories are NEVER created silently.
    """
    path = _require_writable_parent(path, "json log")
    path.write_text(json.dumps(session, indent=2) + "\n", encoding="utf-8")
    return path


def _md_cell(text: str) -> str:
    """Escape a raw string for a GitHub-flavored markdown table cell."""
    return text.replace("|", "\\|").replace("\r\n", "\n").replace("\r", "\n").replace("\n", "<br>")


def _md_row(cells: list[str]) -> str:
    return "| " + " | ".join(cells) + " |"


def _entry_cells(entry: ResultEntry, ref_times: Mapping[str, float] | None) -> list[str]:
    ref = None if ref_times is None else ref_times.get(entry.tag)
    if entry.error:
        error_cell = _md_cell(entry.error)
    elif entry.skip_reason:
        error_cell = _span("yellow", _md_cell(entry.skip_reason))
    else:
        error_cell = ""
    return [
        _md_cell(entry.tag),
        _md_cell(entry.group),
        format_phase_html(entry.phases.get("init")),
        format_phase_html(entry.phases.get("build")),
        format_phase_html(entry.phases.get("run")),
        f"{entry.comp_time:.2f}",
        "-" if ref is None else f"{ref:.2f}",
        "-" if ref is None else format_time_diff(entry.comp_time, ref),
        format_tests_html(entry),
        format_tests_summary_html(entry),
        error_cell,
    ]


def write_markdown_log(
    entries: list[ResultEntry],
    path: Path,
    ref_times: Mapping[str, float] | None = None,
) -> Path:
    """Write the pinned GitHub-flavored markdown summary table; returns the path.

    ``ref_times`` maps case tag -> reference computational time (seconds);
    tags without a reference render ``-`` in the Ref. time / Time diff %
    columns (ported semantics from ``src/tools/testSolver.py``). Same
    parent-directory policy as :func:`write_json_log` (no silent mkdir).
    """
    path = _require_writable_parent(path, "markdown log")
    lines = [
        _md_row(list(MD_COLUMNS)),
        _md_row(["---"] * len(MD_COLUMNS)),
    ]
    for entry in entries:
        lines.append(_md_row(_entry_cells(entry, ref_times)))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _selftest() -> None:
    """Write a synthetic 2-entry session to a temp dir and verify everything."""
    import tempfile

    pass_entry = ResultEntry(
        tag="demo:pass-case",
        group="demo",
        tier="regression",
        phases={
            "init": PhaseResult(status="PASS", returncode=0, seconds=1.25, log="init ok"),
            "build": PhaseResult(status="PASS", returncode=0, seconds=42.5, log=""),
            "run": PhaseResult(status="PASS", returncode=0, seconds=12.5, log="run ok"),
        },
        tests_passed=4,
        tests_total=4,
        comp_time=12.5,
    )
    fail_entry = ResultEntry(
        tag="demo:fail-case",
        group="demo",
        tier="regression",
        phases={
            "init": PhaseResult(status="PASS", returncode=0, seconds=1.0, log=""),
            "build": PhaseResult(status="FAIL", returncode=1, seconds=3.0, log="cicc: error"),
            # no "run" phase -> "-" cell
        },
        tests_passed=0,
        tests_total=4,
        comp_time=0.0,
        error="build failed: cicc: error",
    )
    skip_entry = ResultEntry(
        tag="demo:Skip|case\nmultiline",
        group="demo",
        tier="regression",
        phases={"run": PhaseResult(status="SKIP", returncode=0, seconds=0.0, log="")},
        skip_reason="unfrozen reference",
    )
    entries = [pass_entry, fail_entry]

    # JSON round-trip of every field, nested PhaseResults included.
    for entry in entries + [skip_entry]:
        assert ResultEntry.from_json(entry.to_json()) == entry, entry

    machine = {"gpu": "selftest-gpu", "host": "selftest-host", "ignored": "x"}
    session = session_log(entries, machine, started=1_700_000_000.0)
    assert list(session.keys()) == ["schema", "started", "duration_s", "machine", "entries"]
    assert session["schema"] == 1
    assert isinstance(session["started"], str)
    datetime.fromisoformat(session["started"])  # parses as ISO-8601
    assert isinstance(session["duration_s"], float)
    assert session["machine"] == {"gpu": "selftest-gpu", "host": "selftest-host"}
    assert session["entries"] == [e.to_json() for e in entries]

    with tempfile.TemporaryDirectory() as tmp:
        json_path = write_json_log(session, Path(tmp) / "session.json")
        text = Path(json_path).read_text(encoding="utf-8")
        assert text == json.dumps(session, indent=2) + "\n"
        loaded = json.loads(text)
        assert loaded == session
        reloaded = [ResultEntry.from_json(e) for e in loaded["entries"]]
        assert reloaded == entries  # deep equality incl. nested PhaseResults

        ref_times = {pass_entry.tag: 10.0}
        md_path = write_markdown_log(entries, Path(tmp) / "log.md", ref_times=ref_times)
        md = Path(md_path).read_text(encoding="utf-8").splitlines()
        assert md[0] == (
            "| Case tag | Group | Init | Build | Run | Comp. time |"
            " Ref. time | Time diff % | Tests | Tests result | Error |"
        )
        assert md[1] == "| " + " | ".join(["---"] * len(MD_COLUMNS)) + " |"
        assert len(md) == 2 + len(entries), md

        green_pass = _span("green", "PASS")
        red_fail = _span("red", "FAIL")
        assert green_pass in md[2] and red_fail not in md[2]
        assert _span("green", "4/4") in md[2] and _span("green", "Passed") in md[2]
        # 12.5 vs ref 10.0 -> +25.0 % red span (testSolver rule).
        assert "10.00" in md[2] and _span("red", "25.0 %") in md[2]
        assert red_fail in md[3]
        # fail entry has no "run" phase -> "-" cell; no ref -> "-"; error text present.
        cells = [c.strip() for c in md[3].strip().strip("|").split("|")]
        assert cells[4] == "-" and cells[6] == "-" and cells[7] == "-"
        assert "build failed: cicc: error" in md[3]
        assert _span("red", "0/4") in md[3] and _span("red", "Failed") in md[3]

        # skip_reason renders yellow in the Error cell of its own row.
        skip_cells = _entry_cells(skip_entry, None)
        assert skip_cells[2] == "-" and skip_cells[4] == _span("yellow", "SKIP")
        assert skip_cells[10] == _span("yellow", "unfrozen reference")
        assert skip_cells[0] == "demo:Skip\\|case<br>multiline"  # GFM escaping

        # Failure QA: unwritable/missing parent -> RuntimeError naming the path.
        bad = Path(tmp) / "no-such-dir" / "session.json"
        try:
            write_json_log(session, bad)
        except RuntimeError as exc:
            assert str(bad) in str(exc), exc
        else:
            raise AssertionError("write_json_log must raise RuntimeError on missing parent")
    print("results selftest PASS")


if __name__ == "__main__":
    _selftest()
