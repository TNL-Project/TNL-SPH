"""Informational per-GPU computational-time comparison for the test suite.

Rich-side port of the referential-time comparison from
``src/tools/testSolver.py:128-196`` (same logic in ``src/tools/runAll.py:89-111``)
into the unified test suite -- the OLD pandas/span table columns are already
covered by :func:`suite.results.write_markdown_log`, so this module provides
ONLY the rich/perf side, consumed through the ``perf_fn`` hook of
:func:`suite.reporting.print_final_summary` (plan todo 11; console
integration owned by todo 28).

Semantics (pinned -- INFO-ONLY, NEVER a failure signal of the suite):

* reference times live in
  ``src/Tests/regression/referentialComputationalTimes/referential_comp_times_<gpu>.json``
  keyed by ``case-tag`` strings (e.g. ``WCSPH-DBC/damBreak2D_WCSPH-DBC`` of
  ``configurationsToTest.py``); variant tags without a reference
  (``:MGVT`` etc.) simply render ``-`` cells,
* diff % = ``100 * (t - ref) / ref`` with one decimal, ``[red]`` markup when
  ``> 5``, ``[green]`` when ``< -5``, plain otherwise -- thresholds ported
  from ``testSolver.py:155-160``,
* the GPU name comes from ``nvidia-smi --query-gpu=name --format=csv,noheader``
  (first row, spaces->dashes) or verbatim from an override (the old ``--gpu``
  flag behavior, ``testSolver.py:202-204`` TODO resolved); ANY lookup failure
  degrades to a single skip message, never an exception.

Selftest (synthetic entries + the real A40 reference file, writes nothing,
NEVER exits non-zero)::

    PYTHONPATH=src/Tests python3 -m regression.perf
    PYTHONPATH=src/Tests python3 -m regression.perf --selftest
"""

from __future__ import annotations

import json
import subprocess
from typing import Iterable, Mapping

from suite import paths
from suite.results import ResultEntry

# Diff % thresholds ported from src/tools/testSolver.py:155-160.
_RED_THRESHOLD = 5.0
_GREEN_THRESHOLD = -5.0

_NO_REF_MESSAGE = "perf: no reference times for '{gpu}', perf report skipped"


def detect_gpu(override: str | None = None) -> str:
    """Return the GPU name with spaces -> dashes (``NVIDIA A40`` -> ``NVIDIA-A40``).

    ``override`` (the old ``--gpu`` flag) is returned verbatim; without it the
    first row of ``nvidia-smi --query-gpu=name --format=csv,noheader`` is used.
    Any failure (no driver, no ``nvidia-smi``, empty output) returns ``""`` --
    callers render it as ``-``.
    """
    if override is not None:
        return override
    try:
        proc = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            check=True,
            timeout=15,
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    rows = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
    if not rows:
        return ""
    return rows[0].replace(" ", "-")


def load_reference_times(gpu_name: str) -> dict[str, float] | None:
    """Return ``{case_tag: reference_seconds}`` for ``gpu_name``, else ``None``.

    ANY failure (missing file, unknown GPU, unparsable content) returns
    ``None`` and the caller prints the skip message exactly once -- a missing
    GPU reference file is informational only (pinned by plan todo 11).
    """
    if not gpu_name:
        return None
    ref_path = (
        paths.TESTS_DIR
        / "regression"
        / "referentialComputationalTimes"
        / f"referential_comp_times_{gpu_name}.json"
    )
    try:
        raw = json.loads(ref_path.read_text(encoding="utf-8"))
        return {str(tag): float(value) for tag, value in raw.items()}
    except (OSError, ValueError, TypeError, AttributeError):
        return None


def perf_columns(
    entry: ResultEntry, ref_times: Mapping[str, float] | None
) -> tuple[str, str, str]:
    """Return the rich-markup ``(time, ref_time, diff%)`` cells for ``entry``.

    Rules (port of ``testSolver.py:149-162``): ``t <= 0`` (upstream parse
    failure) renders the time and diff cells as ``-``; a tag without a
    reference renders reference and diff as ``-``; otherwise diff % is
    ``[red]`` above 5, ``[green]`` below -5, plain otherwise.
    """
    t = float(entry.comp_time or 0.0)
    ref = None if ref_times is None else ref_times.get(entry.tag)
    time_str = f"{t:.2f}" if t > 0 else "-"
    if ref is None:
        return time_str, "-", "-"
    ref_str = f"{ref:.2f}"
    if t <= 0 or ref <= 0:
        # diff undefined without a positive measured/reference pair.
        return time_str, ref_str, "-"
    diff = 100.0 * (t - ref) / ref
    diff_str = f"{diff:.1f} %"
    if diff > _RED_THRESHOLD:
        diff_str = f"[red]{diff_str}[/red]"
    elif diff < _GREEN_THRESHOLD:
        diff_str = f"[green]{diff_str}[/green]"
    return time_str, ref_str, diff_str


def attach_perf(
    entries: Iterable[ResultEntry], gpu_name: str | None = None
) -> dict[str, tuple[str, str, str]]:
    """Map ``entry.tag -> perf_columns`` for every entry (NEVER raises).

    Resolves the GPU (``gpu_name`` override verbatim, else ``nvidia-smi``
    autodetect), loads its reference times, and -- when the GPU has none --
    prints the skip message exactly once and degrades all columns to
    info-only ``-`` cells.
    """
    gpu = detect_gpu(gpu_name)
    ref_times = load_reference_times(gpu)
    if ref_times is None:
        print(_NO_REF_MESSAGE.format(gpu=gpu))
    return {entry.tag: perf_columns(entry, ref_times) for entry in entries}


def _selftest() -> None:
    """Exercise every branch with synthetic entries; NEVER exits non-zero."""
    import contextlib
    import io

    from rich.console import Console

    from suite import reporting
    from suite.results import PhaseResult

    # Real A40 reference file (keys are case-tag strings, cross-checked
    # against configurationsToTest.py).
    refs = load_reference_times("NVIDIA-A40")
    assert refs is not None and refs["WCSPH-DBC/damBreak2D_WCSPH-DBC"] == 10.398
    assert "RSPH/damBreak2D_RSPH" in refs

    # Unknown GPU / empty autodetect result degrade to None.
    assert load_reference_times("NVIDIA-BOGUS-XYZ") is None
    assert load_reference_times("") is None

    # detect_gpu: override returns verbatim; autodetect never raises.
    assert detect_gpu("NVIDIA-BOGUS-XYZ") == "NVIDIA-BOGUS-XYZ"
    detected = detect_gpu()
    assert isinstance(detected, str)
    print(f"perf selftest: autodetected gpu = {detected or '-'}")

    def entry(tag: str, comp_time: float) -> ResultEntry:
        return ResultEntry(
            tag=tag,
            group="demo",
            tier="regression",
            phases={"run": PhaseResult(status="PASS", returncode=0, seconds=comp_time)},
            tests_passed=1,
            tests_total=1,
            comp_time=comp_time,
        )

    # perf_columns branches (rich markup, one decimal, testSolver.py:155-160).
    red = perf_columns(entry("WCSPH-DBC/damBreak2D_WCSPH-DBC", 11.5), refs)
    assert red == ("11.50", "10.40", "[red]10.6 %[/red]"), red
    green = perf_columns(entry("RSPH/damBreak2D_RSPH", 49.0), refs)
    assert green == ("49.00", "52.85", "[green]-7.3 %[/green]"), green
    neutral = perf_columns(entry("WCSPH-BI/damBreak2D_WCSPH-BI", 57.0), refs)
    assert neutral == ("57.00", "56.41", "1.0 %"), neutral
    missing = perf_columns(entry("demo:no-such-tag", 5.0), refs)
    assert missing == ("5.00", "-", "-"), missing
    no_time = perf_columns(entry("WCSPH-DBC/damBreak2D_WCSPH-DBC", 0.0), refs)
    assert no_time == ("-", "10.40", "-"), no_time
    zero_ref = perf_columns(
        entry("WCSPH-BI/poiseuilleFlowWithOpenBoundary2D_WCSPH-BI", 12.0), refs
    )
    assert zero_ref == ("12.00", "0.00", "-"), zero_ref
    none_refs = perf_columns(entry("WCSPH-DBC/damBreak2D_WCSPH-DBC", 11.5), None)
    assert none_refs == ("11.50", "-", "-"), none_refs

    demo_entries = [
        entry("WCSPH-DBC/damBreak2D_WCSPH-DBC", 11.5),  # >+5% -> red row
        entry("RSPH/damBreak2D_RSPH", 49.0),  # <-5% -> green row
        entry("WCSPH-BI/damBreak2D_WCSPH-BI", 57.0),  # neutral row
        entry("demo:no-such-tag", 5.0),  # missing ref -> "-" cells
    ]

    # Bogus --gpu XYZ path: one explicit skip message, all info-only dashes.
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        skipped = attach_perf(demo_entries, "XYZ-BOGUS")
    print(buf.getvalue(), end="")
    assert "perf: no reference times for 'XYZ-BOGUS', perf report skipped" in buf.getvalue()
    assert buf.getvalue().count("perf report skipped") == 1
    assert all(cols[1] == "-" and cols[2] == "-" for cols in skipped.values())

    # Autodetect attach path never raises and covers every entry.
    attached = attach_perf(demo_entries)
    assert set(attached) == {e.tag for e in demo_entries}

    # Integration: 3 info-only perf columns through print_final_summary.
    perf_map = attach_perf(demo_entries, "NVIDIA-A40")

    def perf_fn(item: ResultEntry) -> tuple[str, str, str]:
        return perf_map[item.tag]

    original_console = reporting._CONSOLE
    reporting._CONSOLE = Console(record=True, width=200)
    try:
        rc = reporting.print_final_summary({"regression": demo_entries}, perf_fn=perf_fn)
        rendered = reporting._CONSOLE.export_text(styles=False)
    finally:
        reporting._CONSOLE = original_console

    assert rc == 0  # info-only perf never changes the verdict (all PASS here)
    assert "Ref. time" in rendered and "Time diff %" in rendered
    assert "10.6 %" in rendered and "-7.3 %" in rendered and "1.0 %" in rendered
    assert "demo:no-such-tag" in rendered
    assert "info-only perf report" in rendered

    print("perf selftest PASS")


if __name__ == "__main__":
    _selftest()
