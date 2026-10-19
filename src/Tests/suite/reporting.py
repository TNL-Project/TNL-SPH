"""Console reporting of the TNL-SPH unified test suite (rich).

This is the ONLY ``suite`` module allowed to print. All state semantics are
shared with :mod:`suite.results`:

* an entry FAILS when :func:`suite.results.result_failed` is true
  (error text or any FAIL phase),
* an entry is SKIPPED when :func:`suite.results.result_skipped` is true
  (``skip_reason`` set) -- skips are NEVER conflated with failures,
* everything else PASSES.

``rich`` is a HARD dependency of this module, but it is imported lazily
inside :func:`_console` so that importing the module (and all of
:mod:`suite.results`) stays possible without it; a missing install is
announced with an install hint instead of a bare ImportError.

Selftest (prints a PASS + FAIL + SKIP table, writes no files)::

    PYTHONPATH=src/Tests python3 -m suite.reporting --selftest
"""

from __future__ import annotations

from typing import Iterable, Mapping

from . import results
from .results import ResultEntry

_CONSOLE = None


def _console():
    """Return the shared rich Console, announcing a missing rich with a hint."""
    global _CONSOLE
    if _CONSOLE is None:
        try:
            from rich.console import Console
        except ImportError as exc:
            raise SystemExit(
                "suite.reporting requires the 'rich' package (hard dependency).\n"
                "Install it with: pip install rich"
            ) from exc
        _CONSOLE = Console()
    return _CONSOLE


def _escape(text: str) -> str:
    from rich.markup import escape

    return escape(text)


def print_phase(message: str, ok: bool = True, error_detail: str | None = None) -> None:
    """Print one colored phase line; failed phases may show an indented detail."""
    console = _console()
    if ok:
        console.print(f"[bold green]  ✓  {_escape(message)}[/bold green]")
    else:
        console.print(f"[bold red]  ✗  {_escape(message)}[/bold red]")
        if error_detail:
            console.print(f"      [red]{_escape(error_detail)}[/red]")


def _result_cell(entry: ResultEntry) -> str:
    if results.result_failed(entry):
        return "[red]Failed[/red]"
    if results.result_skipped(entry):
        return "[yellow]Skipped[/yellow]"
    return "[green]Success[/green]"


def _tests_cells(entry: ResultEntry) -> tuple[str, str]:
    """Rich-markup port of src/tools/multirunner.py:194-204."""
    passed = entry.tests_passed
    total = entry.tests_total
    if total == 0:
        return "-", "-"
    if passed == total:
        return f"[green]{passed}/{total}[/green]", "[green]Passed[/green]"
    return f"[red]{passed}/{total}[/red]", "[red]Failed[/red]"


def summary_table(entries: Iterable[ResultEntry], title: str = "Test results") -> None:
    """Print the rich summary table (ported from multirunner.print_summary_table)."""
    console = _console()
    from rich.table import Table

    table = Table(title=title, show_lines=True)
    table.add_column("Case tag")
    table.add_column("Result")
    table.add_column("Comp. time")
    table.add_column("Tests")
    table.add_column("Tests result")
    table.add_column("Error")

    for entry in entries:
        tests_str, tests_result_str = _tests_cells(entry)
        if entry.error:
            err = entry.error
        else:
            err = entry.skip_reason
        err_str = _escape(err[:60]) if err else ""
        table.add_row(
            entry.tag,
            _result_cell(entry),
            f"{entry.comp_time:.2f}",
            tests_str,
            tests_result_str,
            err_str,
        )

    console.print(table)


def print_final_summary(entries_by_tier: Mapping[str, Iterable[ResultEntry]]) -> int:
    """Print aggregate pass/fail/skip counts per tier and overall.

    Returns the suggested process exit code: ``1`` when at least one entry
    failed, else ``0``. Skipped entries are counted separately and never
    conflated with failures.
    """
    console = _console()
    total_pass = total_fail = total_skip = 0
    for tier, tier_entries in entries_by_tier.items():
        tier_entries = list(tier_entries)
        n_skip = sum(1 for e in tier_entries if results.result_skipped(e))
        n_fail = sum(
            1
            for e in tier_entries
            if results.result_failed(e) and not results.result_skipped(e)
        )
        n_pass = len(tier_entries) - n_fail - n_skip
        total_pass += n_pass
        total_fail += n_fail
        total_skip += n_skip
        console.print(
            f"{tier}: [green]{n_pass} passed[/green], "
            f"[red]{n_fail} failed[/red], "
            f"[yellow]{n_skip} skipped[/yellow]"
        )
    console.print(
        f"total: [green]{total_pass} passed[/green], "
        f"[red]{total_fail} failed[/red], "
        f"[yellow]{total_skip} skipped[/yellow]"
    )
    if total_fail:
        console.print("[bold red]Some tests FAILED -- exit code must be non-zero.[/bold red]")
        return 1
    return 0


def _selftest() -> None:
    from .results import PhaseResult

    pass_entry = ResultEntry(
        tag="demo:pass-case",
        group="demo",
        tier="demo",
        phases={
            "init": PhaseResult(status="PASS", returncode=0, seconds=1.0),
            "build": PhaseResult(status="PASS", returncode=0, seconds=40.0),
            "run": PhaseResult(status="PASS", returncode=0, seconds=12.5),
        },
        tests_passed=4,
        tests_total=4,
        comp_time=12.5,
    )
    fail_entry = ResultEntry(
        tag="demo:fail-case",
        group="demo",
        tier="demo",
        phases={
            "init": PhaseResult(status="PASS", returncode=0, seconds=1.0),
            "build": PhaseResult(status="FAIL", returncode=1, seconds=3.0, log="cicc: error"),
        },
        tests_passed=0,
        tests_total=4,
        comp_time=0.0,
        error="build failed: cicc: error",
    )
    skip_entry = ResultEntry(
        tag="demo:skip-case",
        group="demo",
        tier="demo",
        phases={"run": PhaseResult(status="SKIP", returncode=0, seconds=0.0)},
        skip_reason="unfrozen reference",
    )
    entries = [pass_entry, fail_entry, skip_entry]

    print_phase("init demo/pass-case", ok=True)
    print_phase("build demo/fail-case", ok=False, error_detail="cicc: error")
    summary_table(entries, title="reporting selftest")

    rc = print_final_summary({"demo": entries})
    assert rc == 1  # 1 passed, 1 failed, 1 skipped -> failure -> hint + 1
    rc_ok = print_final_summary({"demo": [pass_entry, skip_entry]})
    assert rc_ok == 0  # skips alone must not fail the summary
    print("reporting selftest PASS")


if __name__ == "__main__":
    _selftest()
