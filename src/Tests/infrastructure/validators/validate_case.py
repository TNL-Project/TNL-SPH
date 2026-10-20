"""
validate_case.py

Validate a case's ``results/init_report.json`` against the expectations
derived from its own configuration artifacts.

Usage
-----
    PYTHONPATH=src/Tests python3 -m infrastructure.validators.validate_case \
        --case WCSPH-DBC/damBreak2D_WCSPH-DBC [--config-name name] \
        [--distributed N] [--ref <path-or-case-tag>]

``--case`` is resolved relative to the repository's ``examples/`` first,
then relative to CWD.  ``--config-name name`` selects a non-default config:
``sources/<name>.jsonc`` if present, else ``sources/<name>/config.jsonc``,
else the single non-auxiliary ``*.jsonc`` under ``sources/<name>/``.
Exit code is 0 iff no FAIL finding was produced.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import compare, expectations as exp_mod
from .invariants import (
    DECOMP_TILING,
    Finding,
    RULE_NAMES,
    RULES,
)

_REPO_ROOT = Path(__file__).resolve().parents[4]

BOLD_GREEN = "\033[1;32m"
BOLD_RED = "\033[1;31m"
DIM = "\033[2m"
RESET = "\033[0m"

_AUX_CONFIG_PREFIXES = ("config-distributed-domain", "config-open-boundary",
                        "config-periodic-boundary", "config-measuretool")


def _resolve_case_dir(case_arg: str) -> Path:
    candidate = Path(case_arg)
    if candidate.is_absolute() and candidate.is_dir():
        return candidate.resolve()
    for base in (_REPO_ROOT / "examples",
                 Path.cwd(),
                 _REPO_ROOT / "src" / "Tests" / "infrastructure" / "fixtures"):
        resolved = base / case_arg
        if resolved.is_dir():
            return resolved.resolve()
    raise SystemExit(f"validate_case: case directory not found: {case_arg}")


def _resolve_config(case_dir: Path, config_name: str | None) -> Path:
    sources = case_dir / "sources"
    if config_name is None:
        default = sources / "config.jsonc"
        if default.exists():
            return default
        raise SystemExit(f"validate_case: default config not found: {default}")
    candidates = [
        sources / f"{config_name}.jsonc",
        sources / config_name / "config.jsonc",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    subdir = sources / config_name
    if subdir.is_dir():
        concrete = sorted(
            p for p in subdir.glob("*.jsonc")
            if not p.name.startswith(_AUX_CONFIG_PREFIXES))
        if len(concrete) == 1:
            return concrete[0]
    raise SystemExit(
        f"validate_case: config '{config_name}' not resolvable under {sources}")


def _load_report(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise SystemExit(f"validate_case: report not found: {path}"
                         f" (run the case binary past init first)") from None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="infrastructure.validators.validate_case",
        description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--case", required=True,
                        help="case path relative to examples/ (or literal)")
    parser.add_argument("--config-name", default=None,
                        help="non-default config name under sources/")
    parser.add_argument("--distributed", type=int, default=None, metavar="N",
                        help="read the N per-rank init_report_subdomain<i>.json "
                             "files and run the decomposition tiling rule")
    parser.add_argument("--ref", default=None, metavar="PATH_OR_TAG",
                        help="deep-compare the report against a frozen "
                             "reference JSON (path or references/ case tag)")
    args = parser.parse_args(argv)

    case_dir = _resolve_case_dir(args.case)
    config_path = _resolve_config(case_dir, args.config_name)
    expectations = exp_mod.from_config(config_path)

    results_dir = case_dir / str(expectations.get("output_directory") or "results")
    report = _load_report(results_dir / "init_report.json")

    expectations["vtk_counts"] = exp_mod.vtk_counts(
        case_dir / "sources",
        fluid_vtk=expectations.get("fluid_vtk"),
        boundary_vtk=expectations.get("boundary_vtk"))

    rank_findings: list[Finding] = []
    if args.distributed is not None:
        rank_reports = []
        for rank in range(args.distributed):
            rank_path = results_dir / f"init_report_subdomain{rank}.json"
            if not rank_path.exists():
                rank_findings.append(Finding(
                    DECOMP_TILING,
                    f"rank file missing: {rank_path} "
                    f"(expected {args.distributed} per-rank reports)"))
                continue
            rank_reports.append(_load_report(rank_path))
        expectations["rank_reports"] = rank_reports
        expectations["distributed"] = True

    print(f"case   : {case_dir}")
    print(f"config : {config_path}")
    print(f"report : {results_dir / 'init_report.json'}")
    vtk_counts = expectations["vtk_counts"]
    if vtk_counts.get("fluid") is not None or vtk_counts.get("error"):
        print(f"vtk    : fluid={vtk_counts.get('fluid')} "
              f"boundary={vtk_counts.get('boundary')}"
              + (f"  ({vtk_counts['error']})" if vtk_counts.get("error") else ""))

    findings: list[Finding] = []
    for rule in RULES:
        name = RULE_NAMES[rule]
        if name == DECOMP_TILING and args.distributed is None:
            print(f"{DIM}  -  {DECOMP_TILING} (skipped: --distributed N not given){RESET}")
            continue
        finding = rule(report, expectations)
        if finding is None:
            print(f"{BOLD_GREEN}  ✓  {name}{RESET}")
        else:
            findings.append(finding)
            print(f"{BOLD_RED}  ✗  {name}{RESET}")
            print(f"{BOLD_RED}       {finding.message}{RESET}")
    findings.extend(rank_findings)

    if args.ref:
        ref_path = Path(args.ref)
        if ref_path.exists():
            reference = json.loads(ref_path.read_text(encoding="utf-8"))
            source = str(ref_path)
        else:
            reference = compare.load_frozen_reference(args.ref)
            source = f"references/{args.ref}"
        if reference is None:
            print(f"{DIM}  -  reference.deep_compare (skipped: no frozen "
                  f"reference for '{args.ref}'){RESET}")
        else:
            diffs = compare.deep_compare(reference, report)
            if not diffs:
                print(f"{BOLD_GREEN}  ✓  reference.deep_compare{RESET}")
            else:
                findings.append(Finding(
                    "reference.deep_compare",
                    f"{len(diffs)} difference(s) vs {source}"))
                print(f"{BOLD_RED}  ✗  reference.deep_compare "
                      f"({len(diffs)} difference(s) vs {source}){RESET}")
                for diff in diffs:
                    print(diff)

    fail_count = sum(1 for f in findings if f.severity == "FAIL")
    print(f"\nfindings: {fail_count} FAIL / {len(findings)} total")
    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
