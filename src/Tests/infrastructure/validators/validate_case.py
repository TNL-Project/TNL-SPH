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


def validate(case_rel: str,
             config_name: str | None = None,
             distributed: int | None = None) -> list[Finding]:
    """Programmatic entry point for the tier runners (infra boot driver).

    Runs the same expectations/invariant pipeline as :func:`main` but
    without any printing and without the ``--ref`` deep-compare; returns the
    collected findings (an empty list means all rules passed).  The CLI in
    :func:`main` below is untouched and stays byte-compatible.
    """
    case_dir = _resolve_case_dir(case_rel)
    config_path = _resolve_config(case_dir, config_name)
    # A case's distributed subdomain layout comes from the nested subdomains
    # loader and is not always mirrored by a top-level numberOfSubdomains
    # key (e.g. damBreak3D_*_distributed inlines it), so the explicit
    # distributed=N flag is the n_subdomains authority for the rule gates.
    case_expectations = exp_mod.from_config(
        config_path,
        init_overrides={"number_of_subdomains": int(distributed)}
        if distributed is not None else None)

    results_dir = case_dir / str(case_expectations.get("output_directory")
                                 or "results")
    # Under MPI every rank writes init_report_subdomain<i>.json and NO
    # plain init_report.json exists (the emission is per-rank only) -- a
    # distributed validation therefore anchors the shared rules at rank 0
    # and lets the tiling rule cover the rest.
    report_path = results_dir / "init_report.json"
    if distributed is not None and not report_path.exists():
        report_path = results_dir / "init_report_subdomain0.json"
    report = _load_report(report_path)

    case_expectations["vtk_counts"] = exp_mod.vtk_counts(
        case_dir / "sources",
        fluid_vtk=case_expectations.get("fluid_vtk"),
        boundary_vtk=case_expectations.get("boundary_vtk"))

    findings: list[Finding] = []
    if distributed is not None:
        rank_reports = []
        for rank in range(distributed):
            rank_path = results_dir / f"init_report_subdomain{rank}.json"
            if not rank_path.exists():
                findings.append(Finding(
                    DECOMP_TILING,
                    f"rank file missing: {rank_path} "
                    f"(expected {distributed} per-rank reports)"))
                continue
            rank_reports.append(_load_report(rank_path))
        case_expectations["rank_reports"] = rank_reports
        case_expectations["distributed"] = True

    for rule in RULES:
        name = RULE_NAMES[rule]
        if name == DECOMP_TILING and distributed is None:
            continue
        finding = rule(report, case_expectations)
        if finding is not None:
            findings.append(finding)
    return findings


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
    expectations = exp_mod.from_config(
        config_path,
        init_overrides={"number_of_subdomains": int(args.distributed)}
        if args.distributed is not None else None)

    results_dir = case_dir / str(expectations.get("output_directory") or "results")
    report_path = results_dir / "init_report.json"
    if args.distributed is not None and not report_path.exists():
        report_path = results_dir / "init_report_subdomain0.json"
    report = _load_report(report_path)

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
    print(f"report : {report_path}")
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
