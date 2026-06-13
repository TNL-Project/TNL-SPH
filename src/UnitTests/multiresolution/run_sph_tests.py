#! /usr/bin/env python3
"""
run_sph_tests.py

Compares parsed SPH log files against reference JSON files.
Can also orchestrate initialization and running of test configurations.

Test cases are defined in a YAML or JSON test-list file (see --tests).
Each entry specifies a log file and its reference JSON.

Usage:
    # Generate reference JSON from a log:
    python run_sph_tests.py --generate <logfile> [--ref <ref.json>]

    # Run all tests defined in tests.yaml:
    python run_sph_tests.py --tests tests.yaml

    # Run a single comparison directly:
    python run_sph_tests.py --log <logfile> --ref <ref.json>

    # Initialize test configurations:
    python run_sph_tests.py --init [--config-name <name>]

    # Run test configurations:
    python run_sph_tests.py --run [--config-name <name>]

    # Initialize and run test configurations:
    python run_sph_tests.py --init-and-run [--config-name <name>]

    # Full pipeline: initialize, run, then compare results:
    python run_sph_tests.py --all
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional

# ── Fancy outout ───────────────────────────────────────────────────────────────
BOLD_GREEN = '\033[1;32m'
BOLD_RED   = '\033[1;31m'
RESET      = '\033[0m'


def pass_msg(msg='PASS'):
    print(f'{BOLD_GREEN}  ✓  {msg}{RESET}')


def fail_msg(msg):
    print(f'{BOLD_RED}  ✗  {msg}{RESET}')

# ── Try yaml, fall back to json-only test lists ───────────────────────────────
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

from parse_sph_log import parse_log

# ── Test configuration directory ──────────────────────────────────────────────
_test_dir = Path(__file__).parent / "testConfigurations" / "dummyMultiresolutionSimulation2D"
_default_tests_yaml = Path(__file__).parent / "tests.yaml"

sys.path.insert(0, str(_test_dir))
from configurations import CONFIGURATIONS
CONFIG_NAMES = list(CONFIGURATIONS.keys())


# ─────────────────────────────────────────────────────────────────────────────
# Deep comparison
# ─────────────────────────────────────────────────────────────────────────────

class Diff:
    def __init__(self, path: str, expected: Any, actual: Any):
        self.path = path
        self.expected = expected
        self.actual = actual

    def __str__(self):
        return (
            f"  {self.path}\n"
            f"    expected : {self.expected!r}\n"
            f"    actual   : {self.actual!r}"
        )


def deep_compare(expected: Any, actual: Any, path: str = '', tol: float = 1e-9) -> list[Diff]:
    """Recursively compare two structures. Returns list of Diff objects."""
    diffs = []

    if isinstance(expected, dict):
        if not isinstance(actual, dict):
            return [Diff(path, expected, actual)]
        all_keys = set(expected) | set(actual)
        for k in sorted(all_keys):
            child = f'{path}.{k}' if path else k
            if k not in expected:
                #diffs.append(Diff(child, '<missing>', actual[k]))
                continue
            elif k not in actual:
                diffs.append(Diff(child, expected[k], '<missing>'))
            else:
                diffs.extend(deep_compare(expected[k], actual[k], child, tol))

    elif isinstance(expected, list):
        if not isinstance(actual, list):
            return [Diff(path, expected, actual)]
        if len(expected) != len(actual):
            diffs.append(Diff(path + '.length', len(expected), len(actual)))
        for i, (e, a) in enumerate(zip(expected, actual)):
            diffs.extend(deep_compare(e, a, f'{path}[{i}]', tol))

    elif isinstance(expected, float) or isinstance(actual, float):
        try:
            e, a = float(expected), float(actual)
            if abs(e - a) > tol * max(1.0, abs(e)):
                diffs.append(Diff(path, expected, actual))
        except (TypeError, ValueError):
            if expected != actual:
                diffs.append(Diff(path, expected, actual))

    else:
        if expected != actual:
            diffs.append(Diff(path, expected, actual))

    return diffs


# ─────────────────────────────────────────────────────────────────────────────
# Single test
# ─────────────────────────────────────────────────────────────────────────────

def run_single_test(log_path: Path, ref_path: Path, tol: float = 1e-9,
                    name: str = '') -> bool:
    label = name or log_path.name
    print(f'\n{"─"*60}')
    print(f'TEST: {label}')
    print(f'  log : {log_path}')
    print(f'  ref : {ref_path}')

    # Validate paths
    if not log_path.exists():
        #print(f'  ✗  Log file not found: {log_path}')
        fail_msg(f'Log file not found: {log_path}')
        return False
    if not ref_path.exists():
        #print(f'  ✗  Reference JSON not found: {ref_path}')
        fail_msg(f'Reference JSON not found: {ref_path}')
        return False

    # Parse log
    try:
        actual = parse_log(log_path.read_text(encoding='utf-8'))
    except Exception as e:
        #print(f'  ✗  Parse error: {e}')
        fail_msg(f'Parse error: {e}')
        return False

    # Load reference
    try:
        expected = json.loads(ref_path.read_text(encoding='utf-8'))
    except Exception as e:
        #print(f'  ✗  Could not load reference JSON: {e}')
        fail_msg(f'Could not load reference JSON: {e}')
        return False

    diffs = deep_compare(expected, actual, tol=tol)
    if not diffs:
        #print('  ✓  PASS')
        pass_msg()
        return True
    else:
        fail_msg(f'FAIL  ({len(diffs)} difference(s)):')
        for d in diffs:
            print(d)
        return False


# ─────────────────────────────────────────────────────────────────────────────
# Load test list
# ─────────────────────────────────────────────────────────────────────────────

def load_test_list(path: Path) -> list[dict]:
    """
    Load a YAML or JSON test list.

    YAML format:
        tolerance: 1e-9          # optional global tolerance

        tests:
          - name: damBreak2D     # optional label
            log:  results/sim.log
            ref:  refs/sim.json
          - name: anotherCase
            log:  results/other.log
            ref:  refs/other.json

    JSON format (same structure):
        {
          "tolerance": 1e-9,
          "tests": [
            { "name": "damBreak2D", "log": "results/sim.log", "ref": "refs/sim.json" }
          ]
        }
    """
    text = path.read_text(encoding='utf-8')
    if path.suffix in ('.yaml', '.yml'):
        if not HAS_YAML:
            sys.exit('PyYAML not installed. Run: pip install pyyaml  — or use a .json test list.')
        data = yaml.safe_load(text)
    else:
        data = json.loads(text)

    if isinstance(data, list):
        # bare list — no global tolerance
        return 1e-9, data
    tol = data.get('tolerance', 1e-9)
    return tol, data.get('tests', [])


# ─────────────────────────────────────────────────────────────────────────────
# Generate reference JSON
# ─────────────────────────────────────────────────────────────────────────────

def generate_reference(log_path: Path, ref_path: Path):
    text = log_path.read_text(encoding='utf-8')
    data = parse_log(text)
    ref_path.parent.mkdir(parents=True, exist_ok=True)
    ref_path.write_text(json.dumps(data, indent=2), encoding='utf-8')
    print(f'Reference JSON written → {ref_path}')
    print(json.dumps(data, indent=2))


# ─────────────────────────────────────────────────────────────────────────────
# Batch test comparison
# ─────────────────────────────────────────────────────────────────────────────

def run_batch_tests(list_path: Path, cli_tol: float = 1e-9) -> int:
    """Run batch test comparison. Returns number of failures."""
    if not list_path.exists():
        sys.exit(f'Test list not found: {list_path}')

    global_tol, cases = load_test_list(list_path)
    tol = cli_tol if cli_tol != 1e-9 else global_tol

    if not cases:
        print('No test cases found in list.')
        return 0

    base_dir = list_path.parent
    passed = failed = 0

    for case in cases:
        log_p = base_dir / case['log']
        ref_p = base_dir / case['ref']
        name  = case.get('name', '')
        case_tol = case.get('tolerance', tol)

        ok = run_single_test(log_p, ref_p, tol=case_tol, name=name)
        if ok:
            passed += 1
        else:
            failed += 1

    print(f'\n{"═"*60}')
    print(f'Results: {passed} passed, {failed} failed  (total {passed+failed})')
    return failed


# ─────────────────────────────────────────────────────────────────────────────
# Init / Run helpers
# ─────────────────────────────────────────────────────────────────────────────

def _subprocess_run(cmd: list, cwd: Path, label: str):
    print(f"\n{'═'*60}")
    print(f"  {label}")
    print(f"{'═'*60}\n")
    result = subprocess.run(cmd, cwd=str(cwd))
    if result.returncode != 0:
        sys.exit(f"{label} failed with return code {result.returncode}")


def do_init(config_name: Optional[str] = None):
    script = _test_dir / "init_2dmr_configuration.py"
    cmd = [sys.executable, str(script)]
    if config_name:
        cmd += ["--config-name", config_name]
    else:
        cmd += ["--all"]
    target = config_name or "all configurations"
    _subprocess_run(cmd, _test_dir, f"Initializing {target}")


def do_run(config_name: Optional[str] = None):
    script = _test_dir / "run.py"
    cmd = [sys.executable, str(script)]
    if config_name:
        cmd += ["--config-name", config_name]
    else:
        cmd += ["--all"]
    target = config_name or "all configurations"
    _subprocess_run(cmd, _test_dir, f"Running {target}")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='SPH multiresolution test runner — init, run, and compare',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--generate', metavar='LOG',
                      help='Parse LOG and write reference JSON (use --ref for output path)')
    mode.add_argument('--tests', metavar='LIST',
                      help='YAML/JSON file listing test cases to run')
    mode.add_argument('--log', metavar='LOG',
                      help='Single log to compare (requires --ref)')
    mode.add_argument('--init', action='store_true',
                      help='Initialize test configurations')
    mode.add_argument('--run', action='store_true',
                      help='Run test configurations')
    mode.add_argument('--init-and-run', action='store_true',
                      help='Initialize and run test configurations')
    mode.add_argument('--all', action='store_true',
                      help='Initialize, run, and compare all test configurations')

    ap.add_argument('--config-name', choices=CONFIG_NAMES,
                    help='Target a single configuration (applies to --init, --run, --init-and-run)')
    ap.add_argument('--ref', metavar='JSON',
                    help='Reference JSON path (used with --generate or --log)')
    ap.add_argument('--tol', type=float, default=1e-9,
                    help='Floating-point tolerance for comparisons (default: 1e-9)')
    args = ap.parse_args()

    # ── Generate mode ─────────────────────────────────────────────────────
    if args.generate:
        log_path = Path(args.generate)
        ref_path = Path(args.ref) if args.ref else log_path.with_suffix('.ref.json')
        generate_reference(log_path, ref_path)
        return

    # ── Single comparison ─────────────────────────────────────────────────
    if args.log:
        if not args.ref:
            ap.error('--log requires --ref')
        ok = run_single_test(Path(args.log), Path(args.ref), tol=args.tol)
        sys.exit(0 if ok else 1)

    # ── Init mode ─────────────────────────────────────────────────────────
    if args.init:
        do_init(args.config_name)
        return

    # ── Run mode ──────────────────────────────────────────────────────────
    if args.run:
        do_run(args.config_name)
        return

    # ── Init-and-run mode ─────────────────────────────────────────────────
    if args.init_and_run:
        do_init(args.config_name)
        do_run(args.config_name)
        return

    # ── All mode: init + run + compare ────────────────────────────────────
    if args.all:
        do_init()
        do_run()
        failed = run_batch_tests(_default_tests_yaml, cli_tol=args.tol)
        sys.exit(0 if failed == 0 else 1)

    # ── Batch from test list ──────────────────────────────────────────────
    list_path = Path(args.tests)
    failed = run_batch_tests(list_path, cli_tol=args.tol)
    sys.exit(0 if failed == 0 else 1)


if __name__ == '__main__':
    main()
