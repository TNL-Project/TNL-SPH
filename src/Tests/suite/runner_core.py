"""Runner primitives of the TNL-SPH unified test suite (init/run/restore/eval).

Ports the per-case mechanics of ``src/tools/multirunner.py:35-92`` onto the
suite's configuration contract (``suite.config_loader``):

- :func:`init_case` invokes the case's init script DIRECTLY
  (``conf.get("init-script", "init.py")`` -- Metis G3/G11: MR/LTS examples
  ship ``init_multi_resolution*.py`` and have no working ``run.py`` auto-init),
  forwarding every non-reserved config item as ``--<key> <value>`` in dict
  order with ``str(value)``.
- :func:`run_case` invokes ``run.py`` and returns its exit code (no streaming
  by default, same tolerance as ``multirunner.run``).
- :func:`snapshot_tracked_templates`/:func:`restore_tracked_templates`
  guard the TRACKED ``template/config.h`` rewritten by init (D8): tier runners
  snapshot before init+build and restore in ``finally`` (covers exceptions;
  the SIGKILL left-over is caught by the :func:`check_templates_clean`
  preflight -- Metis G18/H5).
- :func:`evaluate` resolves the ``evaluation-function`` STRING through
  ``regression.metrics.REGISTRY`` (D5; imported lazily INSIDE the function so
  this module loads before todo 8's registry exists).
- :func:`parse_comp_time` ports ``parse_tnl_sph_output`` (both the multirunner
  and testSolver variants key on ``total``): any failure -> ``0.0`` plus a
  warning naming the file. This is the ONLY function allowed to print.
- :func:`results_entry_for` builds a fresh :class:`suite.results.ResultEntry`
  (``group`` comes from the internal ``_group`` key when present).

Console printing is restricted to :func:`parse_comp_time`'s warning and the
``__main__`` selftest -- everything else stays silent for the tier runners.

Selftest (writes only into a temp dir)::

    PYTHONPATH=src/Tests python3 -m suite.runner_core
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from . import paths
from .config_loader import RESERVED_KEYS
from .results import ResultEntry


def init_case(case_dir: Path, conf: dict, verbose: bool = False) -> None:
    """Initialize a case by invoking its init script with forwarded kwargs.

    Port of ``multirunner.init`` (:35-48): argument order follows the conf's
    dict order, values are stringified with ``str``, and every key not in
    ``RESERVED_KEYS`` is forwarded as ``--<key> <value>``. Underscore-prefixed
    keys (internal suite plumbing such as ``_group``) are skipped as well --
    they must never reach a case's argparse.

    The case's case-tag is included in the raised error on failure.
    """
    case_dir = Path(case_dir)
    init_script = case_dir / conf.get("init-script", "init.py")
    args = ["python3", str(init_script)]
    for key, value in conf.items():
        if key in RESERVED_KEYS or key.startswith("_"):
            continue
        args += [f"--{key}", str(value)]

    if verbose:
        print(args)
    stdout = None if verbose else subprocess.DEVNULL
    try:
        subprocess.run(
            args,
            check=True,
            stdout=stdout,
            stderr=subprocess.STDOUT,
            cwd=case_dir,
        )
    except subprocess.CalledProcessError as exc:
        tag = conf.get("case-tag", "<unknown>")
        raise RuntimeError(
            f"init failed for case '{tag}' in {case_dir} "
            f"(init script {init_script.name}, exit code {exc.returncode})"
        ) from exc


def run_case(case_dir: Path, verbose: bool = False) -> int:
    """Run the case via its ``run.py`` and return the process exit code.

    Port of ``multirunner.run`` (:61-76): output is suppressed unless
    ``verbose`` is set; nothing is streamed or reformatted either way.
    """
    case_dir = Path(case_dir)
    args = ["python3", str(case_dir / "run.py")]
    stdout = None if verbose else subprocess.DEVNULL
    completed = subprocess.run(
        args,
        stdout=stdout,
        stderr=subprocess.STDOUT,
        cwd=case_dir,
        text=True,
    )
    return completed.returncode


def snapshot_tracked_templates(case_dir: Path) -> dict[Path, bytes]:
    """Snapshot the TRACKED template files init may rewrite (currently only
    ``template/config.h``, D8).

    Returns a mapping of path RELATIVE TO ``case_dir`` -> raw bytes, i.e.
    exactly ``{Path("template/config.h"): <bytes>}`` when the file exists.
    Missing files are skipped, so the dict may be empty.
    """
    case_dir = Path(case_dir)
    snapshot: dict[Path, bytes] = {}
    config_h = case_dir / "template" / "config.h"
    if config_h.is_file():
        snapshot[Path("template/config.h")] = config_h.read_bytes()
    return snapshot


def restore_tracked_templates(case_dir: Path, snapshot: dict[Path, bytes]) -> None:
    """Write back the bytes captured by :func:`snapshot_tracked_templates`.

    Snapshot keys are paths RELATIVE TO ``case_dir``; each entry is written
    to ``case_dir / rel``. Tier runners call this in ``finally`` around
    init+build so the tracked files are restored even on exceptions. A no-op
    for an empty snapshot.
    """
    case_dir = Path(case_dir)
    for rel, data in snapshot.items():
        (case_dir / rel).write_bytes(data)


def check_templates_clean(case_rel_paths: list[str] | None = None) -> list[str]:
    """Preflight (Metis G18/H5): return repo-relative paths of TRACKED
    ``template/config.h`` files with uncommitted modifications.

    Runs ``git -C <root> diff --name-only -- examples/**/template/config.h``
    (git's default pathspec wildcards match across directories, so ``**``
    covers every example depth). When ``case_rel_paths`` is a non-empty list,
    the check is restricted to those cases' ``template/config.h`` files; an
    empty list short-circuits to no dirty files.

    This is the ONLY suite function allowed to shell out to git. Tier runners
    abort with a message naming the returned files when it is non-empty.
    """
    if case_rel_paths is not None and len(case_rel_paths) == 0:
        return []
    if case_rel_paths:
        pathspecs = [
            str(Path("examples") / rel / "template" / "config.h")
            for rel in case_rel_paths
        ]
    else:
        pathspecs = ["examples/**/template/config.h"]
    completed = subprocess.run(
        ["git", "-C", str(paths.ROOT), "diff", "--name-only", "--"] + pathspecs,
        check=True,
        capture_output=True,
        text=True,
    )
    return [line for line in completed.stdout.splitlines() if line.strip()]


def evaluate(conf: dict, case_dir: Path) -> tuple[int, int]:
    """Evaluate the case's test metrics via the ``regression.metrics`` REGISTRY.

    ``conf["evaluation-function"]`` absent/``None``/``""`` -> ``(0, 0)`` (port
    of ``multirunner.evaluate_test_metrics``'s None rule). Otherwise the value
    MUST be a registry NAME (D5): it is looked up lazily
    (``from regression import metrics`` is imported INSIDE this function so
    this module loads before todo 8's registry exists) and called with
    ``case_dir``. Unknown names raise ``ValueError`` listing the valid
    registry names; raw callables are rejected with a ``ValueError``.

    Registry call rule (uniform contract, see
    ``regression.metrics.__init__``): every registry function has the
    signature ``fn(case_dir, *, collect=False, case_tag=None)`` and this
    function ALWAYS threads the entry's ``case-tag`` into the call, so a
    metric function serving several tags can resolve its per-tag committed
    reference -- falling back to its hardcoded base reference when the
    per-tag file is absent.  A missing ``case-tag`` arrives as ``None``,
    which selects the function's base reference.
    """
    name = conf.get("evaluation-function")
    if name is None or name == "":
        return 0, 0
    tag = conf.get("case-tag", "<unknown>")
    if not isinstance(name, str):
        raise ValueError(
            f"evaluation-function for case '{tag}' must be a STRING registry "
            f"name; raw callables are rejected (D5), got {type(name).__name__}"
        )
    try:
        from regression import metrics  # noqa: PLC0415 - lazy import (todo 8)
    except ImportError as exc:
        raise RuntimeError(
            f"cannot evaluate case '{tag}': failed to import "
            f"'regression.metrics' ({exc}); is src/Tests on sys.path?"
        ) from exc
    registry = getattr(metrics, "REGISTRY", None)
    if registry is None:
        raise RuntimeError(
            f"cannot evaluate case '{tag}': regression.metrics.REGISTRY is not "
            f"available yet (provided by todo 8)"
        )
    if name not in registry:
        valid = ", ".join(sorted(registry)) or "(empty registry)"
        raise ValueError(
            f"unknown evaluation function '{name}' for case '{tag}'. "
            f"Valid registry names: {valid}"
        )
    passed, total = registry[name](Path(case_dir), case_tag=conf.get("case-tag"))
    return int(passed), int(total)


def parse_comp_time(case_dir: Path) -> float:
    """Read ``results/time_measurements.json`` and return the ``total`` timer.

    Port of ``multirunner.parse_tnl_sph_output`` (:84-92) and the testSolver
    variant (same ``total`` key, same tolerance): ANY failure (missing file,
    invalid JSON, missing/NaN key) returns ``0.0`` and prints a warning naming
    the file -- the only print allowed in this module.
    """
    filename = Path(case_dir) / "results" / "time_measurements.json"
    try:
        with open(filename) as f:
            timers_dictionary = json.load(f)
        return float(timers_dictionary["total"])
    except Exception as exc:
        print(f"parse_comp_time: Could not read {filename}: {exc}")
        return 0.0


def results_entry_for(conf: dict, tier: str) -> ResultEntry:
    """Build a fresh :class:`ResultEntry` for ``conf`` in tier ``tier``.

    ``tag`` = the entry's ``case-tag``; ``group`` = the configuration module
    name threaded through the internal ``_group`` key (empty when absent).
    """
    return ResultEntry(
        tag=str(conf.get("case-tag", "")),
        group=str(conf.get("_group", "")),
        tier=str(tier),
    )


def _selftest() -> None:
    """Selftest: all primitives against synthetic case dirs in a temp dir."""
    import shutil
    import sys
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        # --- 1) init_case: argv construction, RESERVED/_-key exclusion, order.
        case_dir = tmp_path / "case"
        case_dir.mkdir()
        capture = tmp_path / "argv.txt"
        init_py = case_dir / "init.py"
        init_py.write_text(
            "import sys\n"
            "from pathlib import Path\n"
            f"Path({str(capture)!r}).write_text('\\n'.join(sys.argv[1:]))\n",
            encoding="utf-8",
        )
        conf = {
            "case": "DEMO/case",
            "case-tag": "DEMO/case:tag",
            "evaluation-function": "some_name",
            "init-script": "init.py",
            "target": "case_cuda",
            "init-overrides": {"x": 1},
            "timeout_s": 60,
            "requires": ["gpu"],
            "_group": "demo_group",
            "bc-type": "DBC",
            "h-coef": 1.5,
            "dp": 0.01,
        }
        init_case(case_dir, conf)
        forwarded = capture.read_text(encoding="utf-8").splitlines()
        assert forwarded == [
            "--bc-type",
            "DBC",
            "--h-coef",
            "1.5",
            "--dp",
            "0.01",
        ], forwarded  # dict order, str(value), reserved/_ keys excluded

        # custom init-script name is honored
        custom = case_dir / "init_multi_resolution.py"
        custom.write_text(
            "import sys\n"
            "from pathlib import Path\n"
            f"Path({str(capture)!r}).write_text('\\n'.join(sys.argv[1:]))\n",
            encoding="utf-8",
        )
        mr_conf = {"case": "DEMO/mr", "case-tag": "DEMO/mr:tag",
                   "init-script": "init_multi_resolution.py", "dp": 0.02}
        init_case(case_dir, mr_conf)
        assert capture.read_text(encoding="utf-8").splitlines() == ["--dp", "0.02"]

        # --- 2) init_case failure: error names the case-tag; init.py unchanged.
        fail_script = case_dir / "fail_init.py"
        fail_script.write_text("import sys; sys.exit(3)\n", encoding="utf-8")
        bad_conf = {"case": "DEMO/bad", "case-tag": "DEMO/bad:oops",
                    "init-script": "fail_init.py"}
        try:
            init_case(case_dir, bad_conf)
        except RuntimeError as exc:
            assert "DEMO/bad:oops" in str(exc), exc
            assert isinstance(exc.__cause__, subprocess.CalledProcessError)
        else:
            raise AssertionError("init_case must raise RuntimeError on non-zero exit")

        # --- 3) run_case returns the raw exit code (no exception).
        run_py = case_dir / "run.py"
        run_py.write_text("import sys; sys.exit(7)\n", encoding="utf-8")
        assert run_case(case_dir) == 7
        run_py.write_text("import sys; sys.exit(0)\n", encoding="utf-8")
        assert run_case(case_dir) == 0

        # --- 4) snapshot/restore round-trip of template/config.h.
        template_dir = case_dir / "template"
        template_dir.mkdir()
        config_h = template_dir / "config.h"
        original = b"// tracked template config\n#pragma once\n"
        config_h.write_bytes(original)
        snap = snapshot_tracked_templates(case_dir)
        assert snap == {Path("template/config.h"): original}, snap  # keys relative to case_dir
        config_h.write_bytes(b"\x00\xff JUNK junk \xfe\n")
        restore_tracked_templates(case_dir, snap)
        assert config_h.read_bytes() == original
        # missing template -> empty snapshot -> restore is a no-op
        bare_case = tmp_path / "bare_case"
        bare_case.mkdir()
        assert snapshot_tracked_templates(bare_case) == {}
        restore_tracked_templates(bare_case, {})

        # --- 5) parse_comp_time tolerance: valid total, missing file, garbage.
        results_dir = case_dir / "results"
        results_dir.mkdir()
        (results_dir / "time_measurements.json").write_text(
            json.dumps({"total": 12.5, "init": 1.0}), encoding="utf-8"
        )
        assert parse_comp_time(case_dir) == 12.5
        (results_dir / "time_measurements.json").unlink()
        assert parse_comp_time(bare_case) == 0.0  # warning printed, no raise
        (results_dir / "time_measurements.json").write_text("{not json", encoding="utf-8")
        assert parse_comp_time(case_dir) == 0.0
        (results_dir / "time_measurements.json").write_text(
            json.dumps({"other": 1}), encoding="utf-8"
        )
        assert parse_comp_time(case_dir) == 0.0  # missing 'total' key tolerated

        # --- 6) evaluate: absent/None/"" -> (0, 0); raw callable rejected.
        assert evaluate({"case-tag": "x"}, case_dir) == (0, 0)
        assert evaluate({"case-tag": "x", "evaluation-function": None}, case_dir) == (0, 0)
        assert evaluate({"case-tag": "x", "evaluation-function": ""}, case_dir) == (0, 0)
        try:
            evaluate({"case-tag": "x", "evaluation-function": lambda d: (0, 0)}, case_dir)
        except ValueError as exc:
            assert "raw callables are rejected" in str(exc), exc
        else:
            raise AssertionError("evaluate must reject raw callables with ValueError")

        # --- 7) results_entry_for: tag/group/tier mapping.
        entry = results_entry_for(dict(conf), tier="regression")
        assert entry.tag == "DEMO/case:tag"
        assert entry.group == "demo_group"
        assert entry.tier == "regression"
        assert entry.phases == {} and entry.error == ""
        no_group = results_entry_for({"case-tag": "t"}, tier="unit")
        assert no_group.group == "" and no_group.tier == "unit"

    # --- 8) check_templates_clean: empty list short-circuits; real repo check
    #        runs clean on this tree (dirty-window proof lives in the evidence).
    assert check_templates_clean([]) == []
    clean = check_templates_clean()
    assert clean == [], f"repo tree expected clean, got dirty: {clean}"

    _ = shutil, sys  # silence unused-import linters in the selftest
    print("runner_core selftest PASS")


if __name__ == "__main__":
    _selftest()
