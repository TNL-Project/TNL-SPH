"""Build manager of the TNL-SPH unified test suite.

Wraps CMake configure/build invocations for examples and C++ test fixtures.
Hard rules honoured here (AGENTS.md):

* never more than 2 parallel build jobs by default (template-heavy CUDA
  instantiation gets ``cicc``/``cc1plus`` OOM-killed above that),
* always ``cmake --build <dir> --target <t>``, never plain ``make``,
* an existing CMake cache is NEVER reconfigured - the user's build dir is
  authoritative.

CUDA architecture selection mirrors the CI fallback (.gitlab-ci.yml):
``native`` when a GPU is visible via ``nvidia-smi``, ``75`` otherwise
(GPU-less configure).

Usage: ``python3 -m suite.build`` runs a DRY selftest that prints (but does
not execute) the exact commands the orchestrator would run.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from . import paths

OUTPUT_TAIL_LINES = 40

_DEFAULT_JOBS = 2


def set_jobs(n: int) -> None:
    """Set the module-level default job count (AGENTS hard rule: 2)."""
    if n < 1:
        raise ValueError(f"set_jobs: expected >= 1, got {n}")
    global _DEFAULT_JOBS
    _DEFAULT_JOBS = n


def get_jobs() -> int:
    """Return the module-level default job count (2 unless overridden)."""
    return _DEFAULT_JOBS


def _have_gpu() -> bool:
    """Return True iff ``nvidia-smi -L`` succeeds (GPU visible)."""
    try:
        result = subprocess.run(
            ["nvidia-smi", "-L"], capture_output=True,
        )
    except OSError:
        return False
    return result.returncode == 0


def _tail(text: str, lines: int = OUTPUT_TAIL_LINES) -> str:
    """Return the last ``lines`` lines of ``text`` (for error diagnostics)."""
    return "\n".join(text.splitlines()[-lines:])


def _cache_build_type(build_dir: Path) -> str | None:
    """Return the cached CMAKE_BUILD_TYPE value of ``build_dir``, or None.

    Parses the ``CMAKE_BUILD_TYPE:<type>=<val>`` line of
    ``<build_dir>/CMakeCache.txt``; None when the file or entry is absent.
    """
    cache_file = Path(build_dir) / "CMakeCache.txt"
    if not cache_file.is_file():
        return None
    for line in cache_file.read_text().splitlines():
        if line.startswith("CMAKE_BUILD_TYPE:"):
            _key, _sep, value = line.partition("=")
            return value
    return None


def _configure_command(build_dir: Path, build_type: str) -> list[str]:
    """Assemble the cmake configure command used by ``ensure_configured``."""
    arch = "native" if _have_gpu() else "75"
    return [
        "cmake", "-B", str(build_dir), "-S", str(paths.ROOT),
        f"-DCMAKE_BUILD_TYPE={build_type}",
        f"-DCMAKE_CUDA_ARCHITECTURES={arch}",
    ]


def _build_command(build_dir: Path, jobs: int, target: str) -> list[str]:
    """Assemble the cmake build command used by the ``build_*`` helpers."""
    return [
        "cmake", "--build", str(build_dir), "--parallel", str(jobs),
        "--target", target,
    ]


def ensure_configured(build_dir: Path | None = None,
                      build_type: str = "Release") -> Path:
    """Return a configured ``build_dir``, configuring it when necessary.

    Reuse rule: if ``<build_dir>/CMakeCache.txt`` exists AND its cached
    CMAKE_BUILD_TYPE matches ``build_type``, the directory is returned
    unchanged - an existing cache is NEVER reconfigured. Otherwise run
    ``cmake -B <bd> -S <root> -DCMAKE_BUILD_TYPE=<type>
    -DCMAKE_CUDA_ARCHITECTURES=<native|75>`` (native when a GPU is visible,
    75 otherwise). Raises RuntimeError with the output tail on failure.
    """
    bd = paths.BUILD_DIR if build_dir is None else Path(build_dir)
    if _cache_build_type(bd) == build_type:
        return bd
    if bd.exists() and not bd.is_dir():
        raise RuntimeError(f"ensure_configured: {bd} exists and is not a directory")
    cmd = _configure_command(bd, build_type)
    result = subprocess.run(
        cmd, cwd=paths.ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"ensure_configured: cmake configure failed (rc={result.returncode}):\n"
            f"  command: {' '.join(cmd)}\n{_tail(result.stdout or '')}"
        )
    return bd


def _target_name(case_rel_path: str, target: str | None) -> str:
    """Resolve the CMake target for a case given its examples-relative path.

    Default rule: ``<final path segment>_cuda`` (leaf-dir rule), e.g.
    ``examples/WCSPH-DBC/damBreak2D_WCSPH-DBC`` -> ``damBreak2D_WCSPH-DBC_cuda``.
    An explicit ``target`` always wins.

    Metis G12 exception: multiresolution cases live in directories whose name
    does NOT match their binary target, e.g.
    ``examples/WCSPH-BI/damBreak2D_WCSPH-BI_multiresolution`` builds
    ``damBreak2D_WCSPH-BI_MR_cuda`` - such cases MUST pass explicit targets.
    """
    if target is not None:
        return target
    leaf = Path(case_rel_path.rstrip("/")).name
    return f"{leaf}_cuda"


def _examples_rel(case_rel_path: str) -> Path:
    """Normalize ``case_rel_path`` to a path relative to ``examples/``.

    Accepts both ``examples/WCSPH-DBC/damBreak2D_WCSPH-DBC`` and
    ``WCSPH-DBC/damBreak2D_WCSPH-DBC`` so every call site stays uniform.
    """
    rel = Path(case_rel_path)
    if rel.parts and rel.parts[0] == "examples":
        rel = Path(*rel.parts[1:])
    return rel


def _run_build(target: str, jobs: int | None) -> tuple[int, str]:
    """Run ``cmake --build`` for ``target``; return (returncode, output tail)."""
    resolved_jobs = get_jobs() if jobs is None else jobs
    cmd = _build_command(paths.BUILD_DIR, resolved_jobs, target)
    result = subprocess.run(
        cmd, cwd=paths.ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True,
    )
    tail = _tail(result.stdout or "")
    if result.returncode != 0:
        raise RuntimeError(
            f"build target '{target}' failed (rc={result.returncode}):\n"
            f"  command: {' '.join(cmd)}\n{tail}"
        )
    return result.returncode, tail


def build_target(case_rel_path: str, jobs: int | None = None,
                 target: str | None = None) -> tuple[int, str]:
    """Build an example case's CMake target; return (returncode, output tail).

    The target is resolved by ``_target_name`` (leaf-dir rule / explicit).
    Runs with ``cwd=paths.ROOT`` so it is CWD-independent. ``jobs`` defaults
    to the module-level setting (2, AGENTS hard rule). Raises RuntimeError
    with the output tail on non-zero returncode.
    """
    return _run_build(_target_name(case_rel_path, target), jobs)


def build_fixture(fixture_target: str, jobs: int | None = None) -> tuple[int, str]:
    """Build a C++ test fixture target literally; return (returncode, tail).

    Fixture targets carry NO ``_cuda`` suffix (e.g.
    ``dummyMultiresolutionSimulation2D``), hence the name is taken verbatim.
    Raises RuntimeError with the output tail on non-zero returncode.
    """
    return _run_build(fixture_target, jobs)


def case_binary(case_rel_path: str, target: str | None = None) -> Path:
    """Return the expected build-output path of a case binary.

    ``build/examples/<case rel. to examples/>/<resolved target>`` (mirrors
    the source tree). Existence is NOT asserted here - check ``.exists()``
    after a successful ``build_target``.
    """
    return (paths.BUILD_DIR / "examples" / _examples_rel(case_rel_path)
            / _target_name(case_rel_path, target))


def _selftest() -> int:
    """DRY selftest: print (NOT execute) the exact commands; exit 0."""
    tmp_build = Path("/tmp/tnl-sph-build-selftest")
    dbc_case = "examples/WCSPH-DBC/damBreak2D_WCSPH-DBC"
    mr_case = "examples/WCSPH-BI/damBreak2D_WCSPH-BI_multiresolution"
    mr_target = "damBreak2D_WCSPH-BI_MR_cuda"

    print("# suite.build dry-run selftest - prints commands, executes nothing")
    print()
    print(f"# ensure_configured(Path('{tmp_build}')), fresh tmp dir:")
    print("  " + " ".join(_configure_command(tmp_build, "Release")))
    print(f"# cwd: {paths.ROOT}")
    print()
    print(f"# build_target({dbc_case!r}):")
    print("  " + " ".join(
        _build_command(paths.BUILD_DIR, get_jobs(), _target_name(dbc_case, None))))
    print(f"# cwd: {paths.ROOT}")
    print()
    print(f"# build_target({mr_case!r}, target={mr_target!r}):")
    print("  " + " ".join(
        _build_command(paths.BUILD_DIR, get_jobs(),
                       _target_name(mr_case, mr_target))))
    print(f"# cwd: {paths.ROOT}")
    return 0


if __name__ == "__main__":
    unknown = [a for a in sys.argv[1:] if a not in {"--dry", "--selftest"}]
    if unknown:
        print(f"suite.build: unknown arguments: {unknown}", file=sys.stderr)
        raise SystemExit(2)
    raise SystemExit(_selftest())
