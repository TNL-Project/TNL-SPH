"""Per-case boot data map of the infrastructure tier (todo 21).

``BOOT_CONFIGS`` enumerates every discovered example dir explicitly (Metis
G6: 19 dirs under ``examples/{RSPH,SHTC,WCSPH-BI,WCSPH-DBC}/`` carrying both
``case.h`` and ``template/``); :func:`infrastructure.infra_runner.discover_cases`
finds exactly these dirs by scanning for the same two markers.

Field semantics (read by ``infra_runner`` only):

- ``"init-script"``: the init script the boot invokes DIRECTLY (never
  ``run.py`` -- Metis G11: MR/LTS examples have no ``init.py`` and their
  ``run.py`` auto-init is broken).  Verified per case against the dir
  listing: MR dirs ship ``init_multi_resolution.py``, the LTS dir ships
  ``init_multi_resolution_lts.py``.
- ``"target"``: the CMake target built and executed.  Default rule is
  ``<leaf-dir>_cuda``; MR/LTS dirs and both ``*_distributed`` dirs carry
  target != dir name (Metis G12), verified against each case's own
  ``CMakeLists.txt`` ``add_executable`` line.
- ``"init-overrides"``: init-kwarg overrides for boot speed, each justified
  inline.  Only knobs verified as real argparse options of the case's own
  init script are used.
- ``"timeout_s"``: wall-clock cap of the BINARY boot run (init/build phases
  are uncapped).  180 s default, 360 s for 3D cases (bigger init+report).
- ``"requires"``: environment prerequisites checked by the preflight BEFORE
  anything else; ``["gencase"]`` means the case's init script shells out to
  ``src/tools/3rdparty/GenCase_linux64`` (verified by grepping the init
  scripts: exactly the four 3D dam-break cases below).
- ``"distributed"``: only on the two ``*_distributed`` dirs; ``np`` is the
  ``mpirun -np`` rank count, ``extra-init-overrides`` forces the subdomain
  split (Metis G5: both distributed init scripts default to
  ``--subdomains-x 3``, so a 2-rank boot must override it).
"""

# GenCase cases (verified by ``grep -in gencase examples/*/*/init*.py``):
#   WCSPH-DBC/damBreak3D_WCSPH-DBC              init.py
#   WCSPH-DBC/damBreak3D_WCSPH-DBC_distributed  init.py
#   WCSPH-BI/damBreak3D_WCSPH-BI                init.py
#   WCSPH-BI/damBreak3D_WCSPH-BI_multiresolution  init_multi_resolution.py

BOOT_CONFIGS: dict[str, dict] = {
    # ------------------------------------------------------------------ RSPH
    "RSPH/damBreak2D_RSPH": {
        "init-script": "init.py",
        "target": "damBreak2D_RSPH_cuda",
        "init-overrides": {},  # 2D small case: default dp boots fast enough
        "timeout_s": 180,
        "requires": [],
    },
    # ------------------------------------------------------------------ SHTC
    "SHTC/lidDrivenCavity2D_SHTC": {
        "init-script": "init.py",
        "target": "lidDrivenCavity2D_SHTC_cuda",
        "init-overrides": {},  # 2D small case: default dp boots fast enough
        "timeout_s": 180,
        "requires": [],
    },
    # -------------------------------------------------------------- WCSPH-BI
    "WCSPH-BI/damBreak2D_WCSPH-BI": {
        "init-script": "init.py",
        "target": "damBreak2D_WCSPH-BI_cuda",
        "init-overrides": {},  # 2D small case: default dp boots fast enough
        "timeout_s": 180,
        "requires": [],
    },
    "WCSPH-BI/damBreak2D_WCSPH-BI_multiresolution": {
        "init-script": "init_multi_resolution.py",  # pinned by T7: the dir
        # also ships init_single_resolution.py -- do NOT use it
        "target": "damBreak2D_WCSPH-BI_MR_cuda",  # per CMakeLists.txt:2
        "init-overrides": {},  # keep default dp 0.002: already fixture-scale
        # (init's own argparse default, verified at init_multi_resolution.py:324)
        "timeout_s": 180,
        "requires": [],
    },
    "WCSPH-BI/damBreak2D_WCSPH-BI_multiresolution_local-timestepping": {
        "init-script": "init_multi_resolution_lts.py",  # pinned by T7:
        # NO init.py exists in this dir at all
        "target": "damBreak2D_WCSPH-BI_MR_LTS_cuda",  # per CMakeLists.txt:2
        "init-overrides": {},  # keep default dp 0.002: already fixture-scale
        # (init's own argparse default, verified at init_multi_resolution_lts.py:314)
        "timeout_s": 180,
        "requires": [],
    },
    "WCSPH-BI/damBreak3D_WCSPH-BI": {
        "init-script": "init.py",
        "target": "damBreak3D_WCSPH-BI_cuda",
        "init-overrides": {
            # dp: 2x coarser than the init default 0.02 (init.py:171) for
            # boot speed; 0.04 is the task-pinned minimum for 3D dam breaks
            "dp": 0.04,
        },
        "timeout_s": 360,  # 3D case
        "requires": ["gencase"],  # generate_geometry_with_dualsphysics_gencase
        # shells out to src/tools/3rdparty/GenCase_linux64 (init.py:11-13)
    },
    "WCSPH-BI/damBreak3D_WCSPH-BI_multiresolution": {
        "init-script": "init_multi_resolution.py",  # pinned by T7
        "target": "damBreak3D_WCSPH-BI_MR_cuda",  # per CMakeLists.txt:2
        "init-overrides": {},  # keep default coarse dp 0.04: already the
        # task-pinned 3D coarsening (init_multi_resolution.py:266)
        "timeout_s": 360,  # 3D case
        "requires": ["gencase"],  # run_gencase()
        # (init_multi_resolution.py:43-47)
    },
    "WCSPH-BI/movingSquare2D_WCSPH-BI": {
        "init-script": "init.py",
        "target": "movingSquare2D_WCSPH-BI_cuda",
        "init-overrides": {},  # 2D small case: default dp boots fast enough
        "timeout_s": 180,
        "requires": [],
    },
    "WCSPH-BI/openChannel3D_WCSPH-BI": {
        "init-script": "init.py",
        "target": "openChannel3D_WCSPH-BI_cuda",
        "init-overrides": {
            # dp: default 0.005 (init.py:358) doubled to 0.01 for boot speed
            "dp": 0.01,
        },
        "timeout_s": 360,  # 3D case
        "requires": [],
    },
    "WCSPH-BI/poiseuilleFlowWithOpenBoundary2D_WCSPH-BI": {
        "init-script": "init.py",
        "target": "poiseuilleFlowWithOpenBoundary2D_WCSPH-BI_cuda",
        "init-overrides": {},  # 2D small case: default dp boots fast enough
        "timeout_s": 180,
        "requires": [],
    },
    "WCSPH-BI/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-BI": {
        "init-script": "init.py",
        "target": "poiseuilleFlowWithPeriodicBoundary2D_WCSPH-BI_cuda",
        "init-overrides": {},  # 2D small case: default dp boots fast enough
        "timeout_s": 180,
        "requires": [],
    },
    # ------------------------------------------------------------- WCSPH-DBC
    "WCSPH-DBC/damBreak2D_WCSPH-DBC": {
        "init-script": "init.py",
        "target": "damBreak2D_WCSPH-DBC_cuda",
        "init-overrides": {},  # 2D small case: default dp boots fast enough
        "timeout_s": 180,
        "requires": [],
    },
    "WCSPH-DBC/damBreak2D_WCSPH-DBC_distributed": {
        "init-script": "init.py",
        "target": "damBreak2D_WCSPH-DBC_mpi_cuda",  # per CMakeLists.txt:2;
        # NOTE: examples/WCSPH-DBC/CMakeLists.txt:2 comments out the
        # add_subdirectory for this dir, so the target is absent from the
        # current build tree (boot handles this as SKIP with that reason)
        "init-overrides": {},  # 2D: keep default dp 0.002 (init.py:495)
        "timeout_s": 180,
        "requires": [],
        "distributed": {
            "np": 2,
            "extra-init-overrides": {
                # both distributed init scripts default --subdomains-x 3
                # (Metis G5; init.py:491), a 2-rank boot must override it
                "subdomains-x": 2,
                "subdomains-y": 1,
            },
        },
    },
    "WCSPH-DBC/damBreak3D_WCSPH-DBC": {
        "init-script": "init.py",
        "target": "damBreak3D_WCSPH-DBC_cuda",
        "init-overrides": {
            # dp: 2x coarser than the init default 0.02 (init.py:160) for
            # boot speed; 0.04 is the task-pinned minimum for 3D dam breaks
            "dp": 0.04,
        },
        "timeout_s": 360,  # 3D case
        "requires": ["gencase"],  # generate_geometry_with_dualsphysics_gencase
        # shells out to src/tools/3rdparty/GenCase_linux64 (init.py:11-13)
    },
    "WCSPH-DBC/damBreak3D_WCSPH-DBC_distributed": {
        "init-script": "init.py",
        "target": "damBreak3D_WCSPH-DBC_mpi_cuda",  # per CMakeLists.txt:2
        "init-overrides": {
            # dp: 2x coarser than the init default 0.02 (init.py:122) for
            # boot speed; 0.04 is the task-pinned minimum for 3D dam breaks
            "dp": 0.04,
        },
        "timeout_s": 360,  # 3D case
        "requires": ["gencase"],  # generate_geometry()
        # shells out to src/tools/3rdparty/GenCase_linux64 (init.py:20-24)
        "distributed": {
            "np": 2,
            "extra-init-overrides": {
                # init.py:117 defaults --subdomains-x 3 (Metis G5)
                "subdomains-x": 2,
                "subdomains-y": 1,
            },
        },
    },
    "WCSPH-DBC/movingSquare2D_WCSPH-DBC": {
        "init-script": "init.py",
        "target": "movingSquare2D_WCSPH-DBC_cuda",
        "init-overrides": {},  # 2D small case: default dp boots fast enough
        "timeout_s": 180,
        "requires": [],
    },
    "WCSPH-DBC/openChannel3D_WCSPH-DBC": {
        "init-script": "init.py",  # the ONLY scripted variant exercised;
        # init_wide.py/init_pipe.py geometry variants stay EXCLUDED
        "target": "openChannel3D_WCSPH-DBC_cuda",
        "init-overrides": {
            # dp: default 0.005 (init.py:393) doubled to 0.01 for boot speed
            "dp": 0.01,
        },
        "timeout_s": 360,  # 3D case
        "requires": [],
    },
    "WCSPH-DBC/poiseuilleFlowWithOpenBoundary2D_WCSPH-DBC": {
        "init-script": "init.py",
        "target": "poiseuilleFlowWithOpenBoundary2D_WCSPH-DBC_cuda",
        "init-overrides": {},  # 2D small case: default dp boots fast enough
        "timeout_s": 180,
        "requires": [],
    },
    "WCSPH-DBC/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-DBC": {
        "init-script": "init.py",
        "target": "poiseuilleFlowWithPeriodicBoundary2D_WCSPH-DBC_cuda",
        "init-overrides": {},  # 2D small case: default dp boots fast enough
        "timeout_s": 180,
        "requires": [],
    },
}

# Init-script variants deliberately excluded from the boot matrix (handled
# inline in the rows above: openChannel3D's init_wide.py/init_pipe.py
# geometry variants -- the default init.py is the only exercised one).
EXCLUDE: list[str] = []
