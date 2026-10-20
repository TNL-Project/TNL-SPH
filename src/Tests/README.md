# TNL-SPH test suite

One orchestrator (`runTests.py`), three serial tiers (`unit`, `infrastructure`,
`regression`), plus CTest glue for the `build/` tree; one command, one summary,
one session log per run (`logs/log_<stamp>.{json,md}`). Run `unit` and
`infrastructure` locally; the full `regression` matrix is cluster-scale by design.

## Quick start

```bash
python3 src/Tests/runTests.py --list                                    # crumbs only, nothing executes
python3 src/Tests/runTests.py --tier unit                               # python-only, no GPU
python3 src/Tests/runTests.py --tier infrastructure                     # boots every example + MR fixtures
python3 src/Tests/runTests.py --tier regression --conf configurations/qa_smoke.py
python3 src/Tests/runTests.py --report src/Tests/logs/log_<stamp>.json # reprint a stored session
ctest --test-dir build -L unit                                          # unit tier + DistributedParticleSystemTest
```

## Tiers

- `unit`: pure-python tooling and geometry checks (decomposition, config writer, VTK roundtrip, config schema) plus the C++ `DistributedParticles` test.
- `infrastructure`: boots every discovered example (incl. a 2-rank distributed boot) and the 21 MR fixtures; validates init geometry recorded in `results/init_report*.json`.
- `regression`: config-driven solver runs; named metrics from `regression/metrics/` compared against frozen references (`value` + `atol`/`rtol`).

## Layout

- `suite/` shared library: paths, bootstrap, build, results, reporting, config loader, CLI
- `unit/` python check groups + `cpp/` (DistributedParticles)
- `infrastructure/` boot validators, MR fixture machinery + `fixtures/`, boot driver (incl. distributed boot)
- `regression/` `configurations/`, `metrics/`, `references/` (+ `_staging/`), `runner.py`, `perf.py`
- `logs/` session logs `log_<stamp>.{json,md}` (gitignored)
- `runTests.py` thin entry point (all logic in `suite.cli`)

## Current state / known caveats

- Builds default to `-j2` (AGENTS.md hard rule; template-heavy, OOM above it) and an existing `build/` cache is never reconfigured. Raise with `--jobs N` carefully.
- `unit` is the only tier runnable without a GPU.
- The full `regression` matrix is cluster-scale by design. Local sanity runs are capped at ≤ 10 time steps: set a shorter end-time in the generated `sources/config.jsonc`.
- Tags with a reference only under `references/_staging/` report SKIP `unfrozen reference (staging only)`. Printed, never a silent failure.
- A few pre-existing branch issues are listed in AGENTS.md (TESTING section); see there instead of this file.
- `ctest -L unit` runs the unit tier plus `DistributedParticleSystemTest`.

## Freezing references (manual, cluster-side)

1. Cluster run: `python3 src/Tests/runTests.py --tier regression --conf suite_all --store-results` (backs up `results/` to `results_<tag>/` on pass).
2. Eyeball the `results_<tag>/` output.
3. `PYTHONPATH=src/Tests python3 -m regression.runner --conf configurations/suite_all.py --generate-references [--from-stored]` writes `references/_staging/<sanitized-tag>.json`.
4. Eyeball staged values, set `atol` per metric, drop `tolerance_pending`, set `verified: true`.
5. `git mv` the file from `references/_staging/` to `references/` and commit; the gate is active from then on.

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | all entries passed (SKIPs never fail a session) |
| 1 | at least one FAIL entry, or a tests mismatch (`tests_passed != tests_total`) |
| 2 | harness error: bad flag/group/conf, missing prebuilt binary under `--no-build`, unexpected exception |
