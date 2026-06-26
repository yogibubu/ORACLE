# ORACLE Developer Workflow

Every new workflow should be added in this order:

1. Define a non-Qt service API with dataclasses for normalized inputs/outputs.
2. Add CLI access under `tools/oracle_run.py` or the future `python -m oracle`.
3. Write an `oracle.run.v1` manifest with input/output checksums.
4. Add focused tests using small fixtures.
5. Add GUI wiring only after the service and CLI are tested.
6. Document file contracts and benchmark commands.
7. Run `oracle-test-all` or `python -m pytest` before committing.

Use the sourceable shell helpers when working interactively:

```bash
source /Users/vincenzobarone/ORACLE/scripts/oracle_env.sh
oracle-set
oracle-run --help
oracle-test-all
```

Rules:

- GUI classes must not contain scientific algorithms.
- Gaussian, Fortran and other external formats are adapters, not internal data
  models.
- Prefer typed ORACLE errors from `oracle_core.errors`.
- Important numerical fits should use shared ORACLE numerical primitives for
  damped normal equations, step limiting and rank/condition diagnostics unless
  a workflow has a documented reason to use a specialized solver.
- Use `tests/fixtures/test_molecules` for demanding GICForge/parser
  regressions instead of inventing ad hoc molecule inputs.
- Store new project outputs under `inputs/`, `runs/`, `outputs/`, `reports/`,
  `cache/` or `logs/`.
