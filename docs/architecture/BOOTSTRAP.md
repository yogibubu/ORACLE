# Bootstrap

This skeleton is intentionally light. It should prove package boundaries before
scientific code is moved.

## Local Smoke Test

```bash
cd /Users/vincenzobarone/ORACLE
PYTHONPATH=packages/oracle-core/src python -m pytest
```

## Shell Helpers

ORACLE mirrors the Merlino `merlino-set`, `merlino-run`, `merlino-run-check`
and test helpers with sourceable shell functions:

```bash
source /Users/vincenzobarone/ORACLE/scripts/oracle_env.sh
oracle-set
oracle-run --help
oracle-run-check
oracle-test-all
```

See `docs/architecture/ORACLE_ENVIRONMENT.md`.

## Create A Workspace

```bash
cd /Users/vincenzobarone/ORACLE
PYTHONPATH=packages/oracle-core/src python tools/oracle_run.py init /tmp/oracle-demo
```

## Delegate To Merlino During Migration

When the Merlino repository is importable in `PYTHONPATH`, the temporary CLI can
delegate:

```bash
PYTHONPATH=/Users/vincenzobarone/merlino3.0 python tools/oracle_run.py merlino --help
```
