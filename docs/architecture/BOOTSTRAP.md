# Bootstrap

This skeleton is intentionally light. It should prove package boundaries before
scientific code is moved.

## Local Smoke Test

```bash
cd /Users/vincenzobarone/ORACLE
PYTHONPATH=packages/oracle-core/src python -m pytest
```

## Shell Helpers

ORACLE provides sourceable shell helpers for environment setup, launch checks
and tests:

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

## Inspect The CLI

The temporary CLI is available directly from the repository:

```bash
python tools/oracle_run.py --help
```
