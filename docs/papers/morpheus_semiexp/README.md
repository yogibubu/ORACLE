# MORPHEUS semiexperimental-structure paper

This directory contains the current MATRIX copy of the MORPHEUS manuscript and
supporting information.

Primary sources:

- `main_short.tex`: main manuscript.
- `supporting_information.tex`: supporting information.
- `references.bib`: bibliography used by the current manuscript.
- `figures/`: paper figures.
- `generated/`: generated tables included by the manuscript and SI.

Build commands:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main_short.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error supporting_information.tex
```

The committed `main_short.pdf` and `supporting_information.pdf` are compiled
from these sources.

The MORPHEUS benchmark summary tables can be regenerated with:

```bash
python -m matrix semiexp-benchmark --outdir docs/papers/morpheus_semiexp/generated
```

Some paper-specific generated tables in `generated/` come from dedicated
workflows and should be refreshed together with their corresponding input data.

See `reproducibility_audit.md` for the current command-by-command status of the
paper results.
