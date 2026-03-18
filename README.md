# HelixHarbor

HelixHarbor is a web application for transmembrane-region analysis, list-vs-background comparison, list-vs-list comparison, and export of raw values underlying generated plots.

## Repository scope

This repository contains the HelixHarbor application code, HTML templates, tests, Docker configuration, and background datasets required by the app.

Local runtime artifacts, generated reports, notebook checkpoints, virtual environments, and bundled TMbed model weights are intentionally excluded from version control.

## Main files

- `app.py`: Flask application and routes
- `function_set.py`: analysis and plotting logic
- `templates/`: web templates
- `examples/`: example input files used for manual validation
- `backgrounds_tsv/`: background datasets used by the app
- `tests/`: regression tests
- `Dockerfile`: container build definition

## TMbed note

Do not commit local TMbed model weights or the vendored `tmbed/.git` directory to this repository.

If TMbed is required for deployment, document its installation or container build separately rather than storing large model files in GitHub.
