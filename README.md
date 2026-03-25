# HelixHarbor

HelixHarbor is a web application for the analysis of transmembrane proteins and their topological organization. It supports single-sequence topology analysis, comparison of UniProt-derived protein sets against curated background datasets, direct comparison of two UniProt accession lists, and position-specific amino-acid composition analysis across transmembrane segments.

The tool is also accessible through:

```text
https://service2.bioinformatik.uni-saarland.de/HelixHarbor/
```

## What HelixHarbor does

HelixHarbor currently provides four analysis modes:

- single-sequence topology and transmembrane-segment analysis
- comparison of a UniProt accession list against curated background datasets
- direct comparison of two UniProt accession lists
- position-specific amino-acid composition analysis of transmembrane segments

The web interface also supports export of the raw values used to generate comparison plots.

## Documentation

- [User guide](docs/USER_GUIDE.md)
- [Methods](docs/METHODS.md)

The user guide combines the practical manual, output interpretation notes, and the method-focused background that used to be spread across separate draft files.

## Repository layout

- `app.py`: Flask routes and request handling
- `function_set.py`: core analysis, plotting, and export logic
- `templates/`: web templates
- `backgrounds_tsv/`: curated background datasets used by the application
- `tests/`: regression tests
- `DCS/`: default custom-scale template
- `Dockerfile`: container build definition

## Running with Docker

Pull the latest published image:

```bash
docker pull mmofty/helixharbor:latest
```

Run the container:

```bash
docker run --restart=always -d --name helixharbor -p 5005:5005 mmofty/helixharbor:latest
```

Open the application at:

```text
http://localhost:5005/HelixHarbor
```

If port `5005` is already in use, bind a different host port:

```bash
docker run --restart=always -d --name helixharbor -p 5006:5005 mmofty/helixharbor:latest
```

## Running from source

The repository contains the application code, tests, curated background datasets, and a vendored copy of the TMbed source code.

Important: local TMbed model weights are intentionally excluded from version control. If you want to run sequence mode from source, you must provision the required TMbed assets separately.

Install Python dependencies:

```bash
pip install -r requirements.txt
```

Start the Flask/Gunicorn app according to your deployment setup, or build a Docker image from this repository.

## Testing

The regression suite currently covers the list-based analysis modes, including:

- compare-two-lists
- list-against-background
- empty-result handling
- raw dataframe normalization across mixed data sources

Run the tests with the published container environment:

```bash
docker run --rm --entrypoint python -v "$PWD":/app -w /app mmofty/helixharbor:latest -m unittest discover -s tests -v
```

## Contact

For questions related to HelixHarbor, contact `mohamed.elmofty2hu-berlin.de`.
