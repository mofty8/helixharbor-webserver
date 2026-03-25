# HelixHarbor User Guide

## Overview

HelixHarbor is a web application for analyzing transmembrane helices and related topological regions in transmembrane proteins. The current interface supports four main workflows:

1. Sequence analysis
2. List against background comparison
3. Compare two lists
4. TMH position-specific amino-acid composition

The public web service is available at:

```text
https://service2.bioinformatik.uni-saarland.de/HelixHarbor/
```

## Analysis modes

### 1. Sequence

Use this mode when you want to analyze a single amino-acid sequence.

Input:

- one amino-acid sequence

Output:

- region annotation
- segment coordinates
- inferred topology labels
- per-segment feature values such as surface area, volume, bulkiness, and hydrophobicity

Implementation note:

- sequence mode uses TMbed-based topology prediction

### 2. List Against Background

Use this mode when you want to compare a UniProt accession list against a curated reference set.

Input:

- a list of UniProt accessions
- organism selection
- background dataset selection
- TMH comparison preference: `First TMH` or `All TMHs`
- feature selection

Supported background options in the current interface:

- UniProt-derived background
- TopDB
- TMalpha

Supported features:

- Volume
- Surface Area
- Bulkiness
- Hydrophobicity
- Amino Acid Composition Against Background
- Custom Scale

Output:

- comparison plot between the list of interest and the selected background
- optional non-transmembrane warning list when relevant
- downloadable raw plot values

### 3. Compare 2 Lists

Use this mode when you want to compare two UniProt accession lists directly.

Input:

- list 1: UniProt accessions
- list 2: UniProt accessions
- TMH comparison preference: `First TMH` or `All TMHs`
- feature selection

Output:

- comparison plot for the selected feature
- downloadable raw plot values

Current behavior update:

- this mode now uses two separate text boxes as expected
- handling of mixed background-derived and UniProt-derived records has been stabilized

### 4. TMH Position Specific AAC

Use this mode when you want a position-wise amino-acid composition summary across transmembrane segments.

Input:

- a list of UniProt accessions
- orientation preference: `inside` or `outside`

Output:

- amino-acid composition heatmap across TMH positions

## Downloads

Depending on the selected analysis, HelixHarbor can generate:

- sequence annotation table
- background dataset download
- custom-scale template download
- raw plot data export


## Raw plot export

For list-based comparisons, HelixHarbor can export the values used to generate the plot.

The current raw plot export format is:

```text
group, id, type, begin, end, value
```

Field meanings:

- `group`: comparison group, for example `List 1`, `List 2`, `List of Interest`, or `Background`
- `id`: UniProt accession
- `type`: segment type recorded for the plotted region
- `begin`, `end`: 1-based coordinates of the contributing region
- `value`: numeric feature value used in plotting

This file is suitable for downstream use in other plotting software after reshaping by group.

## Plot interpretation

For the list-based comparison modes, HelixHarbor visualizes distributions using kernel density estimates.

Important interpretation points:

- the x-axis represents the selected feature value
- the y-axis represents estimated probability density, not raw protein count
- each plotted point corresponds to one transmembrane segment or one first-TMH value, depending on the selected mode

## Background datasets

HelixHarbor uses curated datasets from multiple sources, including:

- UniProt-derived annotations
- TOPDB
- TMalpha

These datasets are used to provide comparison baselines and annotated transmembrane-region coordinates.

## Physicochemical features

HelixHarbor computes per-segment averages from residue-level scales.

### Surface area

- based on residue accessibility reference values
- summarized as a mean per residue across the selected segment

Reference:

- Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO. *PLoS ONE* 8(11):e80635, 2013

### Volume

- summarized as mean residue volume across the selected segment

Reference:

- Chothia C. *Nature* 254(5498):304-308, 1975

### Bulkiness

- summarized as mean residue bulkiness across the selected segment

Reference:

- Zimmerman JM, Eliezer N, Simha R. *Journal of Theoretical Biology* 21(2):170-201, 1968

### Hydrophobicity

- based on the Kyte-Doolittle hydropathy scale
- summarized as a mean per residue across the selected segment

Reference:

- Kyte J, Doolittle RF. *Journal of Molecular Biology* 157(1):105-132, 1982

### Custom scale

- user-supplied amino-acid values
- summarized as a mean per residue across the selected segment

## Method notes

### Sequence mode

- TMbed is used to assign per-residue classes
- contiguous runs are collapsed into annotated segments with start and end coordinates

### UniProt-based modes

- UniProt accessions are resolved against curated background data and, where needed, fetched feature annotations
- transmembrane segment sequences are extracted using annotated coordinates

### Statistics and visualization

- list-based distributions are visualized as kernel density estimates
- the codebase uses the SciPy and Matplotlib ecosystem, with Seaborn used for relevant visualizations

## Contact

For questions or support, contact `mohamed.elmofty2hu-berlin.de`.
