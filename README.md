# pLDDT Change Analysis Across Structure Prediction Methods

## Overview

This repository contains a Python analysis script used to compare
**protein structure confidence (pLDDT)** between two conditions
("before" vs "after") across multiple structure prediction methods:

-   Protenix\
-   Protenix_1\
-   ESM\
-   AlphaFold 3 (AF3)

The analysis focuses on **changes in the number of high-confidence
residues (pLDDT ≥ 80)** per protein and summarizes results for a curated
set of **461 reference proteins**.

The main output is a **stacked bar chart** showing how many proteins
improve, worsen, or remain unchanged for each method.

------------------------------------------------------------------------

## What the Script Does

1.  **Loads and harmonizes input datasets**
    -   Accepts either two separate TSV files (before/after) or a single
        TSV with a `source` column (`_old` / `_new`)
    -   Extracts and normalizes gene/protein IDs
    -   Merges before/after entries by gene ID
2.  **Processes pLDDT statistics**
    -   Converts residue-level pLDDT counts to numeric values
    -   Computes cumulative residue counts for thresholds:
        -   ≥50, ≥60, ≥70, ≥80, ≥90
3.  **Restricts analysis to 461 reference proteins**
    -   Uses `461.csv` as the reference list
    -   Normalizes IDs for robust matching
4.  **Quantifies per-protein changes**
    -   Compares ≥80 residue counts before vs after
    -   Classifies proteins as:
        -   Increased confidence
        -   Decreased confidence
        -   No change
5.  **Generates a visualization**
    -   Stacked bar plot summarizing changes by method
    -   Counts are annotated directly on the bars
6.  **Diagnostics and validation**
    -   Reports which reference proteins are missing from ESM
    -   Allows inspection of individual proteins for debugging

------------------------------------------------------------------------

## Outputs

### Files

-   `stacked_bar_461_proteins.png`\
    A publication-ready stacked bar chart summarizing confidence
    changes.

### Console Output

-   Number of merged proteins per dataset
-   Counts of reference proteins found per method
-   Up / down / unchanged protein counts
-   Lists of proteins missing from ESM
-   Optional full-row output for specific proteins

------------------------------------------------------------------------

## Input Requirements

### Software

-   Python 3.8 or later

### Python packages

``` bash
pip install pandas numpy matplotlib pypandoc
```

### Input files

-   TSV files containing:
    -   `gene_name` or `gene` columns
    -   `residue_plddt_50`, `60`, `70`, `80`, `90`
-   `461.csv` containing a column named `old_id`

File paths are defined directly in the script and can be edited as
needed.

------------------------------------------------------------------------

## How to Run

1.  Update file paths in the script if necessary

2.  Run the analysis:

    ``` bash
    python compare_plddt_changes.py
    ```

3.  Inspect:

    -   `stacked_bar_461_proteins.png`
    -   Console summaries for diagnostics

------------------------------------------------------------------------

## Typical Use Case

This analysis is intended for: - Comparing structural confidence
improvements across prediction pipelines - Benchmarking new methods
against existing ones - Assessing coverage and confidence changes for a
defined protein set

------------------------------------------------------------------------

## Notes

-   Gene IDs with duplicates are excluded to avoid ambiguous merges
-   ESM-specific naming conventions are handled automatically
-   Absolute residue counts are used (consistent with Protenix output)

------------------------------------------------------------------------

## Contact / Maintenance

If you modify the input formats or add new prediction methods, ensure: -
Gene ID extraction remains consistent - pLDDT threshold columns follow
the expected naming convention
