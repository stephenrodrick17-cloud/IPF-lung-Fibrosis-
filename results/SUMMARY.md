# GEO-derived summary for GSE24206

## Verified sample counts from the GEO series matrix

The workflow was run against the actual GEO matrix file cached locally in this environment. The metadata extraction script verified the following:

- Total samples: 23
- IPF: 17
- Control: 6
- Platform: GPL570
- GEO accession: GSE24206

## Thesis numbers vs actual GEO numbers

The thesis text quotes the cohort as 9 IPF and 6 controls. The raw GEO series matrix does not support that exact split; it contains 17 IPF and 6 controls in the full public series. The raw series matrix also contains early and advanced IPF labels, which were merged into a single IPF group for the analysis pipeline.

## GEO2R-compatible interpretation

This dataset is consistent with a GEO2R / limma design in which the primary comparison is:

- group = IPF vs Control
- model: ~ group
- method: limma with empirical Bayes moderation

## Project note

This environment had a blocked GEOquery/GEO web access path, so the pipeline uses the directly downloaded series matrix file as the trusted source of record. This is the same data that GEO2R would ingest from the GEO repository.
