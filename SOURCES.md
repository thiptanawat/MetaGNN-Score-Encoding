# Sources

## Third-party inputs, not redistributed

Each must be retrieved once and placed where `reproduce.sh` expects it. The expected SHA-256 is
given so that a retrieved copy can be checked against the one this study used.

| Input | Route | SHA-256 |
|---|---|---|
| Recon3D reconstruction | `http://bigg.ucsd.edu/static/models/Recon3D.json` | `aba925f17547a42f9fdb4c1f685d89364cbf4979bbe7862e9f793af7169b26d5` |
| CellMiner processed expression archive | `https://discover.nci.nih.gov/cellminer/download/processeddataset/nci60_RNA__Affy_HG_U133%28A_B%29_GCRMA.zip` | `334dca10d425eaf9fa3c7fc386f8d413ef32d8408eb4fc5508ca75e6cf6b734e` |
| Corrected consumption and release workbook, Jain et al. 2012 | Supplementary Database S1 of the article of record, DOI `10.1126/science.1218595`. Retrieve through an ordinary browser session; no programmatic route is available | `9a089755b35e16fd17d61b50865ca695bbb33ec8bc8616b285cd5e1fdbfa2ce4` |

## Independent evaluation inputs, retrieved by the pipeline

Both are public repositories and are cloned at a pinned commit by `copeland/copeland_prepare.py`.

| Input | Route | Commit |
|---|---|---|
| Research compendium, Copeland et al. 2023 | `https://github.com/oldhamlab/Copeland.2023.hypoxia.flux` | `354bd99459367f373c7b8a641bca4667bf96c779` |
| Companion expression data package | `https://github.com/wmoldham/rnaseq.lf.hypoxia.molidustat` | `1b6563bac8ba50019e90f7e6c6ed7d2f49a90ea5` |

The article of record is Copeland et al. 2023, eLife 12:e82597, DOI `10.7554/eLife.82597`. The
underlying sequencing is deposited under BioProject `PRJNA721596`.

## Culture medium formulations

The historical environment follows a representative current formulation of RPMI-1640; the
independent evaluation follows supplier formulation tables for MCDB 131. Exact lot information
is not published for either original experiment. Each medium module records its sources and
labels every component by the evidence behind its availability, separately from the assumed
uptake capacity given to it.

## Derived outputs in this repository

The larger derived outputs live under `data/`: per-arm predictions and feasible ranges for all
three reserved scenarios, the 336,336-row agreement ledger, the nested-range ledgers for the
development pilot and the reserved cohort, the cross-encoding regret, distance and self-check
tables, the endpoint reconstruction ledgers, the independent-evaluation predictions and crossed
control, and the reaction-evidence matrix. They are released under CC-BY-4.0 together with the
result tables, and are covered by the same archived record as the code.

## Archived record

Release v1.1.0, the version that accompanies the manuscript revision, is archived at
https://doi.org/10.5281/zenodo.22790149; release v1.0.0 is archived at
https://doi.org/10.5281/zenodo.22745339. The version-independent record covering all releases
is https://doi.org/10.5281/zenodo.22745338. Each identifier is minted when its release is
archived, so the archived copy of this file names the identifiers of the releases before it.
