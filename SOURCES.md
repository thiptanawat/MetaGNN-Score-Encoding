# Sources

## Third-party inputs, not redistributed

The locked study reads four third-party files. Each must be retrieved once and placed under
`data/raw/` at the path given; `data/raw/README.md` describes the routes in more detail and
`src/verify_inputs.py` checks every digest before anything is built. The digests are the ones in
`protocol/PROTOCOL_LOCK.json` and `manifests/raw_sources.tsv`.

| Input | Path | Route | SHA-256 |
|---|---|---|---|
| Recon3D reconstruction (BiGG Models) | `data/raw/Recon3D.json` | `http://bigg.ucsd.edu/static/models/Recon3D.json` | `aba925f17547a42f9fdb4c1f685d89364cbf4979bbe7862e9f793af7169b26d5` |
| CellMiner processed NCI-60 expression archive, Affymetrix HG-U133(A-B) GCRMA | `data/raw/cellminer.zip` | `https://discover.nci.nih.gov/cellminer/download/processeddataset/nci60_RNA__Affy_HG_U133%28A_B%29_GCRMA.zip` | `334dca10d425eaf9fa3c7fc386f8d413ef32d8408eb4fc5508ca75e6cf6b734e` |
| Corrected consumption and release workbook, Jain et al. 2012, Supplementary Database S1 | `data/raw/core/NIHMS419088-supplement-Database_S1.xls` | PubMed Central record PMC3526189, supplementary materials; article of record DOI `10.1126/science.1218595`. The files are served to browsers, not to programmatic requests | `9a089755b35e16fd17d61b50865ca695bbb33ec8bc8616b285cd5e1fdbfa2ce4` |
| Supplementary methods of the same article, whose Table S1 supplies the assay chemistry (`src/prepare_data.py` extracts it with `pdftotext`, which must be installed) | `data/raw/core/NIHMS419088-supplement.pdf` | The same PMC3526189 record | `1bf01aa2715af80924aa3dbc74c00507259f9517877f23e3157344966cc0f187` |

The two direct downloads were re-retrieved on 16 September 2026 and matched the recorded digests.
The endpoint reconstruction (`endpoint/`) reads the workbook as well; the eleven-line list it
compares against is `endpoint/Source_Quality_Provenance.json`, which records the Zenodo archive
and member digest it was verified from.

## Independent evaluation inputs, retrieved by a script

The independent condition-response evaluation reads two public git repositories. They are not
redistributed; `copeland/fetch_inputs.sh` clones both, checks out the pinned commits, refuses to
continue if either checkout resolves to a different commit, and checks the digests of the two
files that `copeland/copeland_prepare.py` reads. The preparation step records the commit of each
checkout and the digest of the expression package in `results/copeland/prepared/prepare_summary.json`.

| Input | Route | Pinned commit | File read | SHA-256 |
|---|---|---|---|---|
| Research compendium, Copeland et al. 2023 | `https://github.com/oldhamlab/Copeland.2023.hypoxia.flux` | `354bd99459367f373c7b8a641bca4667bf96c779` | `data/growth_rates.rda` | `d051b317077f40e6e5d7047a07a32ceb6adeb87c07cd31e903949122ca4757b6` |
| Companion expression data package | `https://github.com/wmoldham/rnaseq.lf.hypoxia.molidustat` | `1b6563bac8ba50019e90f7e6c6ed7d2f49a90ea5` | `data/lf_hyp_bay_rnaseq.rda` | `4088ad11271e5a6def001018b8911fec8e0d7868462a1f390b86f9a146e4901f` |

The article of record is Copeland et al. 2023, eLife 12:e82597, DOI `10.7554/eLife.82597`. The
underlying sequencing is deposited under BioProject `PRJNA721596`. Reading the `.rda` files needs
the `rdata` package (`environment_extensions.txt`).

## Culture medium formulations

The historical environment follows a representative current formulation of RPMI-1640; the
independent evaluation follows supplier formulation tables for MCDB 131. Exact lot information is
not published for either original experiment. Each medium module records its sources and labels
every component by the evidence behind its availability, separately from the assumed uptake
capacity given to it.

## Derived data in this repository

The derived inputs the lock names (`data/`), the completed run outputs (`results/development_v3/`,
`results/reserved_v3/`), and the outputs of the post-lock analyses under `results/` are released
under CC-BY-4.0 together with the result tables, and are covered by the same archived record as
the code. The per-arm records of the two study runs are release assets whose digests are in
`RELEASE_MANIFEST.md`.

## Archived record

Release v1.2.0 is archived at https://doi.org/10.5281/zenodo.22796608 (the Zenodo copy of the
tagged tree, `thiptanawat/MetaGNN-Score-Encoding-v1.2.0.zip`, 29,372,581 bytes, MD5
`49302287a4ee49e71f4aa8da9a7234d6`, holds the same 357 files as the tag). Release v1.1.0 is
archived at https://doi.org/10.5281/zenodo.22790149 and release v1.0.0 at
https://doi.org/10.5281/zenodo.22745339. The version-independent record covering all releases is
https://doi.org/10.5281/zenodo.22745338. Each identifier is minted when its release is archived,
so the archived copy of this file names the identifiers of the releases before it; the repository's
`main` branch records the identifier of the latest release in the commit that follows its tag.
