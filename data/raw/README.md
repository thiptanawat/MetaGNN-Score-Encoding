# Third-party raw inputs (not redistributed)

The locked study reads four third-party files from this directory. They are not part of the
repository because their terms of use differ from the repository's licenses; each must be retrieved
once from its source and placed at the exact path below. `src/verify_inputs.py` checks the SHA-256
of every file before anything is built, and `reproduce.sh` (any stage) runs that check first, so a
wrong, updated or HTML-substituted download is reported before it can be read. The digests are
the ones recorded in `protocol/PROTOCOL_LOCK.json` and `manifests/raw_sources.tsv`.

| Path under `data/raw/` | What it is | Route | SHA-256 |
|---|---|---|---|
| `Recon3D.json` | The pinned Recon3D reconstruction from BiGG Models | `http://bigg.ucsd.edu/static/models/Recon3D.json` (direct download) | `aba925f17547a42f9fdb4c1f685d89364cbf4979bbe7862e9f793af7169b26d5` |
| `cellminer.zip` | The processed NCI-60 Affymetrix HG-U133(A-B) GCRMA archive from CellMiner | `https://discover.nci.nih.gov/cellminer/download/processeddataset/nci60_RNA__Affy_HG_U133%28A_B%29_GCRMA.zip` (direct download) | `334dca10d425eaf9fa3c7fc386f8d413ef32d8408eb4fc5508ca75e6cf6b734e` |
| `core/NIHMS419088-supplement-Database_S1.xls` | The corrected consumption-and-release workbook of Jain et al. (2012), Supplementary Database S1 | PubMed Central record PMC3526189 (`https://pmc.ncbi.nlm.nih.gov/articles/PMC3526189/`), supplementary materials; the article of record is DOI 10.1126/science.1218595 | `9a089755b35e16fd17d61b50865ca695bbb33ec8bc8616b285cd5e1fdbfa2ce4` |
| `core/NIHMS419088-supplement.pdf` | The supplementary methods of the same article, whose Table S1 supplies the assay chemistry that `src/prepare_data.py` extracts with `pdftotext` | The same PMC3526189 record | `1bf01aa2715af80924aa3dbc74c00507259f9517877f23e3157344966cc0f187` |

Notes on retrieval.

- The two PMC files are served to browsers; programmatic requests to the PMC file URLs return an
  HTML access page rather than the file. Download the attachments from the record page in a
  browser, then check the digest (`shasum -a 256 <file>` on macOS, `sha256sum <file>` on Linux).
  A file of a few kilobytes that begins with `<!DOCTYPE html>` is an access page, not the input.
- `pdftotext` (Poppler) must be on `PATH` for the `data` stage of `reproduce.sh`; it is not needed
  to run the analyses on the deposited derived data.
- The two direct downloads were re-retrieved on 16 September 2026 and matched the recorded digests.

Once the four files are in place:

```
python -m src.verify_inputs --root . --check raw        # the four digests
python -m src.verify_inputs --root . --check preflight  # raw inputs, reference table and environment
python -m src.verify_inputs --root . --check lock       # the complete protocol lock (see verify/check_release.py)
```
