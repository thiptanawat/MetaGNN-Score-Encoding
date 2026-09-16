# Clean-checkout verification, release v1.2.0

Performed 16 September 2026 on a machine that had never held the study: a fresh clone of commit
`71ffc4458f5b1a9372419f18541ff83e2d23d606` (the tree of release v1.2.0 without this log and the
checksum line it adds), a Python 3.12.3 virtual environment created from `environment.lock` on
Linux x86-64 (2 CPU threads), the four third-party raw inputs retrieved as `data/raw/README.md`
describes (the two direct downloads matched their digests; the two PubMed Central attachments
were taken from the authors' copies after a browser download and matched theirs), and the two
release-asset archives. The workstation the study ran on is macOS arm64 (`reports/hardware.json`),
so this is also a cross-platform check. The complete console record (`clean_checkout.log`) is
kept with the manuscript's review package; the results are summarized here. Every step passed.

| Step | Command | Result |
|---|---|---|
| Environment | `pip list` against `environment.lock` | the installed packages match the lock exactly |
| Clone | `git clone`, `git checkout 71ffc44` | 356 tracked files |
| 1. Release integrity before the raw inputs | `verify/check_release.py` | 355 files listed in `CHECKSUMS.sha256`, 355 match; lock: 51 paths named, 47 match, 4 absent (the raw inputs), 0 differ |
| 2. Raw inputs placed | `src.verify_inputs --check raw`; `check_release.py --strict` | all four digests verified; lock: 51 of 51 match |
| 3. Locked driver | `reproduce.sh --stage preflight`; `--stage preflight --partition test`; `src.verify_inputs --check lock`, `--check prepared`, `--check controls` | preflight verified (raw inputs, reference table, environment digest `6bdb8429...`); the complete protocol lock verified; the data-build receipt verified against its 17 outputs; the production-controls record verified |
| 4. Tests | `reproduce.sh --stage tests`; `reproduce_extensions.sh tests` | 74 passed (locked suite); 33 + 30 + 11 unit tests and 9 pytest checks passed (mechanism, independent evaluation, coverage, cross-encoding ranges) |
| 5. Reported numbers | `verify/verify_manuscript_numbers.py --manuscript ... --supplement ...` | 284 checks, 0 problems (values and text) |
| 6. Release assets | `verify/fetch_release_assets.sh all`; `check_release.py --strict`; `src.verify_inputs --check complete` and `--check audit` on `results/reserved_v3`, `--check complete` on `results/development_v3` | both archive digests verified; 849 reserved and 201 development flux arrays match the digests in their records; the completeness and audit checks of the locked verifier passed on both deposited runs |
| 7. Bounded reproduction, evaluation | `src.evaluation --partition test --lock-file protocol/PROTOCOL_LOCK.json` for the three scenarios, from the deposited `predictions.tsv`, into a new directory | each of the three evaluation records reproduced the deposited one on all 5,469 shared leaves with a maximum numeric difference of 0, including the 2,000-draw origin bootstrap |
| 8. Bounded reproduction, coverage | `reproduce_extensions.sh coverage` into a new directory | the six scenario and arm-set rows agree with the deposited summary on every count |
| 9. Figures | `reproduce_extensions.sh figures` into a new directory | figures 1 to 3 (from the reserved evaluation) and 4 to 6 (from the post-lock outputs) rendered; the pixel dimensions differ from the deposited files by up to 8 percent because the tight bounding box depends on the installed font metrics and matplotlib build, the content is the same |
| 10. Bounded reproduction of the locked pipeline | `reproduce.sh --stage smoke --workers 2 --output-dir <new>` in a second clone | data receipt and panel verified, 74 tests passed, the eight production controls passed (357 s), 21 optimizations (one development profile, six encodings, three scenarios and the three uniform-cost references) completed with status optimal, the saved-vector audit passed on all 21 arms, and the run summary was written |

Comparison of the 21 smoke arms with the deposited development arms of the same profile,
encoding and scenario (the arm tokens coincide because the input fingerprint is the same):

- the 96 reported exchange points agree within 7.4e-9 model units in every arm (20 arms within
  1e-9), and the 96 reported ranges within 3.1e-8; the classification threshold is 1e-5;
- the primary objective optima agree within 1.3e-11 and the secondary optima within 4.0e-8;
- the full 10,600-reaction flux vectors are not byte-identical on this platform: 4 to 29
  reactions per arm differ by more than 1e-6, by up to 10 units. These are alternative optima of
  the same linear programs, the degeneracy the study's readout is designed to project out;
  `reports/fresh_arm_reproduction_v3/comparison.json` records that on the original machine a fresh
  arm matched its saved record byte for byte. Bitwise portability of the interior vectors across
  operating systems and solver builds is not claimed; the reported exchange coordinates and their
  classification are what the study reports, and those reproduced.

Earlier the same day, on the workstation with the study's own environment (Python 3.12.14), the
assembled tree passed the same preflight, complete lock check, 74 locked tests, post-lock tests
and the 284-check number verification, and repeats of the coverage audit, the interval-position
analysis and the cross-encoding comparison from the deposited ledgers and per-arm records
reproduced the deposited outputs (`RELEASE_MANIFEST.md`, section 4).

Not performed: a complete rerun of the development or reserved study, the RIPTiDe run, the
second-solver check and the independent evaluation's optimizations on the new machine. Their
deposited outputs, run manifests and input digests are what the number checker and the
release checker verify.
