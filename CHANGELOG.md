# Changelog

All notable changes to NC-MutScan. Format follows Keep a Changelog; versions are pinned by git tag and Zenodo DOI.

## [3.3.0] - RELEASE DATE (fill in at tagging)

Frozen release. Core build SHA-256 `b276d0d440453bf0c1c98647f8825b039a8775738b483fa4aad8cb1119edb2c4`
(short `b276d0d44045`); development set `dev_chr1_11.vcf` SHA-256 `7cf78f07c6701e8bf79eaa15a9f2d03745e52b3504d26ce64bff4a44afee29f7`.

### Added
- `GET /version` (via `release_app.py`; frozen core untouched): version, git tag/commit, pinned and live core hash, `build_verified`, DOI.
- `RELEASE.json` as the single source of release metadata.
- Live-API replay of the golden regression set and a blind subset (`tests/live_replay.py`), run by the `live-parity` CI job on every tag, weekly and on demand.
- External validation on SpliceVarDB (hg38 track export; 23,868 SNVs, calibration overlap n = 1 removed):
  AUROC 0.700, sensitivity 0.452, specificity 0.958, PPV 0.923, MCC 0.468. No parameter changed.
- `CITATION.cff` and `.zenodo.json` for archiving.

### Changed (relative to v3.1; confirm against your revision scripts before release)
- 9th TF PWM added (GATA3); TF motif decision threshold tau_TF = 4.2 bits (v3.1: <= -4.00 bits).

### Known limitations
- Sensitivity on SpliceVarDB is modest (0.452); location-stratified performance not yet analysed.
- Only SNVs evaluated in the external splice test.

## [3.1] - earlier release
- See project history.
