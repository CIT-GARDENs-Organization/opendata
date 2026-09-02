# satops

Public satellite-operations data for KASHIWA, SAKURA, YOMOGI, and BOTAN.

- [Operations report](report/operations/report.md)
- [Image capture metadata (撮影時刻・撮影条件)](report/operations/image_capture.md)
- [Attitude analysis (沿磁力線制御の成立性)](report/attitude/attitude.md)

## Layout

- `data/<sat>/`: curated inputs — daily HK (`hk_daily.csv`), orbit altitude from TLE (`orbit.csv`), orbital elements for beta angle (`tle_elements.csv`), on-orbit images (`images/<instrument>/`), and for BOTAN the daily downlink counts (`downlink_daily.csv`).
  - Images are grouped **by instrument (camera unit)**, not by mission: `01_kashiwa/images/cam/`, `02_sakura/images/ecam/` + `scam/`, `03_yomogi/images/cam/`, `04_botan/images/cam/`. The mission / product tag (BOTAN: CORN, AURORA, PUMICE; YOMOGI: AFR, AKS) stays as the filename prefix and in the `mission` column of `report/operations/images.csv`.
  - `report/operations/images.csv` lists every image with size, brightness, status, SHA-256, and the capture metadata recovered from the operation records: `capture_start` / `capture_end` (JST range), `time_basis`, `iso`, `shutter`, `exposure_mode`, `resolution_cmd`, `other_settings`, `capture_confidence`, `capture_source`, `capture_notes`. File names and the `date` column carry the recovered capture datetime (JST, `<prefix>_<yyyymmdd>T<hhmm>_<sha10>`); the original downlink / restoration date is kept in `downlink_date`. Derivation: [report/operations/image_capture.md](report/operations/image_capture.md).
  - Instruments: KASHIWA cam = B0286 / IMX219 (F2.1); SAKURA ecam = SC0023 / IMX219 (F2.0), scam = B0068 / OV5642 / Arducam M2516ZH01 (F2.0); YOMOGI cam = B0103 / IMX219 / Arducam M2506ZH04 (F2.0); BOTAN cam = B0103 / IMX219 / Arducam M40320M06S (F2.0)
- `analysis/`: reproducible analysis, one sub-folder per topic. See [analysis/README.md](analysis/README.md).
  - `attitude/`: scripts (`attitude_analysis.py`, `attitude_transient.py`) and their tables in `attitude/output/`
- `report/`: written reports with their figures.
  - `operations/`: operation periods, orbit, HK, ground station, images
  - `attitude/`: magnetic-alignment analysis, beta angle, settling, spin, alpha/epsilon

Internal source documents, raw HK workbooks, decoder tools, private paths, and red-grid intermediate images are excluded.
