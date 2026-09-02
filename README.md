# opendata

Public satellite-operations data for KASHIWA, SAKURA, YOMOGI, and BOTAN.

- [Operations report](report/operations/report.md)
- [Image capture metadata (撮影時刻・撮影条件)](report/operations/image_capture.md)
- [Attitude analysis (沿磁力線制御の成立性)](report/attitude/attitude.md)

## Layout

データは由来で分けて管理する：衛星からダウンリンクしたもの（`data/downlink/`），地上で記録したもの（`data/recorded/`），それらから解析で生み出した2次データ（`data/derived/`），そして解析結果の報告（`report/`）．

- `data/downlink/<sat>/`: 衛星からダウンリンクしたデータ．
  - `hk/hk.csv`: 全期間の生HK（90秒周期，`time_utc` 付き）．KASHIWAは41バイト16進ログの復号値で，`hshk.csv`（高速HK，6秒周期）も持つ．SAKURAは `hshk_<日付>.csv` を持つ．BOTANの `hshk.csv` はアンテナ展開後の高速HKの未復号16進フレーム．
  - `aprs/aprs.csv`: APRSデジピートミッションで地上局が受信したパケット（date, time_jst, from, to, via, info, source_log）．`aprs/msg.csv` はMSGミッション（衛星のメッセージ蓄積メモリ）のダウンリンク16進フレーム（SAKURA・YOMOGI・BOTAN）．BOTANは `aprs/missionlog.csv`（ミッションログの16進DL）も持つ．
  - `mog/mog.csv`: KASHIWAのMoGミッションのダウンリンク16進フレーム．`mog_20240618.mp3` は受信音声．
  - `gyro/gyro.csv`: BOTANのジャイロミッションのダウンリンク16進フレーム（未復号）．
  - 16進系CSVは共通の列（downlink_date, time_jst, kind, hex, source_log．kind: cmd=コマンド行, frame=受信フレーム）を持つ．
  - `images/<instrument>/`: 復元済みの on-orbit 画像．Images are grouped **by instrument (camera unit)**, not by mission: `01_kashiwa/images/cam/`, `02_sakura/images/ecam/` + `scam/`, `03_yomogi/images/cam/`, `04_botan/images/cam/`. The mission / product tag (BOTAN: CORN, AURORA, PUMICE; YOMOGI: AFR, AKS) stays as the filename prefix and in the `mission` column of `images.csv`.
  - `images.csv`: その衛星の全画像のメタデータ．lists every image with size, brightness, status, SHA-256, capture-time range, exposure, and confidence. File names and the `date` column carry the recovered capture datetime (JST, `<prefix>_<yyyymmdd>T<hhmm>_<sha10>`); the original downlink / restoration date is kept in `downlink_date`. Internal source paths, command bytes, and processing notes are excluded. Derivation: [report/operations/image_capture.md](report/operations/image_capture.md).
- `data/recorded/<sat>/`: 地上で記録したデータ．TLE取得履歴（`tle.csv`：取得時刻，取得元（地上局DB・運用ログ・n2yo），生のTLE 2行）．β角計算の軌道要素はここから取り出す．
- `data/derived/<sat>/`: 解析で生み出した2次データ（表）．日次HK（`hk_daily.csv`：生HKの日別中央値），軌道高度（`orbit.csv`：TLEから算出），姿勢解析の日別指標（`beta_daily.csv`, `attitude_daily.csv`, `eclipse_events.csv`），BOTANのみ地上局の日次受信件数（`downlink_daily.csv`）とMSGミッションのデコード済みメッセージ（`aprs_messages.csv`）．号機間比較の表は `data/derived/00_compare/`．
- `analysis/`: 2次データを再現する解析スクリプト，one sub-folder per topic. See [analysis/README.md](analysis/README.md).
  - `attitude/`: `attitude_analysis.py`（日次データ）と `attitude_transient.py`（生HK）
- `report/`: 解析結果の報告と図．
  - `operations/`: operation periods, orbit, HK, ground station, images
  - `attitude/`: magnetic-alignment analysis, beta angle, settling, spin, alpha/epsilon

Internal source documents, decoder tools, private paths, and red-grid intermediate images are excluded.
