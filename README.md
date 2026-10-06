# opendata

Public satellite-operations data for KASHIWA, SAKURA, YOMOGI, and BOTAN.

- [Operations report (運用期間・軌道・HK・地上局・APRS・画像)](first_stage_analysis/all/operations.md)
- [Image capture metadata (撮影時刻・撮影条件の導出)](first_stage_analysis/all/images/images_capture.md)
- [Raw HK over a few orbits (生HKの数周回プロット)](first_stage_analysis/all/hk/hk_orbit.py)
- [On-orbit image sheets (軌道上画像の一覧)](first_stage_analysis/all/images/images_sheet.py)

## Layout

トップレベルは3つだけ．

| フォルダ | 内容 | 公開 |
|---|---|---|
| `data/` | 実データ．衛星からダウンリンクしたもの（`downlink/`）と地上で記録したもの（`recorded/`） | する |
| `first_stage_analysis/` | 1次解析．実データを集計・換算・描画しただけの表と図，その生成コードと報告 | する |
| `second_stage_analysis/` | 2次解析．モデル・推定・復号を伴う解析（姿勢，ジャイロ）．正しさの検証が済んでいない | しない（`.gitignore`） |

`data/` と `first_stage_analysis/` は号機ごとのフォルダ `01_kashiwa/`，`02_sakura/`，`03_yomogi/`，`04_botan/` に分け，その中をテーマ `<topic>/` で分ける．号機横断のもの（比較図，号機別集計，生成コード，報告）は `first_stage_analysis/all/` に置く．ファイル名は `<topic>_<内容>.<拡張子>`．スクリプトは `all/<topic>/` にあり，号機別の出力は各号機フォルダに，号機横断の出力は自分のフォルダに書く．

### data/

- `downlink/<sat>/`: 衛星からダウンリンクしたデータ．
  - `hk/hk.csv`: 全期間の生HK（90秒周期，`time_utc` 付き）．KASHIWAは41バイト16進ログの復号値で，`hshk.csv`（高速HK，6秒周期）も持つ．SAKURAは `hshk_<日付>.csv` を持つ．BOTANの `hshk.csv` はアンテナ展開後の高速HKの未復号16進フレーム．
  - `aprs/aprs.csv`: APRSデジピートミッションで地上局が受信したパケット（date, time_jst, from, to, via, info, source_log）．`aprs/msg.csv` はMSGミッション（衛星のメッセージ蓄積メモリ）のダウンリンク16進フレーム（SAKURA・YOMOGI・BOTAN）．BOTANは `aprs/missionlog.csv`（ミッションログの16進DL）も持つ．
  - `mog/mog.csv`: KASHIWAのMoGミッションのダウンリンク16進フレーム．`mog_20240618.mp3` は受信音声．
  - `gyro/gyro.csv`: BOTANのジャイロミッション（IMU姿勢）のダウンリンク16進フレーム（未復号）．
  - 16進系CSVは共通の列（downlink_date, time_jst, kind, hex, source_log．kind: cmd=コマンド行, frame=受信フレーム）を持つ．
  - `images/<instrument>/`: 復元済みの on-orbit 画像．Images are grouped **by instrument (camera unit)**, not by mission: `01_kashiwa/images/cam/`, `02_sakura/images/ecam/` + `scam/`, `03_yomogi/images/cam/`, `04_botan/images/cam/`. The mission / product tag (BOTAN: CORN, AURORA, PUMICE; YOMOGI: AFR, AKS) stays as the filename prefix and in the `mission` column of `images.csv`.
  - `images.csv`: その衛星の全画像のメタデータ．lists every image with size, brightness, status, SHA-256, capture-time range, exposure, and confidence. File names and the `date` column carry the recovered capture datetime (JST, `<prefix>_<yyyymmdd>T<hhmm>_<sha10>`); the original downlink / restoration date is kept in `downlink_date`. Internal source paths, command bytes, and processing notes are excluded. Derivation: [images_capture.md](first_stage_analysis/all/images/images_capture.md).
- `recorded/<sat>/`: 地上で記録したデータ．TLE取得履歴（`tle.csv`：取得時刻，取得元（地上局DB・運用ログ・n2yo），生のTLE 2行）．

### first_stage_analysis/

See [first_stage_analysis/README.md](first_stage_analysis/README.md).

- `<sat>/hk/`: `hk_daily.csv`，`hk_daily.png`，`hk_surface.png`（生HKの日別中央値），`hk_orbit.png`（数周回分の生HK）
- `<sat>/orbit/`: `orbit.csv`（TLEの軌道要素から算出した高度）
- `<sat>/aprs/`: `aprs_stations.csv`（局ごとのパケット数）
- `<sat>/images/`: `images_<instrument>[_<mission>].jpg`（画像のサムネイル一覧）
- `04_botan/downlink/`: `downlink_daily.csv`，`downlink_daily.png`（地上局の日次受信件数）
- `all/`: `operations.md`（運用報告），`hk/`（号機間比較図，`hk_orbit.py`，採用区間），`orbit/orbit_perigee.png`，`aprs/`（`aprs_analysis.py`，`aprs_summary.csv`，`aprs_activity.png`），`images/`（`images_sheet.py`，`images_capture.md`）

### second_stage_analysis/（非公開）

`attitude/`（沿磁力線制御の成立性，β角，α/ε），`gyro/`（クォータニオン復元），`aprs/`（MSGメッセージのデコード結果）．

Internal source documents, decoder tools, private paths, red-grid intermediate images, and unverified secondary analysis are excluded.
