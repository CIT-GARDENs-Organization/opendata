# opendata

[![CC BY-NC-ND 4.0](https://licensebuttons.net/l/by-nc-nd/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

衛星 KASHIWA・SAKURA・YOMOGI・BOTAN の運用データと 1次解析．

## 報告

| 文書 | 内容 |
|---|---|
| [operations.md](first_stage_analysis/all/operations.md) | 運用期間・軌道・HK・地上局・APRS・画像 |
| [images_capture.md](first_stage_analysis/all/images/images_capture.md) | 画像の撮影時刻と撮影条件の導出 |

## 構成

| フォルダ | 内容 |
|---|---|
| `data/downlink/` | 衛星からダウンリンクしたデータ |
| `data/recorded/` | 地上で記録したデータ |
| `first_stage_analysis/` | 実データの集計・換算・描画による表と図，生成コード．詳細は [first_stage_analysis/README.md](first_stage_analysis/README.md) |

- 号機ごとのフォルダ: `01_kashiwa/`，`02_sakura/`，`03_yomogi/`，`04_botan/`
- 号機フォルダの下: テーマごとのフォルダ `<topic>/`
- ファイル名: `<topic>_<内容>.<拡張子>`

## data/downlink/

| ファイル | 号機 | 内容 |
|---|---|---|
| `hk/hk.csv` | 全号機 | 全期間の HK．90秒周期，時刻は `time_utc` |
| `hk/hshk.csv` | KASHIWA | 高速 HK．6秒周期 |
| `hk/hshk_<日付>.csv` | SAKURA | 高速 HK |
| `hk/hshk.csv` | BOTAN | アンテナ展開後の高速 HK．未復号の 16進フレーム |
| `aprs/aprs.csv` | 全号機 | APRS デジピートミッションで地上局が受信したパケット |
| `aprs/msg.csv` | SAKURA・YOMOGI・BOTAN | MSG ミッションのメッセージ蓄積メモリ．16進フレーム |
| `aprs/missionlog.csv` | BOTAN | ミッションログ．16進フレーム |
| `mog/mog.csv` | KASHIWA | MoG ミッション．16進フレーム |
| `mog/mog_20240618.mp3` | KASHIWA | MoG ミッションの受信音声 |
| `gyro/gyro.csv` | BOTAN | ジャイロミッションの IMU 姿勢．未復号の 16進フレーム |
| `images/<instrument>/` | 全号機 | 復元済みの軌道上画像．カメラごとのフォルダ |
| `images.csv` | 全号機 | 全画像のメタデータ |

### 16進フレームの列

| 列 | 意味 |
|---|---|
| `downlink_date` | ダウンリンク日 |
| `time_jst` | 時刻（JST） |
| `kind` | `cmd`: コマンド行，`frame`: 受信フレーム |
| `hex` | 16進のバイト列 |
| `source_log` | 元の運用ログ名 |

### 画像

| 号機 | カメラ | ファイル名の接頭辞 |
|---|---|---|
| KASHIWA | `cam/` | `kashiwa` |
| SAKURA | `ecam/`，`scam/` | `sakura_ecam`，`sakura_scam` |
| YOMOGI | `cam/` | ミッション名 `afr`，`aks` |
| BOTAN | `cam/` | ミッション名 `aurora`，`corn`，`pumice` |

- ファイル名: `<接頭辞>_<yyyymmdd>T<hhmm>_<SHA-256 の先頭10桁>.jpg`．時刻は推定した撮影時刻（JST）
- 撮影時刻が不明な画像のファイル名: `<接頭辞>_unknown_<SHA-256 の先頭10桁>.jpg`
- `downlink_date` 列: ダウンリンク日・復元日
- `images.csv` の列: 画素数，バイト数，輝度，コントラスト，状態，撮影時刻の範囲，ISO 感度，シャッター速度，撮影時刻の確度，SHA-256
- 撮影時刻の導出: [images_capture.md](first_stage_analysis/all/images/images_capture.md)

## data/recorded/

| ファイル | 内容 |
|---|---|
| `tle.csv` | TLE の取得履歴．取得時刻，取得元（地上局 DB・運用ログ・n2yo），TLE 2行 |

## 公開の範囲

含めないもの:

- 内部の設計文書と復号ツール
- 運用者の端末のパス
- 画像復元の中間画像
- モデル・推定・復号を伴う 2次解析

## 実行

```bash
python3 -m venv --system-site-packages venv && venv/bin/pip install -r requirements.txt
venv/bin/python first_stage_analysis/all/hk/hk_orbit.py
```

他のスクリプトは [first_stage_analysis/README.md](first_stage_analysis/README.md#実行)．

## ライセンス

[CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/)．全文は [LICENSE](LICENSE)．

利用時のクレジット表記: 「千葉工業大学」または「ChibaTech」．
