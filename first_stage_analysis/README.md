# first_stage_analysis

1次解析．`data/` の実データを集計・換算・描画しただけの表と図，その生成コードと報告を置く．モデル・推定・復号は含まない（それらは `second_stage_analysis/`，非公開）．

## 規則

- `data/` と同じく号機ごとのフォルダ `01_kashiwa/`，`02_sakura/`，`03_yomogi/`，`04_botan/` に分け，その中をテーマ `<topic>/` で分ける．号機横断のもの（比較図，号機別集計，生成コード，報告）は `all/<topic>/` に置く．
- ファイル名は `<topic>_<内容>.<拡張子>`．号機フォルダの中では号機名を付けない（`01_kashiwa/hk/hk_daily.csv`）．
- スクリプトは `all/<topic>/<topic>_<内容>.py`．入力は `data/`（と1次解析の表）だけ．号機別の出力は `<sat>/<topic>/` に，号機横断の出力は自分のフォルダに書く．
- テーマをまたぐ報告は `all/operations.md`．テーマ固有の報告は `all/<topic>/<topic>_<内容>.md`．

## 内容

### `<sat>/`（号機ごと）

| topic | ファイル | 内容 |
|---|---|---|
| `hk/` | `hk_daily.csv`，`hk_daily.png`，`hk_surface.png` | 生HKの日別中央値（バッテリー，6面温度） |
| `hk/` | `hk_orbit.png` | 放出後7日以降で受信欠落のない最長区間の先頭から最大5周回分の生HK（`all/hk/hk_orbit.py` が生成） |
| `orbit/` | `orbit.csv` | TLEの軌道要素から算出した遠地点・近地点高度 |
| `aprs/` | `aprs_stations.csv` | 局ごとのAPRSパケット数と区分（`all/aprs/aprs_analysis.py` が生成） |
| `images/` | `images_<instrument>[_<mission>].jpg` | 軌道上画像のサムネイル一覧．複数ミッションのカメラはミッション別（`all/images/images_sheet.py` が生成） |
| `downlink/` | `downlink_daily.csv`，`downlink_daily.png` | 地上局の日次受信件数（BOTANのみ） |

### `all/`（号機横断）

| topic | ファイル | 内容 |
|---|---|---|
| — | `operations.md` | 運用報告（運用期間，軌道，HK，地上局，APRS，画像） |
| `hk/` | `hk_compare.png`，`hk_surface_compare.png` | 日次中央値の号機間比較 |
| `hk/` | `hk_orbit.py`，`hk_orbit_windows.csv` | 数周回プロットの生成コードと採用区間 |
| `orbit/` | `orbit_perigee.png` | 近地点高度の号機間比較 |
| `aprs/` | `aprs_analysis.py`，`aprs_summary.csv`，`aprs_activity.png` | APRS集計の生成コード，号機別の中継実績 |
| `images/` | `images_sheet.py`，`images_capture.md` | 画像一覧の生成コード，撮影時刻・撮影条件の導出 |

日次HK，軌道高度，日次受信件数の表と図は運用時の内部ツールで作成したもので，ここに生成コードはない．

`hk_orbit.png` の食の帯（5面の太陽電池電圧がすべて 2 V 未満）と周回の目盛（TLEの平均運動）は表示上の目安であり，解析値ではない．画像一覧では，部分受信で末尾の欠けた画像は復号できた範囲まで表示し，復号できないものは灰色枠に `images.csv` の status を添える．

## 実行

```bash
python3 -m venv --system-site-packages venv && venv/bin/pip install -r requirements.txt
venv/bin/python first_stage_analysis/all/hk/hk_orbit.py
venv/bin/python first_stage_analysis/all/images/images_sheet.py
venv/bin/python first_stage_analysis/all/aprs/aprs_analysis.py
```
