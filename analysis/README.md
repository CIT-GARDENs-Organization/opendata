# analysis

再現可能な解析をテーマごとのサブフォルダに置く．スクリプトのみを置き，出力の表は `data/derived/` に，図は `report/<テーマ>/` に書き出す．

| フォルダ | 内容 | 入力 | 出力 |
|---|---|---|---|
| `attitude/` | 沿磁力線制御の成立性，β角，姿勢の整定，スピン，α/ε | `data/derived/<sat>/{hk_daily,orbit}.csv`，`data/recorded/<sat>/tle.csv`，`data/downlink/<sat>/hk/`（生HK） | `data/derived/<sat>/`，`data/derived/00_compare/`，`report/attitude/*.png` |
| `aprs/` | APRSデジピータ・MSGミッションの運用集計 | `data/downlink/<sat>/aprs/{aprs,msg}.csv` | `data/derived/<sat>/aprs_stations.csv`，`data/derived/00_compare/aprs_summary.csv`，`report/operations/aprs_activity.png` |
| `gyro/` | BOTANジャイロのクォータニオン復元とスケール検証 | `data/downlink/04_botan/gyro/gyro.csv` | `data/derived/04_botan/gyro_attitude.csv`，`report/gyro/*.png` |

## 実行

```
python3 -m venv --system-site-packages venv && venv/bin/pip install -r analysis/requirements.txt
venv/bin/python analysis/attitude/attitude_analysis.py    # 日次データの解析
venv/bin/python analysis/attitude/attitude_transient.py   # 生HK（data/downlink/<sat>/hk/）の解析
venv/bin/python analysis/aprs/aprs_analysis.py            # APRS・MSGミッションの集計
venv/bin/python analysis/gyro/gyro_analysis.py            # BOTANジャイロのクォータニオン復元
venv/bin/python analysis/gyro/gyro_animation.py          # 姿勢キューブのアニメーション（要ffmpeg）
```
