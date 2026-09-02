"""APRSデジピータ・メッセージミッションの運用解析．

KASHIWA/SAKURA/YOMOGI/BOTANはアマチュア無線のAPRSデジピータ（中継器）として運用された．
地上の各局が送信したパケットを衛星が中継し，地上局で受信・記録した．

入力
  data/downlink/<sat>/aprs/aprs.csv   受信したAPRSパケット（from, to, via, info）
  data/downlink/<sat>/aprs/msg.csv    MSGミッション（衛星のメッセージ蓄積メモリ）のDL 16進フレーム
  data/derived/04_botan/aprs_messages.csv  BOTANのデコード済みメッセージ

出力
  data/derived/<sat>/aprs_stations.csv   局ごとのパケット数と区分（CIT / 外部アマチュア局）
  data/derived/00_compare/aprs_summary.csv  号機別の中継実績
  report/operations/aprs_activity.png    図
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DL = ROOT / "data" / "downlink"
OUT = ROOT / "data" / "derived"
COMPARE = OUT / "00_compare"
FIG = ROOT / "report" / "operations"

SATS = {"01_kashiwa": "KASHIWA", "02_sakura": "SAKURA",
        "03_yomogi": "YOMOGI", "04_botan": "BOTAN"}
COLORS = {"KASHIWA": "#1f77b4", "SAKURA": "#d62728",
          "YOMOGI": "#2ca02c", "BOTAN": "#9467bd"}

# CITのクラブ局・衛星コールサイン（JS1Y**, JG6Y**）．これ以外を外部アマチュア局とみなす．
CIT = re.compile(r"^(JS1Y|JG6Y)")
# 埋め込みAPRSメッセージ CALL>DEST,PATH:text
EMB = re.compile(r"([A-Z0-9]{3,6}(?:-\d+)?)>([A-Z0-9-]{2,6})((?:,[A-Za-z0-9*-]+)*):([\x20-\x7e]{2,})")


def base_call(c: str) -> str:
    """SSID（-7 等）と * を除いた基幹コールサイン．"""
    return re.sub(r"[-*].*$", "", str(c)).strip().upper()


def analyze_packets(sat: str):
    p = DL / sat / "aprs" / "aprs.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p, dtype=str).fillna("")
    d["call"] = d["from"].map(base_call)
    d = d[d.call.str.len() >= 3]
    stations = []
    for call, g in d.groupby("call"):
        stations.append(dict(station=call, packets=len(g),
                             category="CIT" if CIT.match(call) else "amateur",
                             dates=g[g.date != ""].date.nunique()))
    st = pd.DataFrame(stations).sort_values(["category", "packets"],
                                            ascending=[True, False])
    (OUT / sat).mkdir(parents=True, exist_ok=True)
    st.to_csv(OUT / sat / "aprs_stations.csv", index=False)
    ext = st[st.category == "amateur"]
    return dict(
        sat=SATS[sat], packets=len(d),
        stations=st.station.nunique(),
        amateur_stations=ext.station.nunique(),
        amateur_packets=int(ext.packets.sum()),
        dates=d[d.date != ""].date.nunique(),
    )


def count_messages(sat: str) -> int:
    """MSGミッションのメモリに蓄積された固有メッセージ数（埋め込みAPRSメッセージの重複除外）．"""
    p = DL / sat / "aprs" / "msg.csv"
    if not p.exists():
        return 0
    d = pd.read_csv(p, dtype=str).fillna("")
    msgs = set()
    for h in d[d.kind == "frame"].hex:
        b = bytes(int(x, 16) for x in h.split())
        txt = "".join(chr(c) if 32 <= c < 127 else "\n" for c in b[16:])
        for m in EMB.finditer(txt):
            msgs.add((m.group(1), m.group(4).strip().rstrip("D")))  # 末尾のパディングDを除く
    return len(msgs)


def main():
    rows = []
    for sat in SATS:
        r = analyze_packets(sat)
        if r:
            r["msg_memory"] = count_messages(sat)
            rows.append(r)
            print(f"{r['sat']}: {r['packets']} pkts, {r['stations']} stations "
                  f"({r['amateur_stations']} amateur), {r['dates']} days, "
                  f"MSG memory {r['msg_memory']}")
    summ = pd.DataFrame(rows)
    COMPARE.mkdir(parents=True, exist_ok=True)
    summ.to_csv(COMPARE / "aprs_summary.csv", index=False)

    plt.rcParams.update({"font.size": 9, "axes.grid": True,
                        "grid.alpha": 0.3, "figure.dpi": 130})
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
    x = range(len(summ))
    cols = [COLORS[s] for s in summ.sat]
    cit = summ.packets - summ.amateur_packets
    ax[0].bar(x, cit, color=cols, alpha=0.45)
    ax[0].bar(x, summ.amateur_packets, bottom=cit, color=cols)
    ax[0].set_xticks(list(x))
    ax[0].set_xticklabels(summ.sat, rotation=20)
    ax[0].set_ylabel("Received APRS packets")
    ax[0].bar(0, 0, color="gray", alpha=0.45, label="CIT stations")
    ax[0].bar(0, 0, color="gray", label="External amateurs")
    ax[0].legend(fontsize=7, loc="upper left")
    ax[1].bar(x, summ.amateur_stations, color=cols)
    ax[1].set_xticks(list(x))
    ax[1].set_xticklabels(summ.sat, rotation=20)
    ax[1].set_ylabel("External amateur stations relayed")
    fig.tight_layout()
    fig.savefig(FIG / "aprs_activity.png")
    plt.close(fig)
    print("wrote", COMPARE / "aprs_summary.csv", "and", FIG / "aprs_activity.png")


if __name__ == "__main__":
    main()
