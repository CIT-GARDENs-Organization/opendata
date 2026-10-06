"""数周回分の生HK（90 s 周期）をそのままプロットする．

1次データだけを入力とし，解析（2次データ）には依存しない．

入力
  data/downlink/<sat>/hk/hk.csv   全期間の生HK（time_utc 付き）
  data/recorded/<sat>/tle.csv     TLE取得履歴（軌道周期の目盛に平均運動を使う）
出力
  first_stage_analysis/<sat>/hk/hk_orbit.png      衛星ごとの数周回分の生HK（太陽電池電圧，外面温度，バッテリー電圧・電流・温度）
  first_stage_analysis/all/hk/hk_orbit_windows.csv    採用した区間（開始・終了，件数，周回数，食の回数と長さ，各量の範囲）

区間の選び方
  放出後 7 日以降で，受信欠落（サンプル間隔 > 300 s）を含まない最長の連続区間を選び，
  先頭から最大 5 周回に切り詰める．軌道周期は，区間開始に最も近いエポックの TLE の平均運動から求める．
  食は 5 面の太陽電池電圧がすべて 2 V 未満の区間とし，灰色の帯で示す（HKの値だけで判定した表示上の目安）．
"""
from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
DL = ROOT / "data" / "downlink"
RECORDED = ROOT / "data" / "recorded"
FIRST = ROOT / "first_stage_analysis"  # 号機別の図は first_stage_analysis/<sat>/hk/
FIG = Path(__file__).resolve().parent  # 号機横断の表はこのフォルダ（all/hk/）

GAP_S = 300.0          # これを超える間隔を受信欠落とみなす
MIN_DAY = 7.0          # 放出後この日数以降から選ぶ
MAX_ORBITS = 5
V_DARK = 2.0           # 全面がこの電圧未満なら食とみなす

# 放出時刻（ISS放出．BOTANは最初のHK時刻）とカタログ番号
SATS = {
    "01_kashiwa": dict(name="KASHIWA", deploy="2024-04-11T10:41:00Z", norad=59508),
    "02_sakura": dict(name="SAKURA", deploy="2024-08-29T09:45:00Z", norad=60954),
    "03_yomogi": dict(name="YOMOGI", deploy="2024-12-09T11:15:00Z", norad=62298),
    "04_botan": dict(name="BOTAN", deploy="2025-10-10T09:50:58Z", norad=65942),
}

# 面の配色: 軸ごとに 1 色（X 青，Y 橙，Z 緑），+面は実線，-面は破線
HUE = {"x": "#2a78d6", "y": "#eb6834", "z": "#1baf7a"}
STYLE = {
    "xp": ("+X", HUE["x"], "-"), "xm": ("-X", HUE["x"], "--"),
    "yp": ("+Y", HUE["y"], "-"), "ym": ("-Y", HUE["y"], "--"),
    "zp": ("+Z", HUE["z"], "-"), "zm": ("-Z", HUE["z"], "--"),
}
INK = "#0b0b0b"
INK2 = "#52514e"
ECLIPSE = "#9a9a96"


# --------------------------------------------------------------------------- load
def load_hk(sat: str) -> pd.DataFrame:
    """号機ごとに列名が違う生HKを統一列名に揃える．
    t (UTC), T_xp..T_zm, T_bpb, T_bat, V_xm..V_zm, bat_v, bat_i, days, dark"""
    p = DL / sat / "hk" / "hk.csv"
    if sat == "01_kashiwa":
        # time_utc は day カウンタと起動エポックから復元済み（data/downlink/01_kashiwa/hk/hk.csv）
        d = pd.read_csv(p)
        d = d[(d.day < 150)].copy()
        out = pd.DataFrame({
            "t": pd.to_datetime(d.time_utc, utc=True), "T_xp": d.x_plus, "T_xm": d.x_minus, "T_yp": d.y_plus, "T_ym": d.y_minus,
            "T_zp": d.z_plus, "T_zm": d.z_minus, "T_bpb": d.bpb, "T_bat": d.bat_temp,
            "V_xm": d.v_xm, "V_yp": d.v_yp, "V_ym": d.v_ym, "V_zp": d.v_zp, "V_zm": d.v_zm, "bat_v": d.bat_v, "bat_i": d.bat_i})
    elif sat == "02_sakura":
        d = pd.read_csv(p)
        t = pd.to_datetime(d.time_utc, utc=True, errors="coerce")
        out = pd.DataFrame({"t": t, "T_xp": d["+X_HK"], "T_xm": d["-X_HK"], "T_yp": d["+Y_HK"], "T_ym": d["-Y_HK"],
                            "T_zp": d["+Z_HK"], "T_zm": d["-Z_HK"], "T_bpb": d["BPB_HK"], "T_bat": d["BAT_HK"],
                            "V_xm": d["-X_V"], "V_yp": d["+Y_V"], "V_ym": d["-Y_V"], "V_zp": d["+Z_V"], "V_zm": d["-Z_V"],
                            "bat_v": d["BAT_V"], "bat_i": d["BAT_I"]})
    elif sat == "03_yomogi":
        d = pd.read_csv(p, low_memory=False)
        t = pd.to_datetime(d.time_utc, utc=True, errors="coerce")
        out = pd.DataFrame({"t": t, "T_xp": d["'+X_Temp"], "T_xm": d["'-X_Temp"], "T_yp": d["'+Y_Temp"],
                            "T_ym": d["'-Y_Temp"], "T_zp": d["'+Z_Temp"], "T_zm": d["'-Z_Temp"], "T_bpb": d["BPB_Temp"],
                            "T_bat": d["BAT_Temp"], "V_xm": d["'-X_V"], "V_yp": d["'+Y_V"], "V_ym": d["'-Y_V"],
                            "V_zp": d["'+Z_V"], "V_zm": d["'-Z_V"], "bat_v": d["BAT_V"], "bat_i": d["BAT_I"]})
    else:
        d = pd.read_csv(p, low_memory=False)
        t = pd.to_datetime(d.time_utc, utc=True, errors="coerce")
        out = pd.DataFrame({"t": t, "T_xp": d["'+X_Temp [℃]"], "T_xm": d["'-X_Temp [℃]"],
                            "T_yp": d["'+Y_Temp [℃]"], "T_ym": d["'-Y_Temp [℃]"], "T_zp": d["'+Z_Temp [℃]"],
                            "T_zm": d["'-Z_Temp [℃]"], "T_bpb": d["BPB_Temp [℃]"], "T_bat": d["BAT_Temp [℃]"],
                            "V_xm": d["'-X_V [V]"], "V_yp": d["'+Y_V [V]"], "V_ym": d["'-Y_V [V]"], "V_zp": d["'+Z_V [V]"],
                            "V_zm": d["'-Z_V [V]"], "bat_v": d["BAT_V [V]"], "bat_i": d["BAT_I [A]"]})
    for c in out.columns:
        if c != "t":
            out[c] = pd.to_numeric(out[c], errors="coerce")
    out = out.dropna(subset=["t"]).sort_values("t").drop_duplicates("t").reset_index(drop=True)
    deploy = pd.Timestamp(SATS[sat]["deploy"].replace("Z", "+00:00"))
    out["days"] = (out.t - deploy).dt.total_seconds() / 86400.0
    out["dark"] = out[["V_xm", "V_yp", "V_ym", "V_zp", "V_zm"]].max(axis=1) < V_DARK
    return out


def tle_epoch(line1: str) -> datetime:
    yy = int(line1[18:20])
    doy = float(line1[20:32])
    return datetime(2000 + yy, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy - 1)


def orbit_period_s(sat: str, when: pd.Timestamp) -> float:
    """区間開始に最も近いエポックの TLE（該当カタログ番号）の平均運動 [rev/day] から周期 [s]．"""
    d = pd.read_csv(RECORDED / sat / "tle.csv", dtype=str)
    d = d[d.line2.str[2:7].astype(int) == SATS[sat]["norad"]].drop_duplicates("line1")
    epoch = d.line1.map(tle_epoch)
    k = (epoch - when.to_pydatetime()).abs().idxmin()
    n_rev_day = float(d.loc[k, "line2"][52:63])
    return 86400.0 / n_rev_day


def pick_window(hk: pd.DataFrame) -> pd.DataFrame:
    d = hk[hk.days >= MIN_DAY].reset_index(drop=True)
    t = d.t.values
    dt = np.diff(t).astype("timedelta64[s]").astype(float)
    brk = np.where(dt > GAP_S)[0]
    starts = np.r_[0, brk + 1]
    ends = np.r_[brk, len(t) - 1]
    dur = (t[ends] - t[starts]).astype("timedelta64[s]").astype(float)
    k = int(np.argmax(dur))
    return d.iloc[starts[k]:ends[k] + 1].copy().reset_index(drop=True)


def eclipse_runs(w: pd.DataFrame) -> list[tuple[float, float]]:
    runs, start = [], None
    for h, dark in zip(w.h, w.dark):
        if dark and start is None:
            start = h
        elif not dark and start is not None:
            runs.append((start, h))
            start = None
    if start is not None:
        runs.append((start, float(w.h.iloc[-1])))
    return runs


def main():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "grid.color": "#b8b7b2",
        "figure.dpi": 130, "axes.edgecolor": "#b8b7b2", "axes.labelcolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False,
    })
    rows = []
    for sat, meta in SATS.items():
        name = meta["name"]
        w = pick_window(load_hk(sat))
        t0 = w.t.iloc[0]
        period = orbit_period_s(sat, t0)
        w = w[(w.t - t0).dt.total_seconds() <= MAX_ORBITS * period].copy()
        w["h"] = (w.t - t0).dt.total_seconds() / 3600.0
        runs = eclipse_runs(w)
        n_orb = (w.h.iloc[-1] * 3600) / period
        t1 = w.t.iloc[-1]
        dt = np.diff(w.t.values).astype("timedelta64[s]").astype(float)
        rows.append(dict(
            sat=name, start_utc=t0.strftime("%Y-%m-%dT%H:%M:%SZ"), end_utc=t1.strftime("%Y-%m-%dT%H:%M:%SZ"),
            day=round(float(w.days.iloc[0]), 1), n=len(w), median_dt_s=float(np.median(dt)), max_dt_s=float(dt.max()),
            period_min=round(period / 60, 1), orbits=round(float(n_orb), 2), n_eclipse=len(runs),
            eclipse_min_median=round(float(np.median([(b - a) * 60 for a, b in runs])), 1) if runs else np.nan,
            bat_v_min=float(w.bat_v.min()), bat_v_max=float(w.bat_v.max()),
            bat_i_min=float(w.bat_i.min()), bat_i_max=float(w.bat_i.max()),
            T_bat_min=float(w.T_bat.min()), T_bat_max=float(w.T_bat.max()),
        ))
        print(f"{name}: {t0} – {t1}  day {w.days.iloc[0]:.1f}  n={len(w)}  {n_orb:.2f} orbits  eclipses={len(runs)}")

        fig, axes = plt.subplots(5, 1, figsize=(11, 11), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.1, 1.4, 1, 1, 1], hspace=0.28))
        fig.subplots_adjust(top=0.94, bottom=0.05, left=0.08, right=0.99)
        for ax in axes:
            for a, b in runs:
                ax.axvspan(a, b, color=ECLIPSE, alpha=0.18, lw=0)
        # 太陽電池電圧（+X はアンテナ面で太陽電池なし）
        ax = axes[0]
        for key in ("xm", "yp", "ym", "zp", "zm"):
            lab, col, ls = STYLE[key]
            ax.plot(w.h, w[f"V_{key}"], color=col, ls=ls, lw=1.2, label=lab)
        ax.set_ylabel("Panel V [V]")
        ax.legend(ncol=5, fontsize=8, loc="lower left", bbox_to_anchor=(0, 1.0), handlelength=2.5)
        # 外面温度
        ax = axes[1]
        for key in ("xp", "xm", "yp", "ym", "zp", "zm"):
            lab, col, ls = STYLE[key]
            ax.plot(w.h, w[f"T_{key}"], color=col, ls=ls, lw=1.2, label=lab)
        ax.plot(w.h, w.T_bpb, color=INK2, ls=":", lw=1.2, label="BPB")
        ax.set_ylabel("Panel T [°C]")
        ax.legend(ncol=7, fontsize=8, loc="lower left", bbox_to_anchor=(0, 1.0), handlelength=2.5)
        # バッテリー
        axes[2].plot(w.h, w.bat_v, color=INK, lw=1.2)
        axes[2].set_ylabel("Bat V [V]")
        axes[3].plot(w.h, w.bat_i, color=INK, lw=1.2)
        axes[3].axhline(0, color="#b8b7b2", lw=0.8)
        axes[3].set_ylabel("Bat I [A]")
        axes[4].plot(w.h, w.T_bat, color=INK, lw=1.2)
        axes[4].set_ylabel("Bat T [°C]")
        axes[4].set_xlabel(f"Time since {t0:%Y-%m-%d %H:%M} UTC [h]")
        axes[4].set_xlim(0, w.h.iloc[-1])
        # 周回の目盛
        for k in range(1, int(n_orb) + 1):
            for ax in axes:
                ax.axvline(k * period / 3600, color="#b8b7b2", lw=0.6, ls=(0, (2, 3)))
        fig.suptitle(f"{name} raw HK, day {w.days.iloc[0]:.0f}, {n_orb:.1f} orbits",
                     x=0.01, y=0.99, ha="left", fontsize=11, color=INK)
        fig.text(0.01, 0.965, f"grey: eclipse (all panels < {V_DARK:g} V)   dotted: orbit period {period/60:.1f} min (TLE)",
                 ha="left", fontsize=8, color=INK2)
        (FIRST / sat / "hk").mkdir(parents=True, exist_ok=True)
        fig.savefig(FIRST / sat / "hk" / "hk_orbit.png", bbox_inches="tight")
        plt.close(fig)

    pd.DataFrame(rows).to_csv(FIG / "hk_orbit_windows.csv", index=False)
    print("wrote", FIG / "hk_orbit_windows.csv")


if __name__ == "__main__":
    main()
