"""BOTANジャイロミッション（BNO055系IMU）のクォータニオン解析．

ダウンリンクフレーム（data/downlink/04_botan/gyro/gyro.csv）の32バイトレコードを復号する．
レコード構成（ビッグエンディアン，運用チームの解析ツールのラベルに基づく）:
  [0]      マーカー 0xd0
  [1:3]    pakt_num（サンプル番号）
  [3:5]    time（ミッション開始からの秒）
  [5:11]   Accel_X/Y/Z（int16）
  [11:19]  Vector_i/j/k/Real（クォータニオン，int16）
  [19:29]  SAP電流 -X/+Y/-Y/+Z/-Z（ADC生値）
  [29]     Reserve， [30:32] フッタ 0xf00f

クォータニオンのスケールは Q14（1.0 = 16384）である．生値の4成分の二乗和は
16384²に一致する（本スクリプトで検証）．運用時の解析ツールは32768で割っていたため
ノルムが0.5となり，そこから計算したオイラー角には系統誤差がある．
ここでは16384で割って正しいクォータニオンに復元する．

出力
  data/derived/04_botan/gyro_attitude.csv  復元したクォータニオン，角速度，SAP電流
  report/gyro/gyro_quaternion.png          セッション別のクォータニオンと回転速度
  report/gyro/gyro_norm.png                スケール検証（生値ノルム）
"""
from __future__ import annotations

import struct
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data" / "downlink" / "04_botan" / "gyro" / "gyro.csv"
OUT = ROOT / "data" / "derived" / "04_botan"
FIG = ROOT / "report" / "gyro"

Q14 = 16384.0


def parse_records() -> pd.DataFrame:
    d = pd.read_csv(SRC, dtype=str).fillna("")
    rows = []
    for _, r in d[d.kind == "frame"].iterrows():
        b = bytes(int(x, 16) for x in r.hex.split())
        if len(b) < 22 or b[16:19] != bytes.fromhex("4211ff"):
            continue
        body = b[22:]
        for k in range(0, len(body) - 31, 32):
            rec = body[k:k + 32]
            if rec[0] != 0xD0 or rec[30:32] != bytes.fromhex("f00f"):
                continue
            u16 = [int.from_bytes(rec[i:i + 2], "big") for i in range(1, 29, 2)]
            i16 = [struct.unpack(">h", rec[i:i + 2])[0] for i in range(1, 29, 2)]
            rows.append(dict(
                session=r.source_log.replace("gyro_", "").replace(".txt", ""),
                downlink_date=r.downlink_date,
                pakt_num=u16[0], time_s=u16[1],
                ax=i16[2], ay=i16[3], az=i16[4],
                qi_raw=i16[5], qj_raw=i16[6], qk_raw=i16[7], qw_raw=i16[8],
                sap_xm=u16[9], sap_yp=u16[10], sap_ym=u16[11],
                sap_zp=u16[12], sap_zm=u16[13],
            ))
    df = pd.DataFrame(rows).drop_duplicates(["session", "pakt_num", "time_s"])
    return df.sort_values(["session", "pakt_num"]).reset_index(drop=True)


def add_quaternion(df: pd.DataFrame) -> pd.DataFrame:
    raw = df[["qi_raw", "qj_raw", "qk_raw", "qw_raw"]].to_numpy(float)
    df["qnorm_raw"] = np.linalg.norm(raw, axis=1)
    # 有効なクォータニオンレコード（生値ノルムがQ14の1.0＝16384に一致）．
    # 一部セッションはクォータニオン欄が空で加速度のみをダウンリンクする．
    df["q_valid"] = (df.qnorm_raw - Q14).abs() < 500
    q = np.where(df.q_valid.to_numpy()[:, None], raw / Q14, np.nan)
    for c, v in zip(["qi", "qj", "qk", "qw"], q.T):
        df[c] = np.round(v, 6)
    return df


def quat_mul(a, b):
    x1, y1, z1, w1 = a.T
    x2, y2, z2, w2 = b.T
    return np.stack([
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2], -1)


def add_rates(df: pd.DataFrame) -> pd.DataFrame:
    """隣接サンプルの相対回転 dq = q1⁻¹⊗q2 から機体系角速度を求める．"""
    for c in ["wx", "wy", "wz", "wmag"]:
        df[c] = np.nan
    for _, idx in df.groupby("session").groups.items():
        g = df.loc[idx]
        q = g[["qi", "qj", "qk", "qw"]].to_numpy(float)
        n = np.linalg.norm(q, axis=1)
        q = q / n[:, None]
        t = g.time_s.to_numpy(float)
        dt = np.diff(t)
        conj = q[:-1] * np.array([-1, -1, -1, 1])
        dq = quat_mul(conj, q[1:])
        # 符号連続化（q と -q は同じ姿勢）
        dq[dq[:, 3] < 0] *= -1
        ang = 2 * np.arctan2(np.linalg.norm(dq[:, :3], axis=1), dq[:, 3])
        with np.errstate(invalid="ignore", divide="ignore"):
            axis = dq[:, :3] / np.linalg.norm(dq[:, :3], axis=1)[:, None]
        ok = (dt >= 1) & (dt <= 4)
        w = np.degrees(ang) / dt
        rows = idx[1:]
        df.loc[rows[ok], "wmag"] = np.round(w[ok], 3)
        for j, c in enumerate(["wx", "wy", "wz"]):
            df.loc[rows[ok], c] = np.round((axis[:, j] * w)[ok], 3)
    return df


def main():
    df = parse_records()
    df = add_quaternion(df)
    df = add_rates(df)
    OUT.mkdir(parents=True, exist_ok=True)
    cols = ["session", "downlink_date", "pakt_num", "time_s",
            "qi", "qj", "qk", "qw", "qnorm_raw",
            "wx", "wy", "wz", "wmag", "ax", "ay", "az",
            "sap_xm", "sap_yp", "sap_ym", "sap_zp", "sap_zm"]
    df[cols].to_csv(OUT / "gyro_attitude.csv", index=False)

    v = df[df.q_valid]
    print(f"records: {len(df)}  sessions: {df.session.nunique()}  valid-quaternion: {len(v)}")
    print("raw quaternion norm (valid): mean %.1f  std %.1f  (Q14 => 16384)" %
          (v.qnorm_raw.mean(), v.qnorm_raw.std()))
    print("  => 正しいスケールは 1/16384．運用ツールは 1/32768 で割り，ノルム0.5の"
          "半分スケールのクォータニオンを生成していた（Euler角・回転に系統誤差）．")
    for s, g in df.groupby("session"):
        w = g.wmag.dropna()
        nv = int(g.q_valid.sum())
        tag = "quaternion" if nv else "加速度のみ(クォータニオン欄=0)"
        msg = f"  {s}: n={len(g)} valid={nv} [{tag}] t=[{g.time_s.min()},{g.time_s.max()}]s"
        if len(w):
            msg += (f" |w| median={w.median():.2f} deg/s"
                    f" [q25={w.quantile(.25):.2f}, q75={w.quantile(.75):.2f}]"
                    f" axis(median |wx|,|wy|,|wz|)="
                    f"({g.wx.abs().median():.2f},{g.wy.abs().median():.2f},{g.wz.abs().median():.2f})")
        print(msg)

    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.grid": True,
                        "grid.alpha": 0.3, "figure.dpi": 130})

    # スケール検証: 正しい÷16384 と 運用ツールの÷32768 でノルムを比較
    fig, ax = plt.subplots(figsize=(6.5, 2.9))
    vv = v.reset_index(drop=True)
    ax.plot(vv.qnorm_raw / Q14, ".", ms=3, color="#2ca02c",
            label="÷16384 (Q14, correct)")
    ax.plot(vv.qnorm_raw / 32768, ".", ms=3, color="#d62728",
            label="÷32768 (operational tool)")
    ax.axhline(1.0, color="k", lw=0.6)
    ax.axhline(0.5, color="k", lw=0.6, ls=":")
    ax.set_ylim(0, 1.2)
    ax.set_xlabel("Valid quaternion sample")
    ax.set_ylabel("Quaternion norm |q|")
    ax.legend(fontsize=8, loc="center right")
    fig.tight_layout()
    fig.savefig(FIG / "gyro_norm.png")
    plt.close(fig)

    # セッション別のクォータニオンと回転速度（有効なセッションのみ）
    sessions = [s for s in df.session.unique() if df[df.session == s].q_valid.any()]
    fig, axes = plt.subplots(2, len(sessions), figsize=(3.2 * len(sessions), 5.2),
                             sharey="row")
    if len(sessions) == 1:
        axes = axes[:, None]
    for j, s in enumerate(sessions):
        g = df[df.session == s]
        for c, lbl, col in [("qi", "i", "#1f77b4"), ("qj", "j", "#d62728"),
                            ("qk", "k", "#2ca02c"), ("qw", "w", "#9467bd")]:
            axes[0, j].plot(g.time_s, g[c], lw=0.8, label=lbl, color=col)
        axes[0, j].set_title(s, fontsize=8)
        axes[0, j].set_ylim(-1.1, 1.1)
        axes[1, j].plot(g.time_s, g.wmag, ".", ms=3, color="#9467bd")
        axes[1, j].set_xlabel("Mission time [s]")
    axes[0, 0].set_ylabel("Quaternion")
    axes[0, 0].legend(fontsize=7, ncol=4)
    axes[1, 0].set_ylabel("|ω| [deg/s]")
    fig.tight_layout()
    fig.savefig(FIG / "gyro_quaternion.png")
    plt.close(fig)
    print("wrote", OUT / "gyro_attitude.csv", "and figures in", FIG)


if __name__ == "__main__":
    main()
