"""BOTANジャイロのクォータニオンから，3軸矢印つきキューブの姿勢アニメーションを生成する．

- 20251012_boot: 放出直後．t≈30秒でクォータニオンの成分が飛ぶが，これは二重被覆（q と −q は
  同一姿勢）による符号反転である．キューブは滑らかに回り続け，不連続が表示上の特異点である
  ことが分かる．見やすさのため実時間の3倍に引き伸ばす．
- 20251023: 整定後．Z軸まわりの定常スピン（約18 °/s，周期20秒前後）を実時間で示す．
  スピンでスカラー部が周期的に符号を跨ぐため，生値には多数の符号反転が現れるが，
  いずれも表示上のものである．

3.5秒を超える受信欠落は動画上でも飛ばす（補間で埋めない）．

出力: report/gyro/gyro_boot_animation.mp4, gyro_20251023_animation.mp4
使い方: python3 gyro_animation.py [セッション名...]（省略時は全セッション）
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FFMpegWriter  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data" / "derived" / "04_botan" / "gyro_attitude.csv"
FIG = ROOT / "report" / "gyro"
GAP_S = 3.5  # これを超える欠落は補間せず動画上も飛ばす

SESSIONS = {
    "20251012_boot": dict(out="gyro_boot_animation.mp4", fps=25, slowmo=3.0,
                          title="BOTAN 2025-10-12 (deployment)"),
    "20251023": dict(out="gyro_20251023_animation.mp4", fps=20, slowmo=1.0,
                     title="BOTAN 2025-10-23 (settled spin)"),
}

AXCOL = {"X": "#e0524d", "Y": "#3fae5a", "Z": "#4f8cd6"}


def quat_to_matrix(q):
    x, y, z, w = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


def slerp(q0, q1, u):
    d = np.dot(q0, q1)
    if d < 0:  # 短い方の弧を通る（二重被覆の解消）
        q1 = -q1
        d = -d
    if d > 0.9995:
        q = q0 + u * (q1 - q0)
        return q / np.linalg.norm(q)
    th0 = np.arccos(np.clip(d, -1, 1))
    th = th0 * u
    q2 = q1 - q0 * d
    q2 /= np.linalg.norm(q2)
    return q0 * np.cos(th) + q2 * np.sin(th)


# --- キューブ（1U，±0.5）．+Z面（CIGS）を強調
V = np.array([[i, j, k] for i in (-.5, .5) for j in (-.5, .5) for k in (-.5, .5)])
FACES = {
    "zp": ([2, 3, 7, 6], "#9467bd"), "zm": ([0, 1, 5, 4], "#cfc3e6"),
    "xp": ([4, 5, 7, 6], "#e9e6f2"), "xm": ([0, 1, 3, 2], "#e9e6f2"),
    "yp": ([1, 3, 7, 5], "#ddd8ec"), "ym": ([0, 2, 6, 4], "#ddd8ec"),
}
AXES = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) * 1.15


def render(session: str, cfg: dict):
    d = pd.read_csv(SRC)
    g = d[d.session == session].dropna(subset=["qi", "qj", "qk", "qw"]).sort_values("time_s")
    if len(g) < 10:
        print(f"skip {session}: no valid quaternions")
        return
    t = g.time_s.to_numpy(float)
    q_raw = g[["qi", "qj", "qk", "qw"]].to_numpy(float)
    q_cont = q_raw.copy()
    for i in range(1, len(q_cont)):
        if q_cont[i] @ q_cont[i - 1] < 0:
            q_cont[i] = -q_cont[i]
    flip_t = [t[i] for i in range(1, len(q_raw)) if q_raw[i] @ q_raw[i - 1] < 0]

    # 一様時間グリッド．欠落（> GAP_S）は区間ごとに区切り，動画上も飛ばす
    fps, slowmo = cfg["fps"], cfg["slowmo"]
    seg_bounds = [0] + [i for i in range(1, len(t)) if t[i] - t[i - 1] > GAP_S] + [len(t)]
    tt = np.concatenate([
        np.arange(t[a], t[b - 1], 1.0 / fps / slowmo)
        for a, b in zip(seg_bounds[:-1], seg_bounds[1:]) if b - a > 1
    ])

    def q_at(time):
        j = np.searchsorted(t, time) - 1
        j = np.clip(j, 0, len(t) - 2)
        u = (time - t[j]) / (t[j + 1] - t[j])
        return slerp(q_cont[j], q_cont[j + 1], np.clip(u, 0, 1))

    plt.rcParams.update({"font.size": 9, "figure.dpi": 130})
    fig = plt.figure(figsize=(10, 4.6))
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    axp = fig.add_subplot(1, 2, 2)

    # 成分グラフ．生値（点線）は符号反転するが同一姿勢の別表現なので，
    # 符号連続に補正した値（実線）を主に示す
    for c, lbl, col in [(0, "i", "#1f77b4"), (1, "j", "#d62728"),
                        (2, "k", "#2ca02c"), (3, "w", "#9467bd")]:
        axp.plot(t, q_raw[:, c], lw=0.8, color=col, ls=":", alpha=0.4)
        axp.plot(t, q_cont[:, c], lw=1.3, color=col, label=lbl)
    for ft in flip_t:
        axp.axvline(ft, color="#d99a2b", lw=0.8, ls="--", alpha=0.6)
    if len(flip_t) == 1:
        axp.text(flip_t[0], 1.03, f"raw sign flip @t={flip_t[0]:.0f}s\n(corrected)",
                 color="#b8791a", fontsize=7.5, ha="center", va="bottom")
    elif flip_t:
        axp.text(0.02, 0.98, f"raw sign flips x{len(flip_t)} (dashed) — display only",
                 color="#b8791a", fontsize=7.5, ha="left", va="top", transform=axp.transAxes)
    axp.set_xlabel("Mission time [s]")
    axp.set_ylabel("Quaternion  (dotted=raw, solid=sign-corrected)")
    axp.set_ylim(-1.1, 1.25)
    axp.grid(alpha=0.3)
    axp.legend(fontsize=8, ncol=4, loc="lower right")
    cursor = axp.axvline(t[0], color="k", lw=1.4)

    writer = FFMpegWriter(fps=fps, bitrate=2400)
    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / cfg["out"]
    wmag = g.wmag.ffill().fillna(0).to_numpy(float)

    with writer.saving(fig, str(out), dpi=130):
        for time in tt:
            R = quat_to_matrix(q_at(time))
            Vr = V @ R.T
            ax3d.clear()
            for idx, col in FACES.values():
                poly = Poly3DCollection([Vr[idx]], alpha=0.9)
                poly.set_facecolor(col)
                poly.set_edgecolor("#3a3350")
                ax3d.add_collection3d(poly)
            Ar = AXES @ R.T
            for name, vec in zip("XYZ", Ar):
                ax3d.quiver(0, 0, 0, *vec, color=AXCOL[name], lw=2.5,
                            arrow_length_ratio=0.18)
                ax3d.text(*(vec * 1.18), name, color=AXCOL[name], fontsize=11,
                          fontweight="bold", ha="center", va="center")
            ax3d.set_xlim(-1.3, 1.3)
            ax3d.set_ylim(-1.3, 1.3)
            ax3d.set_zlim(-1.3, 1.3)
            ax3d.set_box_aspect((1, 1, 1))
            ax3d.set_axis_off()
            ax3d.view_init(elev=18, azim=45)
            q = q_at(time)
            wnow = np.interp(time, t, wmag)
            ax3d.set_title(f"{cfg['title']}  t = {time:5.1f} s   |ω| ≈ {wnow:4.1f} deg/s\n"
                           f"q = ({q[0]:+.2f}, {q[1]:+.2f}, {q[2]:+.2f}, {q[3]:+.2f})",
                           fontsize=9, family="monospace")
            cursor.set_xdata([time, time])
            writer.grab_frame()

    plt.close(fig)
    print(f"wrote {out} ({len(tt)} frames, {len(tt)/fps:.1f}s video, "
          f"{len(flip_t)} raw sign flips)")


if __name__ == "__main__":
    targets = sys.argv[1:] or list(SESSIONS)
    for name in targets:
        render(name, SESSIONS[name])
