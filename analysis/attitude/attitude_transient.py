"""軌道周期を分解した生HKによる姿勢解析（沿磁力線制御の判別，姿勢安定までの時間，
スピン速度，α と ε の分離推定）．

入力（ダウンリンクした生HK．data/downlink/<sat>/hk/）
  01_kashiwa/hk/hk.csv    KASHIWA 41byte HK を復号したもの（90 s）
  01_kashiwa/hk/hshk.csv  KASHIWA 高速HK（6 s）
  02_sakura/hk/hk.csv     SAKURA 全期間HK（HK_All）
  02_sakura/hk/hshk_*.csv SAKURA 高速HK（6 s）
  03_yomogi/hk/hk.csv     YOMOGI 全期間HK
  04_botan/hk/hk.csv      BOTAN 全期間HK
本スクリプトの出力（2次データ）
  data/derived/<sat>/attitude_daily.csv   日別の姿勢指標（磁力線モデル相関，面切替回数，食時間など）
  data/derived/<sat>/eclipse_events.csv   HKで検出した食の中心時刻と継続時間，モデル値
  data/derived/00_compare/attitude_transient_summary.csv, attitude_spin.csv, attitude_alpha_eps_transient.csv, attitude_decay.csv
  report/attitude/attitude_settling.png, attitude_orbit_example.png, attitude_spin.png, attitude_thermal_fit.png, attitude_pointing.png, attitude_decay.png

手法の要点
  1. 食の検出: 5面の太陽電池電圧がすべて 2 V 未満を食とする．
  2. 軌道位相: 食の中心で衛星は反太陽方向にあるとして，円軌道モデルの位相を決める．
     TLE の昇交点赤経・半長径と組み合わせ，各HK時刻の位置，太陽方向，IGRF-13 の磁力線方向を求める．
  3. 磁力線判別: 1日ごとに ΔT_Z = T(+Z) − T(−Z) と ŝ·B̂ の相関を求める．Z軸が磁力線に沿えば
     1軌道に2回の反転が現れ相関が高い．軌道面法線に沿う場合はこの変調がない．
  4. 面切替: 日照中に +Z/−Z のどちらが明るいかの切替回数を軌道あたりで数える．
  5. スピン速度: 高速HK（6 s）の側面電圧を Lomb-Scargle 周期解析する．
  6. α, ε: 各面を1節点とみなし，食中の冷却率から ε，日照中の加熱率から α を同定する．
     熱容量 C_p を仮定するため，ε と α は C_p に比例する．α/ε は C_p に依らない．
"""
from __future__ import annotations

import math
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
ROOT = Path(__file__).resolve().parents[2]
from attitude_analysis import (  # noqa: E402
    ALBEDO, COMPARE, EARTH_IR, EXCLUDE, FACES, FACE_LABEL, FIG, MU, OUT, RE, SATS, SIGMA, SOLAR,
    Orbit, bfield_eci, jd_from_dt, parse_utc, plate_earth_view_factor_table, sun_vector,
)

DL = ROOT / "data" / "downlink"
A_FACE = 0.01  # m^2 (1U)
CP_PANEL = 45.0  # J/K 仮定（FR4 1.6 mm + セル ≈ 45 g × 1000 J/kg/K）
V_DARK, V_LIT = 2.0, 3.5
MASS = {"KASHIWA": None, "SAKURA": 1.089, "YOMOGI": 1.1427, "BOTAN": None}


# --------------------------------------------------------------------------- load
def load_hk(sat: str) -> pd.DataFrame:
    """統一列名: t (UTC), T_xp..T_zm, T_bpb, T_bat, V_xm..V_zm, bat_v, bat_i"""
    if sat == "01_kashiwa":
        # time_utc は day カウンタと起動エポック（A: 2024-03-30 01:37Z, B: 2024-05-07 11:52Z）から復元済み
        d = pd.read_csv(DL / sat / "hk" / "hk.csv")
        d = d[(d.day < 150)].copy()
        out = pd.DataFrame({
            "t": pd.to_datetime(d.time_utc, utc=True), "T_xp": d.x_plus, "T_xm": d.x_minus, "T_yp": d.y_plus, "T_ym": d.y_minus,
            "T_zp": d.z_plus, "T_zm": d.z_minus, "T_bpb": d.bpb, "T_bat": d.bat_temp,
            "V_xm": d.v_xm, "V_yp": d.v_yp, "V_ym": d.v_ym, "V_zp": d.v_zp, "V_zm": d.v_zm, "bat_v": d.bat_v, "bat_i": d.bat_i})
    elif sat == "02_sakura":
        d = pd.read_csv(DL / sat / "hk" / "hk.csv")
        t = pd.to_datetime(d.time_utc, utc=True, errors="coerce")
        out = pd.DataFrame({"t": t, "T_xp": d["+X_HK"], "T_xm": d["-X_HK"], "T_yp": d["+Y_HK"], "T_ym": d["-Y_HK"],
                            "T_zp": d["+Z_HK"], "T_zm": d["-Z_HK"], "T_bpb": d["BPB_HK"], "T_bat": d["BAT_HK"],
                            "V_xm": d["-X_V"], "V_yp": d["+Y_V"], "V_ym": d["-Y_V"], "V_zp": d["+Z_V"], "V_zm": d["-Z_V"],
                            "bat_v": d["BAT_V"], "bat_i": d["BAT_I"]})
    elif sat == "03_yomogi":
        d = pd.read_csv(DL / sat / "hk" / "hk.csv", low_memory=False)
        t = pd.to_datetime(d.time_utc, utc=True, errors="coerce")
        out = pd.DataFrame({"t": t, "T_xp": d["'+X_Temp"], "T_xm": d["'-X_Temp"], "T_yp": d["'+Y_Temp"],
                            "T_ym": d["'-Y_Temp"], "T_zp": d["'+Z_Temp"], "T_zm": d["'-Z_Temp"], "T_bpb": d["BPB_Temp"],
                            "T_bat": d["BAT_Temp"], "V_xm": d["'-X_V"], "V_yp": d["'+Y_V"], "V_ym": d["'-Y_V"],
                            "V_zp": d["'+Z_V"], "V_zm": d["'-Z_V"], "bat_v": d["BAT_V"], "bat_i": d["BAT_I"]})
    else:
        d = pd.read_csv(DL / sat / "hk" / "hk.csv", low_memory=False)
        t = pd.to_datetime(d.time_utc, utc=True, errors="coerce")
        out = pd.DataFrame({"t": t, "T_xp": d["'+X_Temp [℃]"], "T_xm": d["'-X_Temp [℃]"],
                            "T_yp": d["'+Y_Temp [℃]"], "T_ym": d["'-Y_Temp [℃]"], "T_zp": d["'+Z_Temp [℃]"],
                            "T_zm": d["'-Z_Temp [℃]"], "T_bpb": d["BPB_Temp [℃]"], "T_bat": d["BAT_Temp [℃]"],
                            "V_xm": d["'-X_V [V]"], "V_yp": d["'+Y_V [V]"], "V_ym": d["'-Y_V [V]"], "V_zp": d["'+Z_V [V]"],
                            "V_zm": d["'-Z_V [V]"], "bat_v": d["BAT_V [V]"], "bat_i": d["BAT_I [A]"]})
    for c in out.columns:
        if c != "t":
            out[c] = pd.to_numeric(out[c], errors="coerce")
    out = out.dropna(subset=["t"]).sort_values(
        "t").drop_duplicates("t").reset_index(drop=True)
    deploy = pd.Timestamp(SATS[sat]["deploy"].replace("Z", "+00:00"))
    out["days"] = (out.t - deploy).dt.total_seconds() / 86400.0
    out["jd"] = out.t.map(lambda x: jd_from_dt(x.to_pydatetime()))
    vmax = out[["V_xm", "V_yp", "V_ym", "V_zp", "V_zm"]].max(axis=1)
    out["dark"] = vmax < V_DARK
    out["lit"] = vmax > V_LIT
    return out


# --------------------------------------------------------------------------- eclipses
def find_eclipses(hk: pd.DataFrame) -> pd.DataFrame:
    """連続する暗サンプルを食とみなす．前後に日照サンプルがあるものだけ採用する．"""
    t = hk.t.values.astype("datetime64[s]").astype(np.int64)
    dark = hk.dark.values
    ev = []
    i, n = 0, len(hk)
    while i < n:
        if not dark[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and dark[j + 1] and (t[j + 1] - t[j]) <= 240:
            j += 1
        dur = t[j] - t[i]
        ok = (12 * 60 <= dur <= 45 * 60) and i > 0 and j < n - 1 and (t[i] - t[i - 1]) <= 300 and (t[j + 1] - t[j]) <= 300 \
            and hk.lit.values[i - 1] and hk.lit.values[j + 1]
        if ok:
            ev.append(dict(i0=i, i1=j, t_mid=pd.Timestamp(
                (t[i] + t[j]) / 2, unit="s", tz="UTC"), dur_min=(dur + 90) / 60.0))
        i = j + 1
    return pd.DataFrame(ev)


# --------------------------------------------------------------------------- orbit model
class PhasedOrbit:
    """食の中心を位相の基準に用いた円軌道モデル．"""

    def __init__(self, sat: str, ecl: pd.DataFrame):
        self.orbit = Orbit(sat, parse_utc(SATS[sat]["deploy"]))
        self.anchors = []  # (jd_mid, u_mid)
        for _, e in ecl.iterrows():
            jd = jd_from_dt(e.t_mid.to_pydatetime())
            raan = self.orbit.raan_of(jd)
            inc = float(self.orbit.inc_of(jd))
            P, Q, nrm = self.basis(raan, inc)
            s, _, _ = sun_vector(np.array([jd]))
            s = s[0]
            r = -(s - np.dot(s, nrm) * nrm)
            r /= np.linalg.norm(r)
            self.anchors.append((jd, math.atan2(np.dot(r, Q), np.dot(r, P))))
        self.anchor_jd = np.array([a[0] for a in self.anchors])
        self.anchor_u = np.array([a[1] for a in self.anchors])

    @staticmethod
    def basis(raan, inc):
        P = np.array([math.cos(raan), math.sin(raan), 0.0])
        Q = np.array([-math.cos(inc) * math.sin(raan),
                     math.cos(inc) * math.cos(raan), math.sin(inc)])
        return P, Q, np.cross(P, Q)

    def n_rad_s(self, jd):
        a = float(self.orbit.a_of(jd))
        return math.sqrt(MU / a**3)

    def phase_consistency(self):
        """隣接する食の間で位相を伝播したときの残差（deg）．"""
        res = []
        for k in range(len(self.anchor_jd) - 1):
            dt = (self.anchor_jd[k + 1] - self.anchor_jd[k]) * 86400
            if dt > 6 * 3600:
                continue
            n = self.n_rad_s((self.anchor_jd[k] + self.anchor_jd[k + 1]) / 2)
            pred = self.anchor_u[k] + n * dt
            res.append(math.degrees(
                (pred - self.anchor_u[k + 1] + math.pi) % (2 * math.pi) - math.pi))
        return np.array(res)

    def geometry(self, jd: np.ndarray):
        """各時刻の ŝ·B̂，モデル日照，太陽天頂角余弦，各面の地球形態係数を返す．"""
        k = np.searchsorted(self.anchor_jd, jd)
        k = np.clip(k, 0, len(self.anchor_jd) - 1)
        kprev = np.clip(k - 1, 0, len(self.anchor_jd) - 1)
        use_prev = np.abs(
            jd - self.anchor_jd[kprev]) < np.abs(jd - self.anchor_jd[k])
        k = np.where(use_prev, kprev, k)
        dt_anchor = (jd - self.anchor_jd[k]) * 86400
        valid = np.abs(dt_anchor) < 12 * 3600
        a = self.orbit.a_of(jd)
        n = np.sqrt(MU / a**3)
        u = self.anchor_u[k] + n * dt_anchor
        raan = np.array([self.orbit.raan_of(x) for x in jd])
        inc = self.orbit.inc_of(jd)
        P = np.stack([np.cos(raan), np.sin(raan), np.zeros_like(raan)], -1)
        Q = np.stack([-np.cos(inc) * np.sin(raan), np.cos(inc)
                     * np.cos(raan), np.sin(inc)], -1)
        rhat = np.cos(u)[:, None] * P + np.sin(u)[:, None] * Q
        s, _, _ = sun_vector(jd)
        rs = np.sum(rhat * s, -1)
        perp = np.linalg.norm(
            rhat * a[:, None] - (rs * a)[:, None] * s, axis=-1)
        sunlit = (rs > 0) | (perp > RE)
        # IGRF は日単位で呼ぶ
        B = np.zeros_like(rhat)
        day_idx = np.floor(jd).astype(int)
        for dval in np.unique(day_idx):
            m = day_idx == dval
            B[m] = bfield_eci(rhat[m] * a[m, None], jd[m])
        bhat = B / np.linalg.norm(B, axis=-1, keepdims=True)
        c = np.sum(s * bhat, -1)
        nadir = -rhat
        cos_bn = np.sum(bhat * nadir, -1)
        # 形態係数（高度別テーブル）
        thetas, Ftab = plate_earth_view_factor_table(float(np.median(a)) - RE)
        F_bp = np.interp(np.arccos(np.clip(cos_bn, -1, 1)), thetas, Ftab)
        F_bm = np.interp(np.arccos(np.clip(-cos_bn, -1, 1)), thetas, Ftab)
        # 側面: B̂ に垂直な面の方位平均
        e1 = np.cross(bhat, nadir)
        e1n = np.linalg.norm(e1, axis=-1, keepdims=True)
        e1 = np.where(e1n > 1e-6, e1 / np.maximum(e1n, 1e-12),
                      np.array([1.0, 0, 0]))
        e2 = np.cross(bhat, e1)
        F_side = np.zeros(len(jd))
        for phi in np.linspace(0, 2 * math.pi, 12, endpoint=False):
            nrm = math.cos(phi) * e1 + math.sin(phi) * e2
            F_side += np.interp(np.arccos(np.clip(np.sum(nrm *
                                nadir, -1), -1, 1)), thetas, Ftab) / 12
        beta = np.arcsin(np.sum(np.cross(P, Q) * s, -1))
        return pd.DataFrame(dict(valid=valid, u_deg=np.degrees(u) % 360, sB=c, sunlit_model=sunlit, cosz=np.clip(rs, 0, None),
                                 F_bp=F_bp, F_bm=F_bm, F_side=F_side, beta_deg=np.degrees(beta)))


# --------------------------------------------------------------------------- daily metrics
def switches_per_orbit(sign: np.ndarray, t_s: np.ndarray, period_s: float) -> float:
    m = sign != 0
    if m.sum() < 10:
        return np.nan
    sg = sign[m]
    ts = t_s[m]
    changes = np.sum((sg[1:] != sg[:-1]) & ((ts[1:] - ts[:-1]) < 600))
    span_orbits = np.sum(np.minimum(ts[1:] - ts[:-1], 600)) / period_s
    return changes / span_orbits if span_orbits > 0.5 else np.nan


def lagged_corr(x, y, lags):
    best = (np.nan, 0)
    for L in lags:
        if L > 0:
            a, b = x[L:], y[:-L]
        else:
            a, b = x, y
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < 30:
            continue
        r = np.corrcoef(a[m], b[m])[0, 1]
        if not np.isfinite(best[0]) or abs(r) > abs(best[0]):
            best = (r, L)
    return best


def daily_metrics(hk: pd.DataFrame, geo: pd.DataFrame, ecl: pd.DataFrame, po: PhasedOrbit) -> pd.DataFrame:
    df = pd.concat([hk.reset_index(drop=True),
                   geo.reset_index(drop=True)], axis=1)
    df["dTz"] = df.T_zp - df.T_zm
    df["day"] = np.floor(df.days).astype(int)
    df["t_s"] = df.t.values.astype("datetime64[s]").astype(np.int64)
    rows = []
    for day, g in df.groupby("day"):
        g = g.sort_values("t_s")
        period = 2 * math.pi / po.n_rad_s(float(g.jd.mean()))
        lit = g[g.lit & g.valid]
        rB, lag = lagged_corr(lit.dTz.values, lit.sB.values, range(
            0, 7)) if len(lit) else (np.nan, 0)
        kB = np.nan
        if len(lit) > 30:
            x = lit.sB.values[:-lag] if lag > 0 else lit.sB.values
            yv = lit.dTz.values[lag:] if lag > 0 else lit.dTz.values
            mm = np.isfinite(x) & np.isfinite(yv)
            if mm.sum() > 30 and np.std(x[mm]) > 0.05:
                kB = float(np.polyfit(x[mm], yv[mm], 1)[0])
        # 連続サンプルのみを使う（ラグ相関のため）
        rB0 = np.corrcoef(lit.dTz, lit.sB)[0, 1] if len(lit) > 30 else np.nan
        # 軌道面法線モデル: ΔT_Z が sin β に比例 → 1日内ではほぼ一定 → 相関は定義できないので，
        # 代わりに ΔT_Z の日内変動のうち ŝ·B̂ で説明できる分散比を R² として記録する
        sz = np.where(np.abs(g.V_zp - g.V_zm) > 0.5,
                      np.sign(g.V_zp - g.V_zm), 0)
        sy = np.where(np.abs(g.V_yp - g.V_ym) > 0.5,
                      np.sign(g.V_yp - g.V_ym), 0)
        litmask = g.lit.values
        sw_z = switches_per_orbit(
            np.where(litmask, sz, 0), g.t_s.values, period)
        sw_y = switches_per_orbit(
            np.where(litmask, sy, 0), g.t_s.values, period)
        both_z = np.mean((g.V_zp > 4.3) & (g.V_zm > 4.3)
                         ) if litmask.sum() else np.nan
        agree = np.mean(g.sunlit_model.values[g.valid.values] == ~
                        g.dark.values[g.valid.values]) if g.valid.sum() else np.nan
        e = ecl[(ecl.t_mid >= g.t.min()) & (ecl.t_mid <= g.t.max())]
        rows.append(dict(day=day, n=len(g), n_lit=int(litmask.sum()), n_valid=int(g.valid.sum()), n_eclipse=len(e),
                         r_B=rB, k_B=kB, lag_min=lag * 1.5, r_B_lag0=rB0, dTz_std=float(np.nanstd(g.dTz)),
                         dTz_mean=float(np.nanmean(g.dTz)), sw_z=sw_z, sw_y=sw_y, both_z_lit=both_z,
                         shadow_agree=agree, beta_deg=float(g.beta_deg.mean()),
                         ecl_dur_meas=float(
                             e.dur_min.mean()) if len(e) else np.nan,
                         ecl_dur_model=float(e.dur_model.mean()) if len(e) and "dur_model" in e else np.nan))
    return pd.DataFrame(rows), df


def model_eclipse_duration(po: PhasedOrbit, jd: float) -> float:
    a = float(po.orbit.a_of(jd))
    raan = po.orbit.raan_of(jd)
    inc = float(po.orbit.inc_of(jd))
    s, _, _ = sun_vector(np.array([jd]))
    P, Q, nrm = po.basis(raan, inc)
    beta = math.asin(float(np.dot(nrm, s[0])))
    x = math.sqrt(a**2 - RE**2) / (a * math.cos(beta)
                                   ) if abs(math.cos(beta)) > 1e-6 else 2
    if x >= 1:
        return 0.0
    frac = (1 / math.pi) * math.acos(x)
    return frac * 2 * math.pi / po.n_rad_s(jd) / 60.0


# --------------------------------------------------------------------------- spin (HSHK)
def spin_periods(sat: str) -> pd.DataFrame:
    from scipy.signal import lombscargle

    rows = []
    sessions = []
    if sat == "01_kashiwa":
        d = pd.read_csv(DL / sat / "hk" / "hshk.csv")
        for s, g in d.groupby("session"):
            g = g.sort_values("elapsed_s")
            sessions.append((s, g.elapsed_s.values.astype(float), {"V_xm": g.v_xm.values, "V_yp": g.v_yp.values, "V_ym": g.v_ym.values,
                                                                   "V_zp": g.v_zp.values, "V_zm": g.v_zm.values},
                             # 2024-05-08 再起動 → 放出後日数 ≈ day + 27
                             g.day.iloc[0] + 27))
    elif sat == "02_sakura":
        for p in sorted((DL / sat / "hk").glob("hshk_*.csv")):
            g = pd.read_csv(p).sort_values("t_s")
            sessions.append((p.stem[-8:], g.t_s.values.astype(float), {"V_xm": g["-X_V"].values, "V_yp": g["+Y_V"].values,
                                                                       "V_ym": g["-Y_V"].values, "V_zp": g["+Z_V"].values, "V_zm": g["-Z_V"].values},
                             g.t_s.iloc[0] / 86400))
    for name, t, ch, mday in sessions:
        # 連続区間（ギャップ < 20 s）ごとに解析し，側面の周期性が最も強い区間を採用
        gaps = np.where(np.diff(t) > 20)[0]
        segs = [sg for sg in np.split(
            np.arange(len(t)), gaps + 1) if len(sg) >= 25]
        if not segs:
            continue
        periods = np.linspace(13, 600, 3000)
        freqs = 2 * math.pi / periods
        best = None
        for seg in segs:
            ts = t[seg] - t[seg][0]
            res = {}
            for c, v in ch.items():
                y = pd.to_numeric(pd.Series(v[seg]), errors="coerce").values
                m = np.isfinite(y)
                if m.sum() < 25 or np.mean(y[m] > V_LIT) < 0.5 or np.nanstd(y[m]) < 0.05:
                    res[c] = (np.nan, 0.0)
                    continue
                yy = y[m] - y[m].mean()
                pg = lombscargle(ts[m], yy, freqs, normalize=True)
                k = int(np.argmax(pg))
                res[c] = (float(periods[k]), float(pg[k]))
            score = np.nanmean([res[c][1] for c in ("V_xm", "V_yp", "V_ym")])
            if best is None or (np.isfinite(score) and score > best[0]):
                best = (score, seg, res)
        score, seg, res = best
        ts = t[seg] - t[seg][0]
        side = [res[c]
                for c in ("V_xm", "V_yp", "V_ym") if np.isfinite(res[c][0])]
        zp = [res[c] for c in ("V_zp", "V_zm") if np.isfinite(res[c][0])]
        rows.append(dict(sat=SATS[sat]["name"], session=name, mission_day=float(mday), n=len(seg), span_min=(ts[-1]) / 60,
                         side_period_s=float(
                             np.median([p for p, _ in side])) if side else np.nan,
                         side_power=float(
                             np.mean([w for _, w in side])) if side else np.nan,
                         z_period_s=float(
                             np.median([p for p, _ in zp])) if zp else np.nan,
                         z_power=float(
                             np.mean([w for _, w in zp])) if zp else np.nan,
                         **{f"P_{c}": res[c][0] for c in ch}, **{f"pow_{c}": res[c][1] for c in ch}))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- thermal identification
def thermal_fit(df: pd.DataFrame, sat: str, polarity: int) -> pd.DataFrame:
    """各面 1 節点モデル: C dT/dt = -ε A(σT⁴ - q_IR F) - G (T - T_bat) + α A (S q_sun + a S F cosz)"""
    g = df.sort_values("t_s").copy()
    dt = np.diff(g.t_s.values)
    ok = (dt > 60) & (dt < 130)
    out = []
    zp_key, zm_key = ("bp", "bm") if polarity > 0 else ("bm", "bp")
    for f in FACES:
        if f in EXCLUDE.get(sat, set()):
            continue
        col = {"x_plus": "T_xp", "x_minus": "T_xm", "y_plus": "T_yp",
               "y_minus": "T_ym", "z_plus": "T_zp", "z_minus": "T_zm"}[f]
        T = g[col].values + 273.15
        Tb = g.T_bat.values + 273.15
        dTdt = np.full(len(g), np.nan)
        dTdt[1:-1] = (T[2:] - T[:-2]) / (g.t_s.values[2:] - g.t_s.values[:-2])
        good = np.zeros(len(g), bool)
        good[1:-1] = ok[1:] & ok[:-1]
        if f == "z_plus":
            F = g[f"F_{zp_key}"].values
            qs = np.clip(g.sB.values * polarity, 0, None)
        elif f == "z_minus":
            F = g[f"F_{zm_key}"].values
            qs = np.clip(-g.sB.values * polarity, 0, None)
        else:
            F = g.F_side.values
            qs = np.sqrt(np.clip(1 - g.sB.values**2, 0, 1)) / math.pi
        X1 = SIGMA * T**4 - EARTH_IR * F
        X2 = T - Tb
        X3 = SOLAR * qs * g.sunlit_model.values + ALBEDO * SOLAR * F * g.cosz.values
        valid = g.valid.values & np.isfinite(dTdt) & np.isfinite(X2) & good
        dark = valid & g.dark.values & ~g.sunlit_model.values
        lit = valid & g.lit.values & g.sunlit_model.values
        if dark.sum() < 100 or lit.sum() < 100:
            continue
        # 食中: dT/dt = a1 X1 + a2 X2
        A = np.column_stack([X1[dark], X2[dark]])
        (a1, a2), *_ = np.linalg.lstsq(A, dTdt[dark], rcond=None)
        resid_d = dTdt[dark] - A @ np.array([a1, a2])
        # 日照中: dT/dt - a1 X1 - a2 X2 = a3 X3
        y = dTdt[lit] - a1 * X1[lit] - a2 * X2[lit]
        a3 = float(np.sum(X3[lit] * y) / np.sum(X3[lit] ** 2))
        resid_l = y - a3 * X3[lit]
        eps = -a1 * CP_PANEL / A_FACE
        alpha = a3 * CP_PANEL / A_FACE
        G = -a2 * CP_PANEL
        # ブートストラップで不確かさ
        rng = np.random.default_rng(0)
        eb, ab = [], []
        idx_d = np.where(dark)[0]
        idx_l = np.where(lit)[0]
        for _ in range(100):
            sd = rng.choice(idx_d, len(idx_d))
            sl = rng.choice(idx_l, len(idx_l))
            Ab = np.column_stack([X1[sd], X2[sd]])
            (b1, b2), *_ = np.linalg.lstsq(Ab, dTdt[sd], rcond=None)
            yb = dTdt[sl] - b1 * X1[sl] - b2 * X2[sl]
            b3 = np.sum(X3[sl] * yb) / np.sum(X3[sl] ** 2)
            eb.append(-b1 * CP_PANEL / A_FACE)
            ab.append(b3 * CP_PANEL / A_FACE)
        out.append(dict(sat=SATS[sat]["name"], face=FACE_LABEL[f], n_dark=int(dark.sum()), n_lit=int(lit.sum()),
                        eps=eps, eps_sd=float(np.std(eb)), alpha=alpha, alpha_sd=float(np.std(ab)),
                        alpha_over_eps=alpha / eps if eps else np.nan, G_W_per_K=G,
                        rms_resid_dark_mK_s=float(
                            np.sqrt(np.mean(resid_d**2))) * 1000,
                        rms_resid_lit_mK_s=float(
                            np.sqrt(np.mean(resid_l**2))) * 1000,
                        r2_dark=1 - np.var(resid_d) / np.var(dTdt[dark]), r2_lit=1 - np.var(resid_l) / np.var(dTdt[lit])))
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- main
def settle_day(daily: pd.DataFrame, thr=0.5, run=3):
    d = daily[daily.n_lit >= 60].sort_values("day")
    ok = (d.r_B.abs() >= thr).values
    for i in range(len(ok) - run + 1):
        if ok[i:i + run].all():
            return int(d.day.values[i])
    return None


def main():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.grid": True,
                        "grid.alpha": 0.3, "figure.dpi": 130})
    colors = {"KASHIWA": "#1f77b4", "SAKURA": "#d62728",
              "YOMOGI": "#2ca02c", "BOTAN": "#9467bd"}
    summ, dailies, examples, therm_all = [], {}, {}, []
    for sat in SATS:
        name = SATS[sat]["name"]
        print(f"=== {name}")
        hk = load_hk(sat)
        ecl = find_eclipses(hk)
        po = PhasedOrbit(sat, ecl)
        ecl["dur_model"] = [model_eclipse_duration(
            po, jd_from_dt(t.to_pydatetime())) for t in ecl.t_mid]
        cons = po.phase_consistency()
        geo = po.geometry(hk.jd.values)
        daily, df = daily_metrics(hk, geo, ecl, po)
        pol = int(np.sign(np.nanmedian(
            daily.loc[daily.n_lit >= 60, "r_B"]))) or 1
        kk = daily.k_B * pol
        k_ref = float(np.nanpercentile(kk[daily.n_lit >= 60], 90))
        daily["theta_eff_deg"] = np.degrees(
            np.arccos(np.clip(kk / k_ref, -1, 1)))
        daily.loc[daily.n_lit < 60, "theta_eff_deg"] = np.nan
        daily.attrs["k_ref"] = k_ref
        sd = settle_day(daily)
        agree = float(np.nanmean(daily.shadow_agree))
        th = thermal_fit(df, sat, pol)
        therm_all.append(th)
        spin = spin_periods(sat) if sat in (
            "01_kashiwa", "02_sakura") else pd.DataFrame()
        summ.append(dict(sat=name, n_hk=len(hk), day_first=float(hk.days.min()), day_last=float(hk.days.max()),
                         n_eclipse=len(ecl), phase_resid_rms_deg=float(np.sqrt(np.mean(cons**2))) if len(cons) else np.nan,
                         ecl_dur_meas_mean=float(ecl.dur_min.mean()), ecl_dur_model_mean=float(ecl.dur_model.mean()),
                         shadow_agree=agree, polarity="+Z→+B" if pol > 0 else "+Z→-B",
                         r_B_median=float(np.nanmedian(
                             daily.loc[daily.n_lit >= 60, "r_B"])),
                         r_B_q25=float(np.nanpercentile(
                             daily.loc[daily.n_lit >= 60, "r_B"] * pol, 25)),
                         r_B_q75=float(np.nanpercentile(
                             daily.loc[daily.n_lit >= 60, "r_B"] * pol, 75)),
                         settle_day=sd, sw_z_median=float(np.nanmedian(daily.sw_z)), sw_y_median=float(np.nanmedian(daily.sw_y)),
                         sw_z_first5=float(np.nanmedian(daily[daily.day < 5].sw_z)), sw_y_first5=float(np.nanmedian(daily[daily.day < 5].sw_y)),
                         both_z_lit_median=float(np.nanmedian(daily.both_z_lit))))
        print(f"  eclipses={len(ecl)} phase_resid_rms={summ[-1]['phase_resid_rms_deg']:.2f}deg shadow_agree={agree:.3f} "
              f"pol={summ[-1]['polarity']} r_B median={summ[-1]['r_B_median']:+.2f} settle_day={sd} "
              f"sw_z={summ[-1]['sw_z_median']:.2f} sw_y={summ[-1]['sw_y_median']:.2f}")
        print(th.round(3).to_string(index=False))
        if len(spin):
            print(spin[["session", "mission_day", "n", "side_period_s", "side_power",
                  "z_period_s", "z_power"]].round(2).to_string(index=False))
        (OUT / sat).mkdir(parents=True, exist_ok=True)
        daily.round(4).to_csv(OUT / sat / "attitude_daily.csv", index=False)
        ecl[["t_mid", "dur_min", "dur_model"]].assign(t_mid=lambda x: x.t_mid.dt.strftime("%Y-%m-%dT%H:%M:%SZ")).round(2) \
            .to_csv(OUT / sat / "eclipse_events.csv", index=False)
        if len(spin):
            spin.round(3).to_csv(
                COMPARE / f"attitude_spin_{name.lower()}.csv", index=False)
        dailies[sat] = daily
        examples[sat] = (df, pol, sd)
    summ = pd.DataFrame(summ)
    COMPARE.mkdir(parents=True, exist_ok=True)
    summ.round(4).to_csv(COMPARE / "attitude_transient_summary.csv", index=False)
    therm = pd.concat(therm_all, ignore_index=True)
    therm.round(4).to_csv(
        COMPARE / "attitude_alpha_eps_transient.csv", index=False)
    # ---- figures
    # settling
    fig, axes = plt.subplots(3, 4, figsize=(14, 7.5), sharex="col")
    for j, sat in enumerate(SATS):
        name = SATS[sat]["name"]
        d = dailies[sat]
        d = d[d.n_lit >= 60]
        _, pol, sd = examples[sat]
        ax = axes[0, j]
        ax.plot(d.day, d.r_B * pol, "o-", ms=2.5, lw=0.8, color=colors[name])
        ax.axhline(0.5, color="gray", lw=0.6, ls="--")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylim(-1, 1)
        ax.set_title(
            name + (f"   settle: day {sd}" if sd is not None else "   settle: n/a"), fontsize=9, loc="left")
        if sd is not None:
            for a in axes[:, j]:
                a.axvline(sd, color=colors[name], lw=0.8, ls=":")
        ax = axes[1, j]
        ax.plot(d.day, d.sw_z, "o-", ms=2.5, lw=0.8,
                color="k", label="+Z/-Z switches per orbit")
        ax.plot(d.day, d.sw_y, "o-", ms=2.5, lw=0.8,
                color="#17becf", label="+Y/-Y switches per orbit")
        ax.set_yscale("log")
        ax = axes[2, j]
        ax.plot(d.day, d.dTz_std, "o-", ms=2.5, lw=0.8,
                color=colors[name], label="std of dT_Z within day")
        ax.plot(d.day, d.dTz_mean, "o-", ms=2.5,
                lw=0.8, color="gray", label="mean dT_Z")
        ax.set_xlabel("Days after deployment")
    axes[0, 0].set_ylabel("corr(dT_Z, s.B) per day")
    axes[1, 0].set_ylabel("switches / orbit")
    axes[2, 0].set_ylabel("dT_Z [degC]")
    axes[1, 0].legend(fontsize=7)
    axes[2, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "attitude_settling.png")
    plt.close(fig)
    # convergence / damping: exponential fits to pointing-angle index and Z-face wobble index
    from scipy.optimize import curve_fit

    def expdec(t, yinf, A, tau):
        return yinf + A * np.exp(-t / tau)

    decay_rows = []
    fig, axes = plt.subplots(2, 4, figsize=(14, 6), sharex="col")
    for j, sat in enumerate(SATS):
        name = SATS[sat]["name"]
        d = dailies[sat]
        for i, (col, lbl, tr) in enumerate([("theta_eff_deg", "theta_eff [deg]", lambda v: v),
                                            ("sw_z", "Z-face wobble index (switches/orbit - 2)", lambda v: np.clip(v - 2, 0, None))]):
            v = d[np.isfinite(d[col]) & (d.n_lit >= 60) & (d.day <= 70)]
            t = v.day.values.astype(float)
            y = tr(v[col].values.astype(float))
            ax = axes[i, j]
            ax.plot(t, y, "o", ms=3, color=colors[name])
            fit = None
            if len(t) >= 6:
                try:
                    p0 = [np.nanmedian(y[-5:]), max(y[0] -
                                                    np.nanmedian(y[-5:]), 1e-3), 10.0]
                    popt, pcov = curve_fit(expdec, t, y, p0=p0, bounds=(
                        [-np.inf, -np.inf, 0.5], [np.inf, np.inf, 200]), maxfev=20000)
                    perr = np.sqrt(np.diag(pcov))
                    tt = np.linspace(t.min(), t.max(), 200)
                    ax.plot(tt, expdec(tt, *popt), "k--", lw=1)
                    # 過剰分が10%になる時刻
                    t_conv = float(t.min() + popt[2] * math.log(10))
                    fit = dict(sat=name, metric=col, y_inf=popt[0], A=popt[1], tau_day=popt[2], tau_sd=perr[2], t_first=float(t.min()),
                               t_conv_10pct=t_conv, n=len(t))
                    ax.set_title(
                        f"{name}: tau = {popt[2]:.1f} d, y_inf = {popt[0]:.1f}", fontsize=9, loc="left")
                except Exception as e:  # pragma: no cover
                    ax.set_title(f"{name}: fit failed", fontsize=9, loc="left")
            if fit:
                decay_rows.append(fit)
            if j == 0:
                ax.set_ylabel(lbl)
            if i == 1:
                ax.set_xlabel("Days after deployment")
    fig.tight_layout()
    fig.savefig(FIG / "attitude_decay.png")
    plt.close(fig)
    pd.DataFrame(decay_rows).round(3).to_csv(
        COMPARE / "attitude_decay.csv", index=False)
    print(pd.DataFrame(decay_rows).round(2).to_string(index=False))
    # pointing stability: effective angle between Z axis and B
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.2), sharey=True)
    for ax, sat in zip(axes, SATS):
        name = SATS[sat]["name"]
        d = dailies[sat]
        d = d[np.isfinite(d.theta_eff_deg)]
        ax.plot(d.day, d.theta_eff_deg, "o-",
                ms=2.5, lw=0.8, color=colors[name])
        ax.axhline(90, color="k", lw=0.5)
        ax.set_title(
            f"{name}  (k_ref = {dailies[sat].attrs['k_ref']:.1f} degC)", fontsize=9, loc="left")
        ax.set_xlabel("Days after deployment")
        ax.set_yticks([0, 30, 60, 90, 120, 150, 180])
    axes[0].set_ylabel("Effective Z-to-B angle [deg]")
    fig.tight_layout()
    fig.savefig(FIG / "attitude_pointing.png")
    plt.close(fig)
    # orbit examples: early day vs settled day
    fig, axes = plt.subplots(4, 2, figsize=(13, 10))
    for i, sat in enumerate(SATS):
        name = SATS[sat]["name"]
        df, pol, sd = examples[sat]
        d = dailies[sat]
        d = d[d.n_lit >= 60]
        dv = d[np.isfinite(d.r_B)]
        early = int(dv.day.min()) if len(dv) else int(d.day.min())
        late = int(dv.loc[dv.r_B.abs().idxmax(), "day"]) if len(dv) else early
        for k, day in enumerate((early, late)):
            ax = axes[i, k]
            g = df[(df.day == day) & df.valid].sort_values("t_s")
            if len(g) == 0:
                continue
            h = (g.t_s - g.t_s.min()) / 3600
            ax.plot(h, g.dTz, color=colors[name],
                    lw=1.0, label="T(+Z) - T(-Z)")
            ax2 = ax.twinx()
            ax2.plot(h, g.sB * pol, color="k", lw=0.7,
                     alpha=0.6, label="model s.B (x polarity)")
            ax2.set_ylim(-1.1, 1.1)
            ax2.grid(False)
            ax.fill_between(h, -100, 100, where=g.dark,
                            color="gray", alpha=0.15, step="mid")
            ax.set_ylim(np.nanmin(g.dTz) - 3, np.nanmax(g.dTz) + 3)
            r = d.loc[d.day == day, "r_B"].values
            ax.set_title(f"{name}  day {day}  r={r[0]:+.2f}" if len(
                r) else f"{name} day {day}", fontsize=9, loc="left")
            ax.set_xlabel("hours")
            if k == 0:
                ax.set_ylabel("dT_Z [degC]")
            else:
                ax2.set_ylabel("s.B")
    axes[0, 0].legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "attitude_orbit_example.png")
    plt.close(fig)
    # spin
    spins = [pd.read_csv(p) for p in COMPARE.glob("attitude_spin_*.csv")]
    if spins:
        sp = pd.concat(spins)
        fig, ax = plt.subplots(figsize=(7, 3.4))
        for name, g in sp.groupby("sat"):
            g = g.sort_values("mission_day")
            ax.errorbar(g.mission_day, 360 / g.side_period_s, fmt="o-",
                        ms=4, color=colors[name], label=f"{name} side faces")
            ax.plot(g.mission_day, 360 / g.z_period_s, "x", ms=5,
                    color=colors[name], alpha=0.5, label=f"{name} Z faces")
        ax.set_xlabel("Days after deployment")
        ax.set_ylabel("Apparent rate 360/period [deg/s]")
        ax.set_yscale("log")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(FIG / "attitude_spin.png")
        plt.close(fig)
        sp.round(3).to_csv(COMPARE / "attitude_spin.csv", index=False)
    # thermal
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
    w = 0.2
    for i, name in enumerate(["KASHIWA", "SAKURA", "YOMOGI", "BOTAN"]):
        t = therm[therm.sat == name].set_index("face")
        x = np.arange(6) + (i - 1.5) * w
        labels = [FACE_LABEL[f] for f in FACES]
        axes[0].bar(x, [t.eps.get(l, np.nan) for l in labels], w, yerr=[
                    t.eps_sd.get(l, 0) for l in labels], color=colors[name], label=name, capsize=2)
        axes[1].bar(x, [t.alpha.get(l, np.nan) for l in labels], w, yerr=[
                    t.alpha_sd.get(l, 0) for l in labels], color=colors[name], label=name, capsize=2)
    for ax, ttl in zip(axes, ["Emissivity eps  (C_p = 45 J/K assumed)", "Absorptivity alpha  (C_p = 45 J/K assumed)"]):
        ax.set_xticks(np.arange(6))
        ax.set_xticklabels([FACE_LABEL[f] for f in FACES])
        ax.set_title(ttl, fontsize=9)
        ax.set_ylim(0, 1.5)
    axes[0].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "attitude_thermal_fit.png")
    plt.close(fig)
    print(summ.round(3).T.to_string())


if __name__ == "__main__":
    main()
