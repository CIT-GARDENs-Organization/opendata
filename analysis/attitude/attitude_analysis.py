"""沿磁力線制御（Z軸を磁力線に沿わせる受動磁気姿勢制御）の成立性を，
公開用の日次HK（外面パネル温度の日別中央値）と軌道要素から検証する．

入力
  data/derived/<sat>/hk_daily.csv  日次HK（day = 放出後日数）
  data/derived/<sat>/orbit.csv     TLE由来の遠地点・近地点高度
  data/recorded/<sat>/tle.csv      TLE取得ログ（β角計算用．重複を除いてエポック・傾斜角・昇交点赤経を取り出す）

出力
  data/derived/<sat>/beta_daily.csv               日別のβ角，日照率，各面の軌道平均入射量（モデル値）
  data/derived/00_compare/attitude_summary.csv    号機別の統計量
  data/derived/00_compare/attitude_alpha_eps.csv  面別の有効α/ε
  report/attitude/attitude_*.png                  図

モデル
  1. 軌道: 円軌道近似．半長径はorbit.csvの遠地点・近地点平均，
     昇交点赤経はTLEの値を補間し，TLEの外側はJ2摂動の解析式で積分する．
  2. 太陽方向: 天文年鑑の低精度式（誤差 < 0.01°）．
  3. 地磁気: IGRF-13（ppigrf）．利用できない場合は傾斜双極子で代用する．
  4. 姿勢仮説: 衛星Z軸が磁力線 B̂ に一致し，Z軸まわりに回転する．
     ・±Z面の太陽入射 ∝ max(0, ±ŝ·B̂)
     ・側面（±X, ±Y）はスピン平均で sqrt(1-(ŝ·B̂)²)/π
  5. 対立仮説: (a) 慣性空間に固定したスピン軸，(b) 軌道面法線に沿うスピン軸，
     (c) 姿勢と無関係（ランダムタンブリング）．
"""
from __future__ import annotations

import math
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RECORDED = DATA / "recorded"  # 地上局DB等に記録されたTLE（tle.csv）
OUT = DATA / "derived"  # 2次データ（入力の日次HK・軌道と，解析出力の表）
COMPARE = OUT / "00_compare"  # 号機間比較の表
FIG = ROOT / "report" / "attitude"  # 図（報告書と同じ場所）
REPORT = FIG

RE = 6378.137  # km
MU = 398600.4418  # km^3/s^2
J2 = 1.08263e-3
SOLAR = 1361.0  # W/m^2
ALBEDO = 0.30
EARTH_IR = 237.0  # W/m^2
SIGMA = 5.670374e-8

SATS = {
    "01_kashiwa": dict(name="KASHIWA", deploy="2024-04-11T10:41:00Z", end="2024-08-05T05:59:00Z"),
    "02_sakura": dict(name="SAKURA", deploy="2024-08-29T09:45:00Z", end="2024-11-21T12:30:00Z"),
    "03_yomogi": dict(name="YOMOGI", deploy="2024-12-09T11:15:00Z", end="2025-03-25T03:57:00Z"),
    "04_botan": dict(name="BOTAN", deploy="2025-10-10T09:50:58Z", end="2026-02-27T00:00:00Z"),
}
FACES = ["x_plus", "x_minus", "y_plus", "y_minus", "z_plus", "z_minus"]
FACE_LABEL = {"x_plus": "+X", "x_minus": "-X", "y_plus": "+Y",
              "y_minus": "-Y", "z_plus": "+Z", "z_minus": "-Z"}
# 外面構成: +X面はアンテナ付きFR4基板（太陽電池なし）．KASHIWA/SAKURA/YOMOGIは残り5面がトリプルジャンクション太陽電池．
# BOTANは+Z面のみCIGS太陽電池（HKの cigs_temp_median は+Z面の温度と一致する）．
FACE_MATERIAL = {"x_plus": "FR4 (antenna)", "x_minus": "TJ cell", "y_plus": "TJ cell", "y_minus": "TJ cell",
                 "z_plus": "TJ cell (BOTAN: CIGS)", "z_minus": "TJ cell"}
# 温度場の解釈前に校正確認が必要とされた面（report/operations/report.md）．統計から除外する．
EXCLUDE = {"03_yomogi": {"x_plus"}}

try:
    import ppigrf  # type: ignore

    HAVE_IGRF = True
except Exception:  # pragma: no cover
    HAVE_IGRF = False


# --------------------------------------------------------------------------- time
def parse_utc(s: str) -> datetime:
    return datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)


def jd_from_dt(dt: datetime) -> float:
    return dt.timestamp() / 86400.0 + 2440587.5


def gmst_rad(jd: np.ndarray) -> np.ndarray:
    t = (jd - 2451545.0) / 36525.0
    g = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * t**2
    return np.deg2rad(g % 360.0)


# --------------------------------------------------------------------------- sun
def sun_vector(jd: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """太陽方向単位ベクトル（地心慣性系）と赤経・赤緯（rad）を返す．"""
    n = jd - 2451545.0
    L = np.deg2rad((280.460 + 0.9856474 * n) % 360)
    g = np.deg2rad((357.528 + 0.9856003 * n) % 360)
    lam = L + np.deg2rad(1.915) * np.sin(g) + np.deg2rad(0.020) * np.sin(2 * g)
    eps = np.deg2rad(23.439 - 0.0000004 * n)
    s = np.stack([np.cos(lam), np.cos(eps) * np.sin(lam),
                 np.sin(eps) * np.sin(lam)], axis=-1)
    ra = np.arctan2(s[..., 1], s[..., 0])
    dec = np.arcsin(s[..., 2])
    return s, ra, dec


NORAD = {"01_kashiwa": 59508, "02_sakura": 60954, "03_yomogi": 62298, "04_botan": 65942}


def load_tle_elements(sat: str) -> pd.DataFrame:
    """data/recorded/<sat>/tle.csv（取得ログ）から，該当カタログ番号の重複を除いた
    TLEを取り出し，エポック・軌道傾斜角・昇交点赤経に変換する．"""
    d = pd.read_csv(RECORDED / sat / "tle.csv")
    d = d[pd.to_numeric(d.line1.str[2:7], errors="coerce") == NORAD[sat]]
    d = d.drop_duplicates(["line1", "line2"])
    yy = 2000 + d.line1.str[18:20].astype(int)
    doy = d.line1.str[20:32].astype(float)
    epoch = pd.to_datetime(yy.astype(str), format="%Y") + pd.to_timedelta(doy - 1, unit="D")
    return pd.DataFrame({
        "epoch_utc": epoch.dt.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "inc_deg": d.line2.str[8:16].astype(float),
        "raan_deg": d.line2.str[17:25].astype(float),
    }).sort_values("epoch_utc").reset_index(drop=True)


# --------------------------------------------------------------------------- orbit
class Orbit:
    """円軌道近似．a(t) は orbit.csv，Ω(t) は TLE 補間＋J2 積分．"""

    def __init__(self, sat: str, deploy: datetime):
        self.deploy = deploy
        orb = pd.read_csv(OUT / sat / "orbit.csv")
        t = np.array([jd_from_dt(parse_utc(x)) for x in orb.time_utc])
        a = RE + (orb.apogee_km.values + orb.perigee_km.values) / 2.0
        order = np.argsort(t)
        self.t_a, self.a_tab = t[order], a[order]
        tle = load_tle_elements(sat)
        tt = np.array([jd_from_dt(parse_utc(x)) for x in tle.epoch_utc])
        order = np.argsort(tt)
        self.t_tle = tt[order]
        self.inc_tle = np.deg2rad(tle.inc_deg.values[order])
        self.n_tle = len(tt)
        # 昇交点赤経の巻き戻し: TLE間隔が長い場合に備え，J2積分の予測値に最も近い枝を選ぶ
        raw = np.deg2rad(tle.raan_deg.values[order])
        w = np.empty_like(raw)
        w[0] = raw[0]
        for k in range(1, len(raw)):
            pred = self._integrate(self.t_tle[k - 1], w[k - 1], self.t_tle[k])
            w[k] = raw[k] + 2 * math.pi * \
                round((pred - raw[k]) / (2 * math.pi))
        self.raan_tle = w

    def _integrate(self, t0: float, w0: float, t1: float) -> float:
        """Ω̇ を Simpson 法で t0 から t1 まで積分する（刻み 0.05 日）．"""
        steps = max(1, int(abs(t1 - t0) / 0.05))
        h = (t1 - t0) / steps
        w, t = w0, t0
        for _ in range(steps):
            w += h * (self.raan_dot(t) + 4 *
                      self.raan_dot(t + h / 2) + self.raan_dot(t + h)) / 6
            t += h
        return w

    def a_of(self, jd):
        return np.interp(jd, self.t_a, self.a_tab)

    def inc_of(self, jd):
        return np.interp(jd, self.t_tle, self.inc_tle)

    def raan_dot(self, jd):  # rad/day
        a = self.a_of(jd)
        n = math.sqrt(MU / a**3) * 86400.0
        return -1.5 * n * J2 * (RE / a) ** 2 * math.cos(float(self.inc_of(jd)))

    def raan_of(self, jd: float) -> float:
        """TLEの範囲内は補間，外側は最寄りTLEからJ2で積分する．"""
        if self.t_tle[0] <= jd <= self.t_tle[-1]:
            return float(np.interp(jd, self.t_tle, self.raan_tle))
        if jd < self.t_tle[0]:
            return float(self._integrate(self.t_tle[0], self.raan_tle[0], jd))
        return float(self._integrate(self.t_tle[-1], self.raan_tle[-1], jd))

    def closure_check(self) -> float | None:
        """最初のTLEから最後のTLEまで積分したときの昇交点赤経の誤差（deg）．"""
        if self.n_tle < 2:
            return None
        w = self._integrate(self.t_tle[0], self.raan_tle[0], self.t_tle[-1])
        return math.degrees((w - self.raan_tle[-1] + math.pi) % (2 * math.pi) - math.pi)


def beta_angle(raan, inc, ra_sun, dec_sun):
    return np.arcsin(np.cos(dec_sun) * np.sin(inc) * np.sin(raan - ra_sun) + np.sin(dec_sun) * np.cos(inc))


# --------------------------------------------------------------------------- geomagnetic field
G10, G11, H11 = -29350.0, -1410.3, 4545.5  # IGRF-13, 2025.0 (nT)
DIPOLE_M = np.array([G11, H11, G10]) / math.sqrt(G10**2 + G11**2 + H11**2)


def bfield_eci(r_eci: np.ndarray, jd: np.ndarray) -> np.ndarray:
    """衛星位置での地磁気ベクトル（地心慣性系，nT）．"""
    theta = gmst_rad(jd)
    c, s = np.cos(theta), np.sin(theta)
    x, y, z = r_eci[:, 0], r_eci[:, 1], r_eci[:, 2]
    # ECI -> ECEF
    xe = c * x + s * y
    ye = -s * x + c * y
    ze = z
    r = np.sqrt(xe**2 + ye**2 + ze**2)
    lat = np.arcsin(ze / r)
    lon = np.arctan2(ye, xe)
    if HAVE_IGRF:
        date = datetime(2000, 1, 1) + \
            timedelta(days=float(jd.mean()) - 2451544.5)
        be, bn, bu = ppigrf.igrf(np.degrees(
            lon), np.degrees(lat), r - RE, date)
        be, bn, bu = (np.asarray(v).reshape(-1) for v in (be, bn, bu))
        # ENU -> ECEF
        east = np.stack([-np.sin(lon), np.cos(lon), np.zeros_like(lon)], -1)
        north = np.stack([-np.sin(lat) * np.cos(lon), -
                         np.sin(lat) * np.sin(lon), np.cos(lat)], -1)
        up = np.stack([xe, ye, ze], -1) / r[:, None]
        b_ecef = be[:, None] * east + bn[:, None] * north + bu[:, None] * up
    else:
        rhat = np.stack([xe, ye, ze], -1) / r[:, None]
        m = DIPOLE_M
        mr = rhat @ m
        b_ecef = (3 * mr[:, None] * rhat - m[None, :]) * \
            (RE / r[:, None]) ** 3 * 30000.0
    # ECEF -> ECI
    bx = c * b_ecef[:, 0] - s * b_ecef[:, 1]
    by = s * b_ecef[:, 0] + c * b_ecef[:, 1]
    return np.stack([bx, by, b_ecef[:, 2]], -1)


# --------------------------------------------------------------------------- view factor
def plate_earth_view_factor_table(alt_km: float, n_theta: int = 91):
    """平板の法線と天底のなす角 θ に対する，平板から地球への形態係数．"""
    r = RE + alt_km
    rho = math.asin(RE / r)
    psi = np.linspace(0, rho, 80)
    chi = np.linspace(0, 2 * math.pi, 160, endpoint=False)
    PSI, CHI = np.meshgrid(psi, chi, indexing="ij")
    dpsi = psi[1] - psi[0]
    dchi = chi[1] - chi[0]
    thetas = np.linspace(0, math.pi, n_theta)
    F = np.empty(n_theta)
    for k, th in enumerate(thetas):
        cosang = np.cos(PSI) * math.cos(th) + np.sin(PSI) * \
            math.sin(th) * np.cos(CHI)
        F[k] = np.sum(np.clip(cosang, 0, None) * np.sin(PSI)) * \
            dpsi * dchi / math.pi
    return thetas, F


# --------------------------------------------------------------------------- daily model
def daily_geometry(orbit: Orbit, jd_center: float, n_per_orbit: int = 120):
    """1日分の軌道をサンプルし，β角，日照率，各面の入射量を返す．"""
    a = float(orbit.a_of(jd_center))
    inc = float(orbit.inc_of(jd_center))
    n_rad_s = math.sqrt(MU / a**3)
    period_day = 2 * math.pi / n_rad_s / 86400.0
    n_orbits = int(round(1.0 / period_day))
    npts = n_orbits * n_per_orbit
    jd = jd_center - 0.5 + np.arange(npts) / npts
    raan0 = orbit.raan_of(jd_center)
    raan = raan0 + orbit.raan_dot(jd_center) * (jd - jd_center)
    u = n_rad_s * (jd - jd_center) * 86400.0
    P = np.stack([np.cos(raan), np.sin(raan), np.zeros_like(raan)], -1)
    Q = np.stack([-math.cos(inc) * np.sin(raan), math.cos(inc) *
                 np.cos(raan), math.sin(inc) * np.ones_like(raan)], -1)
    rhat = np.cos(u)[:, None] * P + np.sin(u)[:, None] * Q
    r_eci = rhat * a
    s, ra, dec = sun_vector(jd)
    beta = beta_angle(raan, inc, ra, dec)
    # 日照判定（円筒影）
    rs = np.sum(rhat * s, -1)
    perp = np.linalg.norm(rhat * a - (rs * a)[:, None] * s, axis=-1)
    sunlit = (rs > 0) | (perp > RE)
    # 磁力線
    B = bfield_eci(r_eci, jd)
    bhat = B / np.linalg.norm(B, axis=-1, keepdims=True)
    c = np.sum(s * bhat, -1)  # ŝ·B̂
    # 形態係数
    thetas, Ftab = plate_earth_view_factor_table(a - RE)
    nadir = -rhat
    cos_bn = np.sum(bhat * nadir, -1)
    F_bp = np.interp(np.arccos(np.clip(cos_bn, -1, 1)), thetas, Ftab)
    F_bm = np.interp(np.arccos(np.clip(-cos_bn, -1, 1)), thetas, Ftab)
    # 側面: B̂ に垂直な面の方位平均
    e1 = np.cross(bhat, nadir)
    e1n = np.linalg.norm(e1, axis=-1, keepdims=True)
    e1 = np.where(e1n > 1e-6, e1 / np.maximum(e1n, 1e-12),
                  np.array([1.0, 0, 0]))
    e2 = np.cross(bhat, e1)
    F_side = np.zeros(npts)
    alb_side = np.zeros(npts)
    sun_side = np.zeros(npts)
    nphi = 24
    for phi in np.linspace(0, 2 * math.pi, nphi, endpoint=False):
        nrm = math.cos(phi) * e1 + math.sin(phi) * e2
        cn = np.sum(nrm * nadir, -1)
        Fs = np.interp(np.arccos(np.clip(cn, -1, 1)), thetas, Ftab)
        F_side += Fs / nphi
        alb_side += Fs * np.clip(rs, 0, None) / nphi
        sun_side += np.clip(np.sum(nrm * s, -1), 0, None) / nphi
    cosz = np.clip(rs, 0, None)  # 衛星直下点の太陽天頂角余弦
    out = dict(
        beta_deg=float(np.degrees(beta.mean())),
        sunlit_fraction=float(sunlit.mean()),
        alt_km=a - RE,
        raan_deg=math.degrees(raan0) % 360,
        sB_mean_sunlit=float(np.mean(c * sunlit)),
        # 太陽入射（太陽定数を1とした軌道平均）
        q_sun_bp=float(np.mean(np.clip(c, 0, None) * sunlit)),
        q_sun_bm=float(np.mean(np.clip(-c, 0, None) * sunlit)),
        q_sun_side=float(np.mean(sun_side * sunlit)),
        # アルベド（アルベド係数×太陽定数を1とした軌道平均）
        q_alb_bp=float(np.mean(F_bp * cosz)),
        q_alb_bm=float(np.mean(F_bm * cosz)),
        q_alb_side=float(np.mean(alb_side)),
        # 地球赤外（地球放射を1とした軌道平均 = 形態係数平均）
        F_bp=float(np.mean(F_bp)),
        F_bm=float(np.mean(F_bm)),
        F_side=float(np.mean(F_side)),
        # 対立仮説用: 日照中の太陽方向平均，軌道面法線
        sun_x=float(np.mean(s[:, 0] * sunlit)),
        sun_y=float(np.mean(s[:, 1] * sunlit)),
        sun_z=float(np.mean(s[:, 2] * sunlit)),
        sin_beta_sunlit=float(np.mean(np.sin(beta) * sunlit)),
    )
    return out


def build_daily(sat: str) -> pd.DataFrame:
    meta = SATS[sat]
    deploy = parse_utc(meta["deploy"])
    orbit = Orbit(sat, deploy)
    hk = pd.read_csv(OUT / sat / "hk_daily.csv")
    rows = []
    for day in hk.day.values:
        jd_c = jd_from_dt(deploy) + float(day) + 0.5
        g = daily_geometry(orbit, jd_c)
        g["day"] = int(day)
        g["date_utc"] = (
            deploy + timedelta(days=float(day) + 0.5)).strftime("%Y-%m-%d")
        rows.append(g)
    df = pd.DataFrame(rows)
    df.attrs["closure_deg"] = orbit.closure_check()
    df.attrs["n_tle"] = orbit.n_tle
    return df


# --------------------------------------------------------------------------- statistics
def linfit(x, y):
    """y = k x + b の最小二乗，R²，Pearson r，p値，k の標準誤差．"""
    from scipy import stats

    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 4:
        return dict(k=np.nan, b=np.nan, r2=np.nan, r=np.nan, p=np.nan, se=np.nan, n=len(x))
    res = stats.linregress(x, y)
    return dict(k=res.slope, b=res.intercept, r2=res.rvalue**2, r=res.rvalue, p=res.pvalue, se=res.stderr, n=len(x))


def multifit(X, y):
    """多変数最小二乗の自由度調整済みR²（対立仮説の比較用）．"""
    m = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
    X, y = X[m], y[m]
    A = np.column_stack([X, np.ones(len(y))])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    n, p = len(y), A.shape[1]
    r2 = 1 - ss_res / ss_tot
    r2adj = 1 - (1 - r2) * (n - 1) / (n - p)
    return dict(r2=r2, r2adj=r2adj, coef=coef, n=n)


def alpha_over_eps(T_c, q_sun, q_alb, F):
    """一面平板の放射平衡から有効 α/ε を求める．"""
    T = T_c + 273.15
    num = SIGMA * T**4 - EARTH_IR * F
    den = SOLAR * (q_sun + ALBEDO * q_alb)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 5.0, num / den, np.nan)


def analyze(sat: str, df: pd.DataFrame, hk: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    m = hk.merge(df, on="day")
    excl = EXCLUDE.get(sat, set())
    faces = [f for f in FACES if f not in excl]
    T = {f: m[f"{f}_temp_median"].astype(float).values for f in FACES}
    # 主検定: ΔT_Z vs 磁力線モデルの入射差
    dTz = T["z_plus"] - T["z_minus"]
    A_B = SOLAR * ((m.q_sun_bp - m.q_sun_bm) + ALBEDO *
                   (m.q_alb_bp - m.q_alb_bm)).values
    fitB = linfit(A_B, dTz)
    # 対立仮説
    fit_inertial = multifit(m[["sun_x", "sun_y", "sun_z"]].values, dTz)
    fit_normal = linfit(m.sin_beta_sunlit.values, dTz)
    fit_beta = linfit(m.beta_deg.values, dTz)
    # 側面: スピン平均なら ±X, ±Y の差は小さく，モデルとも無相関
    dTx = T["x_plus"] - T["x_minus"]
    dTy = T["y_plus"] - T["y_minus"]
    side_faces = [f for f in ["x_plus", "x_minus",
                              "y_plus", "y_minus"] if f not in excl]
    side_mat = np.column_stack([T[f] for f in side_faces])
    side_spread = np.nanmedian(side_mat.max(1) - side_mat.min(1))
    # 各面とβ角・日照率の相関
    corr_beta = {f: linfit(m.beta_deg.values, T[f])["r"] for f in faces}
    corr_abs_beta = {f: linfit(np.abs(m.beta_deg.values), T[f])[
        "r"] for f in faces}
    corr_sun = {f: linfit(m.sunlit_fraction.values, T[f])["r"] for f in faces}
    mean_T = np.nanmean(np.column_stack([T[f] for f in faces]), axis=1)
    fit_mean_sun = linfit(m.sunlit_fraction.values, mean_T)
    fit_mean_beta = linfit(np.abs(m.beta_deg.values), mean_T)
    # 有効 α/ε（極性は主検定の符号で決める）
    pol = 1 if (np.isfinite(fitB["k"]) and fitB["k"] >= 0) else -1
    zp, zm = ("bp", "bm") if pol > 0 else ("bm", "bp")
    ae = {}
    for f in faces:
        key = zp if f == "z_plus" else zm if f == "z_minus" else "side"
        ae[f] = alpha_over_eps(
            T[f], m[f"q_sun_{key}"].values, m[f"q_alb_{key}"].values, m[f"F_{key}"].values)
    # 等温立方体としての全体 α/ε
    n_face = len(faces)
    q_sun_tot = m.q_sun_bp + m.q_sun_bm + 4 * m.q_sun_side
    q_alb_tot = m.q_alb_bp + m.q_alb_bm + 4 * m.q_alb_side
    F_tot = m.F_bp + m.F_bm + 4 * m.F_side
    scale = 6.0 / n_face  # 除外面がある場合も6面で放射するとみなす
    ae_body = (6 * SIGMA * (mean_T + 273.15) ** 4 - EARTH_IR *
               F_tot.values) / (SOLAR * (q_sun_tot + ALBEDO * q_alb_tot).values)
    ae_df = pd.DataFrame(
        {"day": m.day, **{f"ae_{f}": ae[f] for f in faces}, "ae_body": ae_body})
    summary = dict(
        sat=SATS[sat]["name"],
        n_days=len(m),
        day_first=int(m.day.min()),
        day_last=int(m.day.max()),
        n_tle=df.attrs.get("n_tle"),
        raan_closure_deg=df.attrs.get("closure_deg"),
        beta_min=m.beta_deg.min(),
        beta_max=m.beta_deg.max(),
        sunlit_min=m.sunlit_fraction.min(),
        sunlit_max=m.sunlit_fraction.max(),
        dTz_rms=float(np.sqrt(np.nanmean(dTz**2))),
        dTz_mean=float(np.nanmean(dTz)),
        dTz_std=float(np.nanstd(dTz)),
        dTx_std=float(np.nanstd(dTx)),
        dTy_std=float(np.nanstd(dTy)),
        dTx_mean=float(np.nanmean(dTx)),
        dTy_mean=float(np.nanmean(dTy)),
        side_spread_median=float(side_spread),
        B_k=fitB["k"],
        B_k_se=fitB["se"],
        B_b=fitB["b"],
        B_r=fitB["r"],
        B_r2=fitB["r2"],
        B_p=fitB["p"],
        polarity="+Z→+B" if pol > 0 else "+Z→-B",
        inertial_r2adj=fit_inertial["r2adj"],
        normal_r2=fit_normal["r2"],
        beta_r=fit_beta["r"],
        beta_r2=fit_beta["r2"],
        meanT_vs_sunlit_r=fit_mean_sun["r"],
        meanT_vs_absbeta_r=fit_mean_beta["r"],
        **{f"meanT_{f}": float(np.nanmean(T[f])) for f in FACES},
        **{f"rbeta_{f}": corr_beta.get(f, np.nan) for f in FACES},
        **{f"rabsbeta_{f}": corr_abs_beta.get(f, np.nan) for f in FACES},
        **{f"rsun_{f}": corr_sun.get(f, np.nan) for f in FACES},
        **{f"ae_{f}": float(np.nanmedian(ae[f])) if f in ae else np.nan for f in FACES},
        **{f"ae_iqr_{f}": float(np.nanpercentile(ae[f], 75) - np.nanpercentile(ae[f], 25)) if f in ae else np.nan for f in FACES},
        ae_body=float(np.nanmedian(ae_body)),
        ae_body_iqr=float(np.nanpercentile(ae_body, 75) -
                          np.nanpercentile(ae_body, 25)),
    )
    # BOTAN: CIGS温度センサと+Z面温度の整合
    if "cigs_temp_median" in m:
        cg = m.cigs_temp_median.astype(float).values
        summary["cigs_vs_zplus_r"] = float(np.corrcoef(cg, T["z_plus"])[0, 1])
        summary["cigs_minus_zplus_mean"] = float(np.nanmean(cg - T["z_plus"]))
    m["dTz"] = dTz
    m["A_B"] = A_B
    m["dTx"] = dTx
    m["dTy"] = dTy
    m["meanT"] = mean_T
    m = m.merge(ae_df, on="day")
    return summary, m


# --------------------------------------------------------------------------- plots
def make_plots(results: dict[str, pd.DataFrame], summaries: pd.DataFrame):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.grid": True,
                        "grid.alpha": 0.3, "figure.dpi": 130})
    colors = {"KASHIWA": "#1f77b4", "SAKURA": "#d62728",
              "YOMOGI": "#2ca02c", "BOTAN": "#9467bd"}
    names = {k: SATS[k]["name"] for k in results}

    # 1. β角と日照率
    fig, ax = plt.subplots(2, 1, figsize=(8, 5.5), sharex=True)
    for k, m in results.items():
        ax[0].plot(m.day, m.beta_deg, color=colors[names[k]], label=names[k])
        ax[1].plot(m.day, m.sunlit_fraction * 100, color=colors[names[k]])
    ax[0].axhline(0, color="k", lw=0.5)
    ax[0].set_ylabel("Beta angle [deg]")
    ax[1].set_ylabel("Sunlit fraction [%]")
    ax[1].set_xlabel("Days after deployment")
    ax[0].legend(ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "attitude_beta.png")
    plt.close(fig)

    # 2. ΔT_Z 時系列 vs モデル（号機別）
    fig, axes = plt.subplots(len(results), 1, figsize=(
        8, 2.4 * len(results)), sharex=False)
    for ax, (k, m) in zip(np.atleast_1d(axes), results.items()):
        s = summaries.loc[summaries.sat == names[k]].iloc[0]
        ax.plot(m.day, m.dTz, "o-", ms=2.5, lw=0.8,
                color=colors[names[k]], label="HK: T(+Z) - T(-Z)")
        ax.plot(m.day, s.B_k * m.A_B + s.B_b, "k--", lw=1.0,
                label=f"B-aligned model fit (r={s.B_r:+.2f})")
        ax2 = ax.twinx()
        ax2.plot(m.day, m.beta_deg, color="gray", lw=0.8, alpha=0.6)
        ax2.set_ylabel("Beta [deg]", color="gray")
        ax2.grid(False)
        ax.set_ylabel("dT_Z [degC]")
        ax.set_title(names[k], loc="left", fontsize=9)
        ax.legend(fontsize=7, loc="upper left")
    np.atleast_1d(axes)[-1].set_xlabel("Days after deployment")
    fig.tight_layout()
    fig.savefig(FIG / "attitude_zdiff.png")
    plt.close(fig)

    # 3. 散布図: ΔT_Z vs モデル入射差
    fig, axes = plt.subplots(
        1, len(results), figsize=(3.0 * len(results), 3.0))
    for ax, (k, m) in zip(np.atleast_1d(axes), results.items()):
        s = summaries.loc[summaries.sat == names[k]].iloc[0]
        sc = ax.scatter(m.A_B, m.dTz, c=m.beta_deg,
                        cmap="coolwarm", s=12, vmin=-75, vmax=75)
        xx = np.linspace(np.nanmin(m.A_B), np.nanmax(m.A_B), 10)
        ax.plot(xx, s.B_k * xx + s.B_b, "k--", lw=0.8)
        ax.set_title(f"{names[k]}  r={s.B_r:+.2f}, p={s.B_p:.1e}", fontsize=8)
        ax.set_xlabel("Model q(+B) - q(-B) [W/m2]")
    np.atleast_1d(axes)[0].set_ylabel("T(+Z) - T(-Z) [degC]")
    cb = fig.colorbar(sc, ax=axes, shrink=0.9, pad=0.02)
    cb.set_label("Beta [deg]")
    fig.savefig(FIG / "attitude_zdiff_scatter.png", bbox_inches="tight")
    plt.close(fig)

    # 4. 面別温度 vs β角
    fig, axes = plt.subplots(len(results), 6, figsize=(
        13, 2.2 * len(results)), sharey="row")
    for i, (k, m) in enumerate(results.items()):
        for j, f in enumerate(FACES):
            ax = axes[i, j]
            y = m[f"{f}_temp_median"].astype(float)
            ax.scatter(m.beta_deg, y, s=6, color=colors[names[k]])
            if f in EXCLUDE.get(k, set()):
                ax.set_facecolor("#f4f4f4")
            if i == 0:
                ax.set_title(FACE_LABEL[f])
            if j == 0:
                ax.set_ylabel(f"{names[k]}\nT [degC]")
            if i == len(results) - 1:
                ax.set_xlabel("Beta [deg]")
    fig.tight_layout()
    fig.savefig(FIG / "attitude_faces_beta.png")
    plt.close(fig)

    # 5. 面別の有効 α/ε
    fig, ax = plt.subplots(figsize=(8, 3.6))
    w = 0.2
    for i, (k, m) in enumerate(results.items()):
        s = summaries.loc[summaries.sat == names[k]].iloc[0]
        vals = [s[f"ae_{f}"] for f in FACES] + [s["ae_body"]]
        errs = [s[f"ae_iqr_{f}"] / 2 for f in FACES] + [s["ae_body_iqr"] / 2]
        x = np.arange(7) + (i - 1.5) * w
        ax.bar(x, vals, w, yerr=errs,
               color=colors[names[k]], label=names[k], capsize=2)
    ax.set_xticks(np.arange(7))
    ax.set_xticklabels([FACE_LABEL[f] for f in FACES] + ["body"])
    ax.set_ylabel("Effective alpha/eps [-]")
    ax.legend(ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "attitude_alpha_eps.png")
    plt.close(fig)

    # 6. 側面の対称性: ±X, ±Y 差 vs ±Z 差
    fig, axes = plt.subplots(1, len(results), figsize=(
        3.0 * len(results), 2.8), sharey=True)
    for ax, (k, m) in zip(np.atleast_1d(axes), results.items()):
        ax.plot(m.day, m.dTz, color="k", lw=1.0, label="+Z - -Z")
        ax.plot(m.day, m.dTx, color="#ff7f0e", lw=0.8, label="+X - -X")
        ax.plot(m.day, m.dTy, color="#17becf", lw=0.8, label="+Y - -Y")
        ax.axhline(0, color="gray", lw=0.5)
        ax.set_title(names[k], fontsize=9)
        ax.set_xlabel("Days after deployment")
    np.atleast_1d(axes)[0].set_ylabel("Pair difference [degC]")
    np.atleast_1d(axes)[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "attitude_pair_diff.png")
    plt.close(fig)


# --------------------------------------------------------------------------- main
def main():
    print("IGRF:", "ppigrf" if HAVE_IGRF else "tilted dipole fallback")
    results, summaries = {}, []
    for sat in SATS:
        df = build_daily(sat)
        cols = ["day", "date_utc", "beta_deg", "sunlit_fraction", "alt_km", "raan_deg", "sB_mean_sunlit",
                "q_sun_bp", "q_sun_bm", "q_sun_side", "q_alb_bp", "q_alb_bm", "q_alb_side", "F_bp", "F_bm", "F_side"]
        (OUT / sat).mkdir(parents=True, exist_ok=True)
        df[cols].round(5).to_csv(OUT / sat / "beta_daily.csv", index=False)
        hk = pd.read_csv(OUT / sat / "hk_daily.csv")
        s, m = analyze(sat, df, hk)
        summaries.append(s)
        results[sat] = m
        print(f"{s['sat']}: n={s['n_days']} TLE={s['n_tle']} closure={s['raan_closure_deg']} "
              f"beta[{s['beta_min']:.1f},{s['beta_max']:.1f}] dTz_rms={s['dTz_rms']:.2f} "
              f"k={s['B_k']:.4f}±{s['B_k_se']:.4f} r={s['B_r']:+.3f} p={s['B_p']:.2e} pol={s['polarity']} "
              f"inertial_r2adj={s['inertial_r2adj']:.3f} normal_r2={s['normal_r2']:.3f} beta_r={s['beta_r']:+.3f} "
              f"side_spread={s['side_spread_median']:.2f} dTx_std={s['dTx_std']:.2f} dTy_std={s['dTy_std']:.2f} "
              f"ae_body={s['ae_body']:.2f}")
    summ = pd.DataFrame(summaries)
    COMPARE.mkdir(parents=True, exist_ok=True)
    summ.round(4).to_csv(COMPARE / "attitude_summary.csv", index=False)
    ae_cols = ["sat"] + [f"ae_{f}" for f in FACES] + \
        [f"ae_iqr_{f}" for f in FACES] + ["ae_body", "ae_body_iqr"]
    summ[ae_cols].round(3).to_csv(COMPARE / "attitude_alpha_eps.csv", index=False)
    make_plots(results, summ)
    return summ, results


if __name__ == "__main__":
    main()
