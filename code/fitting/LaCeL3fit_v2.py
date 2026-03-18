import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, find_peaks

# ============================================================
# Paths (EDIT)
# ============================================================
file_exp = "/mnt/data/SYC/xasfitting/exp/LaCe_Ce_L3_merged.txt"
root_dir = "/mnt/data/SYC/xasfitting/full_LaCe_doublelayer/full_mu"

OUT_DIR = "/mnt/data/SYC/xasfitting/newversion/zif8result"
os.makedirs(OUT_DIR, exist_ok=True)
OUT_CSV = os.path.join(OUT_DIR, "CeL3_peak_ranking.csv")

# ============================================================
# Windows
# ============================================================
PLOT_WIN  = (5700, 5810)
WIN1      = (5718, 5728)
WIN2      = (5728, 5742)
WIN3      = (5745, 5770)
ALIGN_WIN = (5715, 5735)

# ============================================================
# Parameters
# ============================================================
SIM_INITIAL_SHIFT = -3.5  # optional coarse shift

SG_WINDOW = 15
SG_POLY = 3
PEAK_DISTANCE = 3

# peak sensitivity
PROM_FRAC_E1 = 0.02
PROM_FRAC_E2 = 0.004  
PROM_FRAC_E3 = 0.006

# rescue gates (anti-false-peak)
CURV_MIN_NORM_E2 = 0.10
CURV_MIN_NORM_E3 = 0.10
BUMP_MIN_FRAC_E2 = 0.020
BUMP_MIN_FRAC_E3 = 0.020

SLOPE_BUMP_FRAC_E2 = 0.015
SLOPE_BUMP_FRAC_E3 = 0.015

# residual rescue (VERY LAST resort; still gated)
RESID_SIGMA_FACTOR_E2 = 2.0
RESID_SIGMA_FACTOR_E3 = 2.0
RESID_CENTER_FRACTION = 0.20  # forbid picking near window edges

# alignment scan
SHIFT_SCAN_MIN = -6.0
SHIFT_SCAN_MAX =  6.0
SHIFT_SCAN_STEP = 0.05
SHIFT_ABS_MAX = 4.0
CORR_MIN = 0.55

# scoring (soft penalty)
PENALTY_MISSING_0 = 0.0
PENALTY_MISSING_1 = 50.0
PENALTY_MISSING_2 = 200.0
PENALTY_MISSING_3 = 1000.0

TOP_N_PLOTS = 200
DPI = 180

# ============================================================
# Utilities
# ============================================================
def smooth(y):
    if SG_WINDOW % 2 == 0 or SG_WINDOW < 5 or SG_WINDOW >= len(y):
        return y.copy()
    return savgol_filter(y, SG_WINDOW, SG_POLY)

def parabola_refine(x, y, i):
    if i <= 0 or i >= len(x) - 1:
        return float(x[i])
    y0, y1, y2 = y[i-1], y[i], y[i+1]
    denom = y0 - 2*y1 + y2
    if abs(denom) < 1e-12:
        return float(x[i])
    dx = x[i+1] - x[i]
    return float(x[i] + 0.5*(y0 - y2)/denom * dx)

def load_xy(path):
    a = np.loadtxt(path)
    return a[:,0], a[:,1]

def clip_xy(x, y, win):
    m = (x >= win[0]) & (x <= win[1])
    return x[m], y[m]

def corrcoef_safe(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 8 or len(b) < 8:
        return -1.0
    a = a - np.mean(a)
    b = b - np.mean(b)
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na < 1e-12 or nb < 1e-12:
        return -1.0
    return float(np.dot(a, b) / (na * nb))

def best_shift_by_corr(x_exp, y_exp_s, x_sim, y_sim_s, win):
    m = (x_exp >= win[0]) & (x_exp <= win[1])
    if not np.any(m):
        return np.nan, -1.0
    xw = x_exp[m]
    yw = y_exp_s[m]

    best_s, best_c = np.nan, -1.0
    for s in np.arange(SHIFT_SCAN_MIN, SHIFT_SCAN_MAX + 0.5*SHIFT_SCAN_STEP, SHIFT_SCAN_STEP):
        ysi = np.interp(xw, x_sim + s, y_sim_s, left=np.nan, right=np.nan)
        ok = np.isfinite(ysi)
        if ok.sum() < 12:
            continue
        c = corrcoef_safe(yw[ok], ysi[ok])
        if c > best_c:
            best_c, best_s = c, float(s)
    return best_s, best_c

# ============================================================
# Peak picking: strict one per window, with 3-step rescue
# ============================================================
def pick_peak_strict_rescue(x, y_s, win, prom_frac,
                            curv_min_norm, bump_min_frac,
                            slope_bump_frac,
                            resid_sigma_factor):
    """
    Return one peak energy in win, or NaN.
    Order:
      1) find_peaks
      2) curvature (2nd derivative min) with gates
      3) slope zero-crossing with gates
      4) residual max after linear baseline with gates (last resort)
    """
    m = (x >= win[0]) & (x <= win[1])
    if not np.any(m):
        return np.nan

    xs, ys = x[m], y_s[m]
    amp = float(np.max(ys) - np.min(ys)) + 1e-12
    prom_min = float(prom_frac) * amp
    base = float(np.median(ys))

    # ---- 1) classical peaks
    pk, props = find_peaks(ys, distance=PEAK_DISTANCE, prominence=prom_min)
    if len(pk) > 0:
        prominences = props.get("prominences", np.ones(len(pk)))
        j = int(pk[np.argmax(prominences)])
        idx = np.where(m)[0][0] + j
        E = parabola_refine(x, y_s, idx)
        if win[0] <= E <= win[1]:
            return float(E)

    # ---- 2) curvature rescue
    d2 = np.gradient(np.gradient(y_s, x), x)
    d2w = d2[m]
    j = int(np.argmin(d2w))
    idx = np.where(m)[0][0] + j
    E = parabola_refine(x, y_s, idx)

    curv_norm = float(abs(d2w[j])) / (float(np.max(np.abs(d2w))) + 1e-12)
    bump = float(y_s[idx] - base) / amp
    if (curv_norm >= curv_min_norm) and (bump >= bump_min_frac) and (win[0] <= E <= win[1]):
        return float(E)

    # ---- 3) slope zero-crossing rescue
    d1 = np.gradient(y_s, x)
    d1w = d1[m]
    sgn = np.sign(d1w)
    for k in range(1, len(sgn)):
        if sgn[k] == 0:
            sgn[k] = sgn[k-1]
    cross = np.where((sgn[:-1] > 0) & (sgn[1:] < 0))[0]
    if len(cross) > 0:
        # choose crossing with strongest bump
        best_idx = None
        best_b = -1e9
        for c in cross:
            ii = int(c + 1)
            b = float(ys[ii] - base) / amp
            if b > best_b:
                best_b = b
                best_idx = ii
        if best_idx is not None and best_b >= slope_bump_frac:
            idx = np.where(m)[0][0] + best_idx
            E = parabola_refine(x, y_s, idx)
            if win[0] <= E <= win[1]:
                return float(E)

    # ---- 4) residual max after linear baseline (last resort)
    # baseline fit on window
    p = np.polyfit(xs, ys, 1)
    resid = ys - (p[0]*xs + p[1])
    sigma = float(np.std(resid)) + 1e-12
    j = int(np.argmax(resid))

    # forbid edges
    wlen = win[1] - win[0]
    left_forbid  = win[0] + RESID_CENTER_FRACTION*wlen
    right_forbid = win[1] - RESID_CENTER_FRACTION*wlen

    if (resid[j] >= resid_sigma_factor * sigma) and (xs[j] >= left_forbid) and (xs[j] <= right_forbid):
        idx = np.where(m)[0][0] + j
        E = parabola_refine(x, y_s, idx)
        if win[0] <= E <= win[1]:
            return float(E)

    return np.nan

def pick_peak_simple(x, y_s, win, prom_frac):
    # used for E1 mostly; no crazy rescue needed
    m = (x >= win[0]) & (x <= win[1])
    if not np.any(m):
        return np.nan
    xs, ys = x[m], y_s[m]
    amp = float(np.max(ys) - np.min(ys)) + 1e-12
    prom_min = float(prom_frac) * amp
    pk, props = find_peaks(ys, distance=PEAK_DISTANCE, prominence=prom_min)
    if len(pk) == 0:
        # fallback: max in window (white-line is usually clear)
        j = int(np.argmax(ys))
    else:
        prominences = props.get("prominences", np.ones(len(pk)))
        j = int(pk[np.argmax(prominences)])
    idx = np.where(m)[0][0] + j
    return parabola_refine(x, y_s, idx)

def score(missing, d2, d3, corr_ok):
    if not corr_ok:
        return 1e9
    base = 0.0
    if np.isfinite(d2): base += abs(d2)
    if np.isfinite(d3): base += abs(d3)
    penalty = {0: PENALTY_MISSING_0, 1: PENALTY_MISSING_1, 2: PENALTY_MISSING_2, 3: PENALTY_MISSING_3}[missing]
    return base + penalty

# ============================================================
# Load experimental
# ============================================================
x_exp, y_exp = load_xy(file_exp)
x_exp, y_exp = clip_xy(x_exp, y_exp, PLOT_WIN)
y_exp_s = smooth(y_exp)

E1_exp = pick_peak_simple(x_exp, y_exp_s, WIN1, 0.02)
E2_exp = pick_peak_strict_rescue(
    x_exp, y_exp_s, WIN2, 0.02,
    curv_min_norm=0.10, bump_min_frac=0.02,
    slope_bump_frac=0.015,
    resid_sigma_factor=2.0
)
E3_exp = pick_peak_strict_rescue(
    x_exp, y_exp_s, WIN3, 0.02,
    curv_min_norm=0.10, bump_min_frac=0.02,
    slope_bump_frac=0.015,
    resid_sigma_factor=2.0
)

print("Locked experimental peaks:")
print(f"E1_exp={E1_exp:.4f}, E2_exp={E2_exp:.4f}, E3_exp={E3_exp:.4f}")

# ============================================================
# Walk simulations
# ============================================================
rows = []
plot_cache = []

for root, _, _ in os.walk(root_dir):
    if not root.endswith("L3"):
        continue
    f = os.path.join(root, "xmu_extracted.dat")
    if not os.path.exists(f):
        continue

    tag = os.path.basename(root)

    try:
        x_sim, y_sim = load_xy(f)
        x_sim = x_sim + SIM_INITIAL_SHIFT
        x_sim, y_sim = clip_xy(x_sim, y_sim, PLOT_WIN)
        if len(x_sim) < 30:
            continue

        y_sim_s = smooth(y_sim)

        shift, corr = best_shift_by_corr(x_exp, y_exp_s, x_sim, y_sim_s, ALIGN_WIN)
        corr_ok = np.isfinite(shift) and (corr >= CORR_MIN) and (abs(shift) <= SHIFT_ABS_MAX)

        if corr_ok:
            x_al = x_sim + shift

            # E1 (for plotting + missing accounting)
            E1_sim = pick_peak_simple(x_al, y_sim_s, WIN1, PROM_FRAC_E1)

            # E2/E3 (hard rescue)
            E2_sim = pick_peak_strict_rescue(
                x_al, y_sim_s, WIN2, PROM_FRAC_E2,
                curv_min_norm=CURV_MIN_NORM_E2, bump_min_frac=BUMP_MIN_FRAC_E2,
                slope_bump_frac=SLOPE_BUMP_FRAC_E2,
                resid_sigma_factor=RESID_SIGMA_FACTOR_E2
            )
            E3_sim = pick_peak_strict_rescue(
                x_al, y_sim_s, WIN3, PROM_FRAC_E3,
                curv_min_norm=CURV_MIN_NORM_E3, bump_min_frac=BUMP_MIN_FRAC_E3,
                slope_bump_frac=SLOPE_BUMP_FRAC_E3,
                resid_sigma_factor=RESID_SIGMA_FACTOR_E3
            )
        else:
            x_al = x_sim
            E1_sim = E2_sim = E3_sim = np.nan

        missing = int(not np.isfinite(E1_sim)) + int(not np.isfinite(E2_sim)) + int(not np.isfinite(E3_sim))
        dE2 = (E2_sim - E2_exp) if (np.isfinite(E2_sim) and np.isfinite(E2_exp)) else np.nan
        dE3 = (E3_sim - E3_exp) if (np.isfinite(E3_sim) and np.isfinite(E3_exp)) else np.nan

        sc = score(missing, dE2, dE3, corr_ok)

        rows.append({
            "tag": tag,
            "file": f,
            "sim_initial_shift_eV": float(SIM_INITIAL_SHIFT),
            "align_shift_eV": float(shift) if np.isfinite(shift) else np.nan,
            "align_corr": float(corr),
            "corr_ok": int(corr_ok),
            "E1_exp_eV": float(E1_exp),
            "E2_exp_eV": float(E2_exp) if np.isfinite(E2_exp) else np.nan,
            "E3_exp_eV": float(E3_exp) if np.isfinite(E3_exp) else np.nan,
            "E1_sim_eV": float(E1_sim) if np.isfinite(E1_sim) else np.nan,
            "E2_sim_eV": float(E2_sim) if np.isfinite(E2_sim) else np.nan,
            "E3_sim_eV": float(E3_sim) if np.isfinite(E3_sim) else np.nan,
            "dE2_eV": float(dE2) if np.isfinite(dE2) else np.nan,
            "dE3_eV": float(dE3) if np.isfinite(dE3) else np.nan,
            "missing": int(missing),
            "score": float(sc),
        })

        plot_cache.append({
            "tag": tag,
            "x_al": x_al,
            "y_sim": y_sim,
            "E1_sim": E1_sim, "E2_sim": E2_sim, "E3_sim": E3_sim,
            "shift": shift, "corr": corr, "corr_ok": corr_ok,
            "missing": missing, "score": sc
        })

    except Exception as e:
        print(f"[WARN] Failed on {f}: {e}")

# ============================================================
# Save ranking
# ============================================================
df = pd.DataFrame(rows)
df = df.sort_values(["missing", "score", "tag"], ascending=[True, True, True]).reset_index(drop=True)
df["rank"] = np.arange(1, len(df) + 1)
df.to_csv(OUT_CSV, index=False)
print("Saved:", OUT_CSV)

# ============================================================
# Plot top-N
# ============================================================
TOP_N = min(TOP_N_PLOTS, len(df))
top_tags = df.head(TOP_N)["tag"].tolist()
rank_map = {t: i+1 for i, t in enumerate(top_tags)}

print(f"Saving top {TOP_N} plots...")

for it in plot_cache:
    tag = it["tag"]
    if tag not in rank_map:
        continue
    r = rank_map[tag]

    x_al = it["x_al"]
    y_sim = it["y_sim"]
    E1_sim, E2_sim, E3_sim = it["E1_sim"], it["E2_sim"], it["E3_sim"]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x_exp, y_exp, label="experiment", lw=1.5)
    ax.plot(x_al, y_sim, label="simulation (aligned)", lw=1.5)

    exp_es = [E1_exp, E2_exp, E3_exp]
    exp_es = [e for e in exp_es if np.isfinite(e)]
    ax.scatter(exp_es, np.interp(exp_es, x_exp, y_exp), s=80, marker="o", label="exp peaks")

    sim_es = [E1_sim, E2_sim, E3_sim]
    sim_es = [e for e in sim_es if np.isfinite(e)]
    if sim_es:
        ax.scatter(sim_es, np.interp(sim_es, x_al, y_sim), s=80, marker="x", label="sim peaks")

    ax.set_xlim(PLOT_WIN)
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("mu(E)")
    ax.set_title(
        f"Rank {r:02d} | {tag}\n"
        f"score={it['score']:.3f} | missing={it['missing']} | corr={it['corr']:.3f} | "
        f"shift={SIM_INITIAL_SHIFT:+.2f}{(it['shift'] if np.isfinite(it['shift']) else 0.0):+.2f} eV"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, f"Ce-L3-{tag}.png"), dpi=DPI)
    plt.close(fig)

print(f"Top-{TOP_N} plots saved to {OUT_DIR}")
