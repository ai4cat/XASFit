import os
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import wofz

file1 = "/mnt/data/SYC/xasfitting/proced/ZIF-8-Ce-L3.txt"
root_dir = "/mnt/data/SYC/xasfitting/full_LaCe_doublelayer/full_mu"  
prominence = 0.03
distance = 5
fit_half_window = 10

def load_xy(path):
    d = np.loadtxt(path)
    return d[:,0], d[:,1]

def gauss(x, a, mu, sigma, c):
    return a*np.exp(-(x-mu)**2/(2*sigma**2)) + c

def voigt(x, amp, cen, sigma, gamma, offset):
    z = ((x - cen) + 1j*gamma) / (sigma*np.sqrt(2))
    return amp*np.real(wofz(z)) / (sigma*np.sqrt(2*np.pi)) + offset

def refine_peak_mu(x, y, idx, half_win):
    i0 = max(0, idx - half_win)
    i1 = min(len(x), idx + half_win + 1)
    xs, ys = x[i0:i1], y[i0:i1]
    if len(xs) < 5:
        return x[idx]
    
    #local_max_idx = np.argmax(ys)
    #mu0 = xs[local_max_idx]

    c0 = np.median(ys)
    a0 = ys.max() - c0
    mu0 = x[idx]
    dx = np.median(np.diff(xs)) if len(xs) > 1 else 1.0
    sigma0 = max(dx * half_win / 6.0, dx*1.5)
    p0 = [a0, mu0, sigma0, c0]
    bounds = ([0, xs.min(), dx*0.5, -np.inf], [np.inf, xs.max(), np.ptp(xs), np.inf])
    try:
        popt, _ = curve_fit(voigt, xs, ys, p0=p0, bounds=bounds, maxfev=2000)
        return float(popt[1])
    except Exception:
        return x[idx]

def detect_peaks_with_refine(x, y):
    kwargs = dict(prominence=prominence)
    if distance is not None: kwargs["distance"] = distance
    pk, props = find_peaks(y, **kwargs)
    mus = [refine_peak_mu(x, y, i, fit_half_window) for i in pk]
    prom = props.get("prominences", np.ones(len(pk)))
    order = np.argsort(mus)
    return np.array(mus)[order], prom[order]

def match_peaks(mus_ref, mus_cmp, max_dist=np.inf):
    if len(mus_ref)==0 or len(mus_cmp)==0:
        return [], [], np.array([]), np.arange(len(mus_ref)), np.arange(len(mus_cmp))
    D = cdist(mus_ref.reshape(-1,1), mus_cmp.reshape(-1,1), metric='euclidean')
    matched_ref, matched_cmp, dists = [], [], []
    used_ref, used_cmp = set(), set()
    pairs = [(i,j,D[i,j]) for i in range(D.shape[0]) for j in range(D.shape[1])]
    pairs.sort(key=lambda t: t[2])
    for i,j,d in pairs:
        if i in used_ref or j in used_cmp:
            continue
        if d <= max_dist:
            used_ref.add(i); used_cmp.add(j)
            matched_ref.append(i); matched_cmp.append(j); dists.append(d)
    unmatched_ref = np.array([i for i in range(len(mus_ref)) if i not in used_ref])
    unmatched_cmp = np.array([j for j in range(len(mus_cmp)) if j not in used_cmp])
    return np.array(matched_ref), np.array(matched_cmp), np.array(dists), unmatched_ref, unmatched_cmp

def detect_micro_fluctuation_single_peak(x, y, start=5730, end=5740, window_size=5, range_thresh=0.001):
    '''(max-min), no d2x because no sensitive to micro peak'''
    mask = (x >= start) & (x <= end)
    xs, ys = x[mask], y[mask]
    if len(ys) < window_size:
        return np.array([])

    max_range = 0
    peak_x = None
    for i in range(len(ys)-window_size+1):
        local_window = ys[i:i+window_size]
        local_range = np.ptp(local_window) 
        if local_range > range_thresh and local_range > max_range:
            max_range = local_range
            local_idx = i + np.argmax(local_window)
            peak_x = xs[local_idx]
    
    if peak_x is not None:
        return np.array([peak_x])
    else:
        return np.array([])

# reference
x1, y1 = load_xy(file1)
lim1 = (x1 >= 5700) & (x1 <= 5810)
x1, y1 = x1[lim1], y1[lim1]
mus1, prom1 = detect_peaks_with_refine(x1, y1)
mus1[1]=mus1[1]+0.4
mus1[2]=mus1[2]-2.5

results = []
i=1

for dirpath, dirnames, filenames in os.walk(root_dir):
    if dirpath.endswith("L3"):  # L3
        file2 = os.path.join(dirpath, "xmu_extracted.dat")
        if os.path.exists(file2):
            try:
                x2, y2 = load_xy(file2)
                x2 = x2 - 3.5
                lim2 = (x2 >= 5700) & (x2 <= 5810)
                xl, yl = x2[lim2], y2[lim2]
                lim2 = (x2 >= 5700) & (x2 <= 5770)
                x2, y2 = x2[lim2], y2[lim2]
                mus2, prom2 = detect_peaks_with_refine(x2, y2)

                micro_peaks = detect_micro_fluctuation_single_peak(x2, y2)
                mus2 = np.concatenate([mus2, micro_peaks])
                mus2 = np.sort(mus2)

                mref, mcmp, dists, unref, uncmp = match_peaks(mus1, mus2, max_dist=np.inf)

                n_ref = len(mus1)
                n_cmp = len(mus2)
                n_matched = len(dists)
                if n_matched > 0:
                    mean_shift = np.mean(np.abs(dists))
                    rmsd = np.sqrt(np.mean(dists**2))
                else:
                    mean_shift, rmsd = np.nan, np.nan
                match_rate = n_matched / n_ref if n_ref > 0 else 0.0

                results.append({
                    "file": file2,
                    "n_ref": n_ref,
                    "n_cmp": n_cmp,
                    "n_matched": n_matched,
                    "unmatched_ref": len(unref),
                    "unmatched_cmp": len(uncmp),
                    "match_rate": match_rate,
                    "mean_shift": mean_shift,
                    "rmsd": rmsd
                })

                plt.figure(figsize=(10,5))
                plt.plot(x1, y1, label='experiment', color='blue')
                plt.plot(xl, yl, label='simulation', color='orange')

                plt.scatter(mus1, y1[np.searchsorted(x1, mus1)], 
                            color='blue', marker='o', s=80, label='exp. peak')
                plt.scatter(mus2, y2[np.searchsorted(x2, mus2)], 
                            color='orange', marker='x', s=80, label='sim. peak')
        
                for i_ref, i_cmp in zip(mref, mcmp):
                    plt.plot([mus1[i_ref], mus2[i_cmp]],
                             [y1[np.searchsorted(x1, mus1[i_ref])],
                              y2[np.searchsorted(x2, mus2[i_cmp])]],
                             'k--', linewidth=0.8)

                plt.title(f"exp. vs sim. peak comparison\n{file2}")
                plt.xlabel("x")
                plt.ylabel("y")
                plt.legend()
                plt.tight_layout()
                plt.savefig(f"/mnt/data/SYC/xasfitting/full_LaCe_doublelayer/zif8/ce/Ce-L3-{os.path.basename(os.path.dirname(file2))}.png")
                plt.close()
                #plt.show()
                i+=1
            except Exception as e:
                print(f"Error processing {file2}: {e}")

df = pd.DataFrame(results)
out_csv = "/mnt/data/SYC/xasfitting/full_LaCe_doublelayer/zif8/Ce_L3_peak_fitting_results.csv"
df.to_csv(out_csv, index=False)
print(f"Results saved to {out_csv}")
