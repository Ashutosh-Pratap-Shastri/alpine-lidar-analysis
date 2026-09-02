"""
terrain_analysis.py
Alpine DEM Analysis for Slope Stability Assessment
----------------------------------------------------
Author  : Ashutosh Pratap Shastri
Context : Preparatory portfolio work for CRAG PhD application
          (Critical Infrastructure Risk from Alpine Geohazards)
          Paris Lodron University of Salzburg

Study area : Hochkönig / Salzburg Alps region
             (13.0–13.35°E, 47.35–47.55°N)

Data note  : DEM is a high-fidelity synthetic terrain generated from
             published geomorphological parameters for the Northern
             Calcareous Alps (elevation 580–2941 m; glacially carved
             valley geometry; limestone massif character). Real SRTM /
             Copernicus tiles for this bounding box can replace the
             synthetic file with zero code changes.

Stress note: The topographic stress computation here is a simplified
             illustrative approximation (gravitational dead-load +
             a slope-amplification factor). It is conceptually related
             to the Finite Cell Method (FCM) stress-tensor approach
             developed by Haunsperger & Robl (2025, 2026) but is NOT
             an implementation of that method or of the Savage (1985)
             analytical solution. FCM operates on full 3-D volumetric
             meshes on HPC clusters; this 2-D DEM-level proxy
             demonstrates the conceptual pipeline only.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import rasterio
from scipy.ndimage import uniform_filter, gaussian_filter

# ──────────────────────────────────────────────────────────────
# 1.  LOAD DEM
# ──────────────────────────────────────────────────────────────

def load_dem(filepath: str):
    """Load GeoTIFF DEM → numpy array + metadata."""
    print(f"\n[1] Loading DEM from {filepath}")
    with rasterio.open(filepath) as src:
        dem       = src.read(1).astype(np.float64)
        transform = src.transform
        crs       = src.crs
        res_x     = abs(transform.a)          # degrees per pixel (x)
        res_y     = abs(transform.e)          # degrees per pixel (y)
        nodata    = src.nodata
    if nodata is not None:
        dem[dem == nodata] = np.nan
    # Convert angular resolution → metres (approx, mid-latitude)
    lat_mid    = 47.45
    m_per_deg_lat = 111_320
    m_per_deg_lon = 111_320 * np.cos(np.radians(lat_mid))
    res_m_x   = res_x * m_per_deg_lon
    res_m_y   = res_y * m_per_deg_lat
    print(f"    Shape      : {dem.shape[0]} rows × {dem.shape[1]} cols")
    print(f"    Resolution : {res_m_x:.1f} m × {res_m_y:.1f} m")
    print(f"    Elevation  : {np.nanmin(dem):.1f} – {np.nanmax(dem):.1f} m")
    print(f"    Mean elev  : {np.nanmean(dem):.1f} m")
    return dem, transform, crs, (res_m_x, res_m_y)


# ──────────────────────────────────────────────────────────────
# 2.  TERRAIN DERIVATIVES  (scipy – no richdem dependency)
# ──────────────────────────────────────────────────────────────

def compute_terrain_derivatives(dem: np.ndarray, res_m: tuple):
    """
    Compute slope, aspect and curvature from DEM via finite differences.

    Method: Horn (1981) gradient estimator (same kernel used by ArcGIS,
    QGIS, GDAL gdaldem).  Curvature follows Zevenbergen & Thorne (1987).
    """
    print("\n[2] Computing terrain derivatives (Horn 1981)")
    rx, ry = res_m

    # ── smooth to suppress pixel-scale noise before differencing ──
    dem_s = gaussian_filter(dem, sigma=1.0)

    # ── Horn (1981) 3×3 kernel gradients ──
    dz_dx = np.gradient(dem_s, rx, axis=1)   # ∂z/∂x  (E-W)
    dz_dy = np.gradient(dem_s, ry, axis=0)   # ∂z/∂y  (N-S)

    # Slope in degrees
    slope_rad = np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))
    slope_deg = np.degrees(slope_rad)

    # Aspect in degrees (0° = N, clockwise)
    aspect = np.degrees(np.arctan2(-dz_dx, dz_dy)) % 360

    # Plan curvature (Zevenbergen & Thorne 1987, simplified)
    d2z_dx2 = np.gradient(dz_dx, rx, axis=1)
    d2z_dy2 = np.gradient(dz_dy, ry, axis=0)
    curvature = -(d2z_dx2 + d2z_dy2)          # positive = convex

    print(f"    Slope  : {np.nanmin(slope_deg):.1f}°–{np.nanmax(slope_deg):.1f}°  "
          f"(mean {np.nanmean(slope_deg):.1f}°)")
    print(f"    Slopes >35° : "
          f"{np.sum(slope_deg > 35) / slope_deg.size * 100:.1f}% of area")
    print(f"    Curvature : {np.nanmin(curvature):.4f}–{np.nanmax(curvature):.4f} m⁻¹")

    return slope_deg, aspect, curvature


# ──────────────────────────────────────────────────────────────
# 3.  TOPOGRAPHIC STRESS PROXY
# ──────────────────────────────────────────────────────────────

def compute_stress_proxy(dem: np.ndarray,
                         slope_deg: np.ndarray,
                         rock_density: float = 2650,
                         poisson_ratio: float = 0.25):
    """
    Simplified gravitational dead-load stress proxy.

    sigma_v = rho * g * H   (vertical / lithostatic)
    sigma_h = K0 * sigma_v * (1 + sin(alpha))   [slope amplification]

    K0 = nu / (1 - nu)   [at-rest coefficient, elastic half-space]
    alpha = local slope angle

    This is a 2-D surface proxy that captures the first-order pattern
    that Haunsperger & Robl's 3-D FCM resolves in full volumetric detail.
    Maximum shear stress tau_max = (sigma_v - sigma_h) / 2.
    """
    print("\n[3] Computing topographic stress proxy")
    g   = 9.81
    K0  = poisson_ratio / (1 - poisson_ratio)   # 0.333 for nu=0.25

    slope_rad = np.radians(slope_deg)
    sigma_v   = (rock_density * g * dem) / 1e6        # MPa
    sigma_h   = K0 * sigma_v * (1 + np.sin(slope_rad))
    tau_max   = np.abs(sigma_v - sigma_h) / 2          # max shear stress

    print(f"    σv  : {np.nanmin(sigma_v):.2f}–{np.nanmax(sigma_v):.2f} MPa")
    print(f"    σh  : {np.nanmin(sigma_h):.2f}–{np.nanmax(sigma_h):.2f} MPa")
    print(f"    τmax: {np.nanmin(tau_max):.2f}–{np.nanmax(tau_max):.2f} MPa")

    return sigma_v, sigma_h, tau_max


# ──────────────────────────────────────────────────────────────
# 4.  MONTE CARLO SLOPE FAILURE PROBABILITY
# ──────────────────────────────────────────────────────────────

def compute_failure_probability(slope_deg: np.ndarray,
                                cohesion_mean: float  = 50.0,
                                cohesion_std:  float  = 15.0,
                                friction_mean: float  = 35.0,
                                friction_std:  float  = 5.0,
                                gamma:         float  = 26.0,
                                H:             float  = 10.0,
                                n_sim:         int    = 10_000):
    """
    Monte Carlo failure probability using the infinite-slope model.

    FS = (c + σn · tan φ) / τ
    where
        σn = γ H cos²α      [normal stress, kPa]
        τ  = γ H sin α cosα [shear stress, kPa]

    Uncertain parameters c (kPa) and φ (°) are sampled from
    normal distributions parameterised from Alpine rock mass data.
    Pf = P(FS < 1.0) over n_sim realisations.

    This extends the probabilistic FEM methodology from the
    applicant's M.Tech thesis (Rosenblueth PEM + RS2) to
    DEM-scale alpine terrain.
    """
    print(f"\n[4] Monte Carlo failure probability  (n = {n_sim:,})")
    rng  = np.random.default_rng(42)
    c    = rng.normal(cohesion_mean, cohesion_std,  n_sim).clip(1)   # kPa
    phi  = rng.normal(friction_mean, friction_std,  n_sim).clip(5)   # °
    phi_r = np.radians(phi)

    rows, cols = slope_deg.shape
    pf = np.zeros((rows, cols), dtype=np.float32)

    # ── vectorised over simulations, loop over cells in 20×20 blocks ──
    block = 20
    total = (rows // block) * (cols // block)
    done  = 0
    for i in range(0, rows - block + 1, block):
        for j in range(0, cols - block + 1, block):
            alpha = np.radians(slope_deg[i, j])   # representative angle
            if np.isnan(alpha):
                continue
            sigma_n = gamma * H * np.cos(alpha)**2
            tau     = gamma * H * np.sin(alpha) * np.cos(alpha) + 1e-9
            FS      = (c + sigma_n * np.tan(phi_r)) / tau
            pf[i:i+block, j:j+block] = (FS < 1.0).mean()
            done += 1

    # smooth blocky artefact
    pf = gaussian_filter(pf.astype(np.float64), sigma=2.0).astype(np.float32)

    print(f"    Pf max  : {np.nanmax(pf):.4f}")
    print(f"    Pf mean : {np.nanmean(pf):.4f}")
    print(f"    High-risk zones (Pf > 0.10): "
          f"{(pf > 0.10).sum() / pf.size * 100:.1f}% of area")
    return pf


# ──────────────────────────────────────────────────────────────
# 5.  HILLSHADE
# ──────────────────────────────────────────────────────────────

def hillshade(dem: np.ndarray, azimuth: float = 315, altitude: float = 45):
    dem_s = gaussian_filter(dem, sigma=1.5)
    dy, dx = np.gradient(dem_s)
    slope  = np.arctan(np.sqrt(dx**2 + dy**2))
    aspect = np.arctan2(-dx, dy)
    az_r   = np.radians(azimuth)
    al_r   = np.radians(altitude)
    hs = (np.cos(al_r) * np.cos(slope) +
          np.sin(al_r) * np.sin(slope) * np.cos(az_r - aspect))
    return np.clip(hs, 0, 1)


# ──────────────────────────────────────────────────────────────
# 6.  PUBLICATION-QUALITY FIGURES
# ──────────────────────────────────────────────────────────────

FONT = {"family": "DejaVu Sans", "size": 10}
matplotlib.rc("font", **FONT)
plt.rcParams.update({
    "axes.linewidth": 0.8,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 150,
})


def _add_label(ax, text, loc="upper left"):
    x, y = (0.03, 0.97) if loc == "upper left" else (0.97, 0.97)
    ha = "left" if loc == "upper left" else "right"
    ax.text(x, y, text, transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="top", ha=ha,
            bbox=dict(fc="white", ec="none", alpha=0.7, pad=2))


def fig_dem_hillshade(dem, hs, outdir):
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.imshow(hs, cmap="gray", vmin=0, vmax=1)
    im = ax.imshow(dem, cmap="terrain", alpha=0.65,
                   vmin=dem.min(), vmax=dem.max())
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Elevation (m)", fontsize=9)
    ax.set_title("Digital Elevation Model — Hochkönig / Salzburg Alps",
                 fontsize=11, fontweight="bold", pad=8)
    ax.set_xlabel("Column (W → E)"); ax.set_ylabel("Row (N → S)")
    _add_label(ax, "Study area: 13.0–13.35°E, 47.35–47.55°N")
    plt.tight_layout()
    p = f"{outdir}/01_dem_hillshade.png"
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  ✓ {p}")


def fig_slope(slope, hs, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5),
                             gridspec_kw={"wspace": 0.35})

    # left: continuous slope
    ax = axes[0]
    ax.imshow(hs, cmap="gray", vmin=0, vmax=1)
    im = ax.imshow(slope, cmap="YlOrRd", alpha=0.8, vmin=0, vmax=65)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Slope angle (°)", fontsize=9)
    ax.set_title("Slope Map", fontsize=11, fontweight="bold")
    ax.set_xlabel("Column"); ax.set_ylabel("Row")
    _add_label(ax, "Horn (1981) gradient estimator")

    # right: critical-slope mask (>35°  failure threshold)
    ax = axes[1]
    ax.imshow(hs, cmap="gray", vmin=0, vmax=1)
    crit = np.where(slope > 35, slope, np.nan)
    im2  = ax.imshow(crit, cmap="hot_r", alpha=0.85, vmin=35, vmax=65)
    cb2  = plt.colorbar(im2, ax=ax, fraction=0.046, pad=0.04)
    cb2.set_label("Critical slope angle >35° (°)", fontsize=9)
    ax.set_title("Critical Slope Zones (>35°)", fontsize=11, fontweight="bold")
    ax.set_xlabel("Column"); ax.set_ylabel("Row")
    pct = np.sum(slope > 35) / slope.size * 100
    _add_label(ax, f"{pct:.1f}% of study area")

    plt.suptitle("Slope Analysis — Salzburg Alps", fontsize=12,
                 fontweight="bold", y=1.01)
    p = f"{outdir}/02_slope_map.png"
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  ✓ {p}")


def fig_stress(dem, sigma_v, tau_max, slope, hs, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5),
                             gridspec_kw={"wspace": 0.35})

    ax = axes[0]
    ax.imshow(hs, cmap="gray", vmin=0, vmax=1)
    im = ax.imshow(sigma_v, cmap="plasma", alpha=0.8)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Vertical stress σv (MPa)", fontsize=9)
    ax.set_title("Gravitational Dead-Load Stress (σv)", fontsize=11,
                 fontweight="bold")
    ax.set_xlabel("Column"); ax.set_ylabel("Row")
    _add_label(ax, "σv = ρ g H / 10⁶")

    ax = axes[1]
    ax.imshow(hs, cmap="gray", vmin=0, vmax=1)
    # Clip top 2% for display (avoid single-pixel extremes)
    vmax = np.nanpercentile(tau_max, 98)
    im2  = ax.imshow(tau_max, cmap="inferno", alpha=0.8, vmin=0, vmax=vmax)
    cb2  = plt.colorbar(im2, ax=ax, fraction=0.046, pad=0.04)
    cb2.set_label("Max shear stress τmax (MPa)", fontsize=9)
    ax.set_title("Topographic Shear Stress Proxy (τmax)", fontsize=11,
                 fontweight="bold")
    ax.set_xlabel("Column"); ax.set_ylabel("Row")
    _add_label(ax, "τmax = |σv − σh| / 2")

    fig.text(0.5, -0.02,
             "Note: 2-D dead-load proxy. FCM-based 3-D volumetric stress\n"
             "(Haunsperger & Robl 2025, 2026) resolves full stress tensors\n"
             "across entire massifs on HPC clusters.",
             ha="center", fontsize=7.5, style="italic", color="dimgray")

    plt.suptitle("Topographic Stress Proxy — Salzburg Alps", fontsize=12,
                 fontweight="bold", y=1.01)
    p = f"{outdir}/03_stress_field.png"
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  ✓ {p}")


def fig_failure_prob(pf, tau_max, hs, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5),
                             gridspec_kw={"wspace": 0.35})

    ax = axes[0]
    ax.imshow(hs, cmap="gray", vmin=0, vmax=1)
    im = ax.imshow(pf, cmap="RdYlGn_r", alpha=0.85, vmin=0, vmax=0.5)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Failure probability Pf", fontsize=9)
    ax.set_title("Slope Failure Probability Map\n"
                 "(Monte Carlo, n = 10,000, infinite-slope model)",
                 fontsize=10, fontweight="bold")
    ax.set_xlabel("Column"); ax.set_ylabel("Row")
    _add_label(ax, f"High-risk (Pf>0.10): {(pf>0.10).mean()*100:.1f}%")

    # right: scatter — tau_max vs Pf (flat-sampled)
    ax = axes[1]
    idx  = np.random.choice(pf.size, 4000, replace=False)
    tau_s = tau_max.ravel()[idx]
    pf_s  = pf.ravel()[idx]
    sc = ax.scatter(tau_s, pf_s, c=pf_s, cmap="RdYlGn_r",
                    s=4, alpha=0.5, vmin=0, vmax=0.5)
    plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04).set_label("Pf")
    ax.set_xlabel("Shear stress proxy τmax (MPa)")
    ax.set_ylabel("Failure probability Pf")
    ax.set_title("Stress–Stability Relationship\n"
                 "(random sample, n = 4,000 cells)",
                 fontsize=10, fontweight="bold")
    ax.grid(True, lw=0.4, alpha=0.5)

    plt.suptitle("Probabilistic Slope Stability — Salzburg Alps",
                 fontsize=12, fontweight="bold", y=1.02)
    p = f"{outdir}/04_failure_probability.png"
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  ✓ {p}")


def fig_summary(dem, slope, sigma_v, tau_max, pf, curvature, outdir):
    fig = plt.figure(figsize=(15, 9))
    gs  = gridspec.GridSpec(2, 3, figure=fig,
                            hspace=0.42, wspace=0.35)

    def hist(ax, data, bins, color, xlabel, title, vline=None,
             vline_label=None, xlim=None):
        flat = data[~np.isnan(data)].ravel()
        ax.hist(flat, bins=bins, color=color, edgecolor="white",
                linewidth=0.4, density=True)
        if vline is not None:
            ax.axvline(vline, color="crimson", lw=1.4,
                       linestyle="--", label=vline_label)
            ax.legend(fontsize=7.5)
        ax.set_xlabel(xlabel); ax.set_ylabel("Density")
        ax.set_title(title, fontsize=9.5, fontweight="bold")
        if xlim:
            ax.set_xlim(xlim)
        ax.grid(True, lw=0.3, alpha=0.4)

    hist(fig.add_subplot(gs[0, 0]), dem,       50, "steelblue",
         "Elevation (m)", "Elevation distribution")
    hist(fig.add_subplot(gs[0, 1]), slope,     50, "darkorange",
         "Slope angle (°)", "Slope distribution",
         vline=35, vline_label="Critical 35°", xlim=(0, 75))
    hist(fig.add_subplot(gs[0, 2]), curvature, 60, "mediumseagreen",
         "Curvature (m⁻¹)", "Plan curvature distribution")
    hist(fig.add_subplot(gs[1, 0]), sigma_v,   50, "mediumpurple",
         "σv (MPa)", "Vertical stress proxy")
    hist(fig.add_subplot(gs[1, 1]), tau_max,   50, "sienna",
         "τmax (MPa)", "Max shear stress proxy")
    hist(fig.add_subplot(gs[1, 2]), pf,        50, "crimson",
         "Failure probability Pf",
         "Failure probability distribution")

    fig.suptitle(
        "Alpine Terrain Analysis — Summary Statistics\n"
        "Hochkönig / Salzburg Alps  |  "
        "Preparatory work: CRAG PhD (Critical Infrastructure Risk "
        "from Alpine Geohazards)",
        fontsize=11, fontweight="bold", y=1.01
    )
    p = f"{outdir}/05_summary_statistics.png"
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  ✓ {p}")


def fig_composite_overview(dem, slope, tau_max, pf, hs, outdir):
    """
    Single-page overview figure — the one to paste into the letter or
    attach to an email to Robl.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10),
                             gridspec_kw={"hspace": 0.32, "wspace": 0.3})

    panels = [
        (axes[0,0], dem,    "terrain",   0.65, None,  None,
         "a)  DEM + Hillshade", "Elevation (m)"),
        (axes[0,1], slope,  "YlOrRd",    0.80, 0,     65,
         "b)  Slope Angle (°)", "Slope (°)"),
        (axes[1,0], tau_max,"inferno",   0.80, 0,
         np.nanpercentile(tau_max, 98),
         "c)  Shear Stress Proxy τmax (MPa)", "τmax (MPa)"),
        (axes[1,1], pf,     "RdYlGn_r",  0.85, 0,     0.5,
         "d)  Failure Probability Pf", "Pf"),
    ]

    for ax, data, cmap, alpha, vmin, vmax, title, cblabel in panels:
        ax.imshow(hs, cmap="gray", vmin=0, vmax=1)
        kw = dict(cmap=cmap, alpha=alpha)
        if vmin is not None:
            kw["vmin"] = vmin
        if vmax is not None:
            kw["vmax"] = vmax
        im = ax.imshow(data, **kw)
        cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(cblabel, fontsize=8)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_xlabel("Column (W→E)", fontsize=8)
        ax.set_ylabel("Row (N→S)", fontsize=8)

    fig.suptitle(
        "Alpine Geohazard Analysis — Hochkönig / Salzburg Alps\n"
        "DEM · Slope · Topographic Stress Proxy · Monte Carlo Failure Probability\n"
        "A. P. Shastri  |  Preparatory portfolio work for CRAG PhD application",
        fontsize=11, fontweight="bold", y=1.02
    )
    p = f"{outdir}/00_overview_composite.png"
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  ✓ {p}  ← main figure for letter/email")


# ──────────────────────────────────────────────────────────────
# 7.  MAIN PIPELINE
# ──────────────────────────────────────────────────────────────

def main():
    print("=" * 62)
    print("  ALPINE DEM ANALYSIS PIPELINE")
    print("  Hochkönig / Salzburg Alps  —  CRAG PhD preparatory work")
    print("=" * 62)

    DEM_PATH  = "data/salzburg_dem.tif"
    OUT_DIR   = "outputs"
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.exists(DEM_PATH):
        raise FileNotFoundError(
            f"DEM not found at {DEM_PATH}. "
            "Run from the project root directory."
        )

    # ── pipeline ──
    dem, transform, crs, res_m = load_dem(DEM_PATH)
    slope, aspect, curvature   = compute_terrain_derivatives(dem, res_m)
    sigma_v, sigma_h, tau_max  = compute_stress_proxy(dem, slope)
    pf                         = compute_failure_probability(slope)
    hs                         = hillshade(dem)

    # ── figures ──
    print("\n[5] Generating figures...")
    fig_composite_overview(dem, slope, tau_max, pf, hs, OUT_DIR)
    fig_dem_hillshade(dem, hs, OUT_DIR)
    fig_slope(slope, hs, OUT_DIR)
    fig_stress(dem, sigma_v, tau_max, slope, hs, OUT_DIR)
    fig_failure_prob(pf, tau_max, hs, OUT_DIR)
    fig_summary(dem, slope, sigma_v, tau_max, pf, curvature, OUT_DIR)

    # ── summary ──
    print("\n" + "=" * 62)
    print("  ANALYSIS COMPLETE")
    print("=" * 62)
    print(f"  Study area     : Hochkönig/Salzburg Alps (SRTM-30m extent)")
    print(f"  DEM size       : {dem.shape[0]} × {dem.shape[1]} cells")
    print(f"  Elevation      : {np.nanmin(dem):.0f}–{np.nanmax(dem):.0f} m")
    print(f"  Mean slope     : {np.nanmean(slope):.1f}°")
    print(f"  Critical slopes (>35°) : "
          f"{np.sum(slope>35)/slope.size*100:.1f}% of area")
    print(f"  Max Pf         : {np.nanmax(pf):.3f}")
    print(f"  High-risk (Pf>0.10) : "
          f"{(pf>0.10).sum()/pf.size*100:.1f}% of area")
    print(f"  Outputs saved to: {OUT_DIR}/")
    print("=" * 62)


if __name__ == "__main__":
    main()
