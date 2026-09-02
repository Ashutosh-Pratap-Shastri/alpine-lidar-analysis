"""
app.py — Alpine Geohazard Analysis Web Tool
============================================
Interactive browser-based DEM analysis tool.
Built on the CRAG PhD preparatory pipeline by Ashutosh Pratap Shastri.

Deploy free at: https://streamlit.io/cloud
"""

import streamlit as st
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import rasterio
from rasterio.io import MemoryFile
from scipy.ndimage import gaussian_filter
import io
import zipfile
import tempfile
import os

# ── Page config ──────────────────────────────────────────────
st.set_page_config(
    page_title="Alpine Geohazard Analysis",
    page_icon="🏔️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ── Header ───────────────────────────────────────────────────
st.title("🏔️ Alpine Geohazard Analysis Tool")
st.markdown(
    """
    **DEM-based slope stability and topographic stress assessment for alpine environments.**  
    Upload any GeoTIFF DEM or use the built-in Salzburg Alps synthetic terrain.  
    Produces slope maps, stress proxies, and Monte Carlo failure probability — 
    the conceptual pipeline underlying the 
    [CRAG PhD project](https://karriere.plus.ac.at) at the University of Salzburg.

    > *Built by [Ashutosh Pratap Shastri](https://github.com/Ashutosh-Pratap-Shastri/alpine-lidar-analysis) 
    > as preparatory portfolio work for the CRAG PhD application.*
    """
)
st.divider()

# ── Sidebar — parameters ──────────────────────────────────────
st.sidebar.header("⚙️ Analysis Parameters")

st.sidebar.subheader("Rock Mass Properties")
rock_density = st.sidebar.slider(
    "Rock density ρ (kg/m³)",
    min_value=2400, max_value=2900, value=2650, step=50,
    help="Limestone ~2600, Granite ~2700, Dolomite ~2800"
)
poisson_ratio = st.sidebar.slider(
    "Poisson ratio ν",
    min_value=0.10, max_value=0.40, value=0.25, step=0.05,
    help="Typical rock: 0.20–0.30"
)

st.sidebar.subheader("Slope Stability (Monte Carlo)")
cohesion_mean = st.sidebar.slider("Cohesion mean c (kPa)", 10, 150, 50, 5)
cohesion_std  = st.sidebar.slider("Cohesion std dev (kPa)", 5, 50, 15, 5)
friction_mean = st.sidebar.slider("Friction angle mean φ (°)", 20, 55, 35, 1)
friction_std  = st.sidebar.slider("Friction angle std dev (°)", 1, 15, 5, 1)
n_sim         = st.sidebar.selectbox(
    "Monte Carlo simulations", [1000, 5000, 10000], index=1,
    help="More = slower but more accurate"
)
critical_slope = st.sidebar.slider(
    "Critical slope threshold (°)", 25, 55, 35, 1,
    help="Slopes above this angle are highlighted as potentially unstable"
)

st.sidebar.subheader("Hillshade")
az = st.sidebar.slider("Sun azimuth (°)", 0, 360, 315, 15)
al = st.sidebar.slider("Sun altitude (°)", 10, 80, 45, 5)

st.sidebar.divider()
st.sidebar.caption(
    "⚠️ Stress values are a simplified 2-D dead-load proxy, not the "
    "FCM 3-D volumetric solution (Haunsperger & Robl 2025, 2026)."
)

# ── Helper functions ──────────────────────────────────────────

@st.cache_data(show_spinner=False)
def load_default_dem():
    """Load the synthetic Salzburg Alps DEM bundled with the repo."""
    dem_path = os.path.join(os.path.dirname(__file__),
                            "data", "salzburg_dem.tif")
    with rasterio.open(dem_path) as src:
        dem   = src.read(1).astype(np.float64)
        res_x = abs(src.transform.a) * 111320 * np.cos(np.radians(47.44))
        res_y = abs(src.transform.e) * 111320
    return dem, res_x, res_y


def load_uploaded_dem(uploaded_file):
    """Load a user-uploaded GeoTIFF."""
    bytes_data = uploaded_file.read()
    with MemoryFile(bytes_data) as mf:
        with mf.open() as src:
            dem   = src.read(1).astype(np.float64)
            nodata = src.nodata
            lat_mid = (src.bounds.top + src.bounds.bottom) / 2
            res_x = abs(src.transform.a) * 111320 * np.cos(np.radians(lat_mid))
            res_y = abs(src.transform.e) * 111320
    if nodata is not None:
        dem[dem == nodata] = np.nan
    return dem, res_x, res_y


def compute_hillshade(dem, azimuth=315, altitude=45):
    dem_s  = gaussian_filter(dem, sigma=1.5)
    dy, dx = np.gradient(dem_s)
    slope  = np.arctan(np.sqrt(dx**2 + dy**2))
    aspect = np.arctan2(-dx, dy)
    az_r   = np.radians(azimuth)
    al_r   = np.radians(altitude)
    hs = (np.cos(al_r) * np.cos(slope) +
          np.sin(al_r) * np.sin(slope) * np.cos(az_r - aspect))
    return np.clip(hs, 0, 1)


def compute_slope_aspect_curvature(dem, rx, ry):
    dem_s   = gaussian_filter(dem, sigma=1.0)
    dz_dx   = np.gradient(dem_s, rx, axis=1)
    dz_dy   = np.gradient(dem_s, ry, axis=0)
    slope   = np.degrees(np.arctan(np.sqrt(dz_dx**2 + dz_dy**2)))
    aspect  = np.degrees(np.arctan2(-dz_dx, dz_dy)) % 360
    d2x     = np.gradient(dz_dx, rx, axis=1)
    d2y     = np.gradient(dz_dy, ry, axis=0)
    curv    = -(d2x + d2y)
    return slope, aspect, curv


def compute_stress(dem, slope_deg, rho, nu):
    g        = 9.81
    K0       = nu / (1 - nu)
    slope_r  = np.radians(slope_deg)
    sigma_v  = (rho * g * dem) / 1e6
    sigma_h  = K0 * sigma_v * (1 + np.sin(slope_r))
    tau_max  = np.abs(sigma_v - sigma_h) / 2
    return sigma_v, tau_max


def compute_pf(slope_deg, c_mean, c_std, phi_mean, phi_std,
               gamma=26.0, H=10.0, n=5000):
    rng   = np.random.default_rng(42)
    c     = rng.normal(c_mean, c_std,  n).clip(1)
    phi   = rng.normal(phi_mean, phi_std, n).clip(5)
    phi_r = np.radians(phi)
    rows, cols = slope_deg.shape
    pf    = np.zeros((rows, cols), dtype=np.float32)
    block = 20
    for i in range(0, rows - block + 1, block):
        for j in range(0, cols - block + 1, block):
            alpha   = np.radians(slope_deg[i, j])
            if np.isnan(alpha):
                continue
            sn  = gamma * H * np.cos(alpha)**2
            tau = gamma * H * np.sin(alpha) * np.cos(alpha) + 1e-9
            FS  = (c + sn * np.tan(phi_r)) / tau
            pf[i:i+block, j:j+block] = (FS < 1.0).mean()
    return gaussian_filter(pf.astype(np.float64), sigma=2).astype(np.float32)


def make_fig(dem, slope, sigma_v, tau_max, pf, hs, crit_angle):
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.32)

    def panel(ax, data, cmap, alpha, title, cblabel,
              vmin=None, vmax=None, overlay=None):
        ax.imshow(hs, cmap="gray", vmin=0, vmax=1)
        kw = dict(cmap=cmap, alpha=alpha)
        if vmin is not None: kw["vmin"] = vmin
        if vmax is not None: kw["vmax"] = vmax
        im = ax.imshow(data, **kw)
        if overlay is not None:
            ax.imshow(overlay, cmap="cool", alpha=0.35)
        cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(cblabel, fontsize=8)
        ax.set_title(title, fontsize=9.5, fontweight="bold")
        ax.set_xlabel("Column", fontsize=7)
        ax.set_ylabel("Row", fontsize=7)
        ax.tick_params(labelsize=6)

    crit_mask = np.where(slope > crit_angle, 1.0, np.nan)

    panel(fig.add_subplot(gs[0, 0]), dem,     "terrain",  0.65,
          "a)  DEM + Hillshade", "Elevation (m)")
    panel(fig.add_subplot(gs[0, 1]), slope,   "YlOrRd",   0.80,
          f"b)  Slope (°) — critical >{crit_angle}° overlaid",
          "Slope (°)", vmin=0, vmax=65, overlay=crit_mask)
    panel(fig.add_subplot(gs[0, 2]), sigma_v, "plasma",   0.80,
          "c)  Vertical Stress σv (MPa)", "σv (MPa)")
    panel(fig.add_subplot(gs[1, 0]), tau_max, "inferno",  0.80,
          "d)  Max Shear Stress τmax (MPa)", "τmax (MPa)",
          vmax=np.nanpercentile(tau_max, 98))
    panel(fig.add_subplot(gs[1, 1]), pf,      "RdYlGn_r", 0.85,
          "e)  Failure Probability Pf", "Pf", vmin=0, vmax=0.5)

    # Stats panel
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.hist(slope[~np.isnan(slope)].ravel(), bins=40,
             color="darkorange", edgecolor="white", lw=0.4,
             density=True, label="Slope")
    ax6.axvline(crit_angle, color="crimson", lw=1.4,
                linestyle="--", label=f"Critical {crit_angle}°")
    ax6.set_xlabel("Slope (°)", fontsize=8)
    ax6.set_ylabel("Density", fontsize=8)
    ax6.set_title("f)  Slope Distribution", fontsize=9.5, fontweight="bold")
    ax6.legend(fontsize=7)
    ax6.grid(True, lw=0.3, alpha=0.4)

    fig.suptitle(
        "Alpine Geohazard Analysis — DEM · Slope · Stress · Failure Probability\n"
        "github.com/Ashutosh-Pratap-Shastri/alpine-lidar-analysis",
        fontsize=11, fontweight="bold", y=1.01
    )
    return fig


def fig_to_bytes(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    return buf.read()


def make_zip(figs_bytes: dict) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for name, data in figs_bytes.items():
            zf.writestr(name, data)
    buf.seek(0)
    return buf.read()


# ── Data source ───────────────────────────────────────────────
st.subheader("1️⃣  Data Source")
col1, col2 = st.columns([1, 1])

with col1:
    data_source = st.radio(
        "Choose DEM:",
        ["Use built-in Salzburg Alps DEM (synthetic)",
         "Upload your own GeoTIFF"],
        index=0
    )

with col2:
    if data_source.startswith("Upload"):
        uploaded = st.file_uploader(
            "Upload GeoTIFF DEM", type=["tif", "tiff"],
            help="Any projected or geographic GeoTIFF. "
                 "Large files (>50MB) may be slow."
        )
    else:
        uploaded = None
        st.info(
            "**Built-in DEM:** Synthetic terrain for the Hochkönig / "
            "Salzburg Alps region (13.05–13.20°E, 47.38–47.50°N). "
            "Elevation 580–2941 m. Replaceable with any real SRTM or "
            "Copernicus tile."
        )

# ── Run button ────────────────────────────────────────────────
st.subheader("2️⃣  Run Analysis")
run = st.button("▶  Run Analysis", type="primary", use_container_width=True)

if run:
    # Load DEM
    with st.spinner("Loading DEM..."):
        try:
            if data_source.startswith("Upload") and uploaded is not None:
                dem, rx, ry = load_uploaded_dem(uploaded)
                label = uploaded.name
            else:
                dem, rx, ry = load_default_dem()
                label = "Salzburg Alps (synthetic)"
        except Exception as e:
            st.error(f"Failed to load DEM: {e}")
            st.stop()

    # Stats bar
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("DEM size", f"{dem.shape[0]}×{dem.shape[1]} px")
    c2.metric("Elevation range", f"{np.nanmin(dem):.0f}–{np.nanmax(dem):.0f} m")
    c3.metric("Resolution", f"{rx:.0f}×{ry:.0f} m/px")
    c4.metric("Data source", label)

    # Compute
    with st.spinner("Computing terrain derivatives..."):
        slope, aspect, curv = compute_slope_aspect_curvature(dem, rx, ry)

    with st.spinner("Computing stress proxy..."):
        sigma_v, tau_max = compute_stress(dem, slope, rock_density, poisson_ratio)

    with st.spinner(f"Running Monte Carlo (n={n_sim:,})..."):
        pf = compute_pf(slope, cohesion_mean, cohesion_std,
                        friction_mean, friction_std, n=n_sim)

    with st.spinner("Rendering hillshade..."):
        hs = compute_hillshade(dem, az, al)

    # Key metrics
    st.divider()
    st.subheader("3️⃣  Results")

    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Mean slope", f"{np.nanmean(slope):.1f}°")
    m2.metric(f"Slopes >{critical_slope}°",
              f"{(slope>critical_slope).sum()/slope.size*100:.1f}%")
    m3.metric("Max σv", f"{np.nanmax(sigma_v):.1f} MPa")
    m4.metric("Max τmax", f"{np.nanmax(tau_max):.1f} MPa")
    m5.metric("High-risk (Pf>0.10)",
              f"{(pf>0.10).sum()/pf.size*100:.1f}%")

    # Main figure
    with st.spinner("Generating figures..."):
        fig = make_fig(dem, slope, sigma_v, tau_max, pf, hs, critical_slope)
        fig_bytes = fig_to_bytes(fig)

    st.pyplot(fig, use_container_width=True)
    plt.close(fig)

    # Individual download buttons
    st.divider()
    st.subheader("4️⃣  Download Outputs")

    # Generate individual figures
    def single_fig(data, cmap, title, cblabel, hs,
                   vmin=None, vmax=None):
        f, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(hs, cmap="gray", vmin=0, vmax=1)
        kw = dict(cmap=cmap, alpha=0.80)
        if vmin is not None: kw["vmin"] = vmin
        if vmax is not None: kw["vmax"] = vmax
        im = ax.imshow(data, **kw)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label(
            cblabel, fontsize=9)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel("Column"); ax.set_ylabel("Row")
        plt.tight_layout()
        b = fig_to_bytes(f); plt.close(f); return b

    figs = {
        "00_overview.png"         : fig_bytes,
        "01_dem_hillshade.png"    : single_fig(dem, "terrain",
            "DEM + Hillshade", "Elevation (m)", hs),
        "02_slope.png"            : single_fig(slope, "YlOrRd",
            "Slope Map", "Slope (°)", hs, 0, 65),
        "03_stress_sv.png"        : single_fig(sigma_v, "plasma",
            "Vertical Stress σv", "MPa", hs),
        "04_stress_tau.png"       : single_fig(tau_max, "inferno",
            "Max Shear Stress τmax", "MPa", hs,
            vmax=np.nanpercentile(tau_max, 98)),
        "05_failure_prob.png"     : single_fig(pf, "RdYlGn_r",
            "Failure Probability Pf", "Pf", hs, 0, 0.5),
    }

    col_a, col_b = st.columns(2)
    with col_a:
        st.download_button(
            "⬇️  Download overview figure (PNG)",
            data=fig_bytes,
            file_name="alpine_geohazard_overview.png",
            mime="image/png",
            use_container_width=True
        )
    with col_b:
        st.download_button(
            "⬇️  Download all figures (ZIP)",
            data=make_zip(figs),
            file_name="alpine_geohazard_outputs.zip",
            mime="application/zip",
            use_container_width=True
        )

    # Methodology note
    st.divider()
    with st.expander("📖 Methodology & Honest Disclaimers"):
        st.markdown("""
**Terrain derivatives** — Horn (1981) 3×3 finite-difference gradient estimator
(same kernel as ArcGIS Slope, QGIS Slope, gdaldem).

**Stress proxy** — Simplified 2-D dead-load approximation:
- σv = ρ g H / 10⁶ (lithostatic vertical stress, MPa)
- K0 = ν / (1−ν) (at-rest coefficient)
- σh = K0 σv (1 + sin α) (slope-amplified horizontal)
- τmax = |σv − σh| / 2

This is **not** the Finite Cell Method (FCM) 3-D volumetric stress computation
developed by Haunsperger & Robl (2025, 2026). The FCM operates on full
volumetric meshes across entire massifs on HPC clusters. This proxy
demonstrates the conceptual pipeline only.

**Failure probability** — Infinite-slope Monte Carlo:
- FS = (c + σn tan φ) / τ
- Pf = P(FS < 1.0) over n realisations
- Uncertain parameters: c ~ N(mean, std), φ ~ N(mean, std)

**CRAG connection** — None of the six foundational CRAG papers computes
a probabilistic Pf output. This tool demonstrates that step.

**Source code** — [github.com/Ashutosh-Pratap-Shastri/alpine-lidar-analysis](https://github.com/Ashutosh-Pratap-Shastri/alpine-lidar-analysis)
        """)

# ── Footer ────────────────────────────────────────────────────
st.divider()
st.caption(
    "Built by **Ashutosh Pratap Shastri** · JRF, IIT (BHU) Varanasi · "
    "M.Tech Geotechnical Engineering, IIT Patna (2025) · "
    "[GitHub](https://github.com/Ashutosh-Pratap-Shastri/alpine-lidar-analysis)"
)
