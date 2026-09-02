# Alpine DEM Analysis — Hochkönig / Salzburg Alps

Terrain analysis pipeline for slope stability and topographic stress
assessment in high-alpine environments.

---

## Context

Preparatory portfolio work for the **CRAG PhD position**
(*Critical Infrastructure Risk from Alpine Geohazards*)
at the Paris Lodron University of Salzburg — Department of Environment
and Biodiversity, DSP DynamitE Vol. 2.

The pipeline demonstrates familiarity with DEM-based terrain analysis
as the foundational step in the multiscale monitoring workflow
(TLS → UAV → InSAR → geomechanical model) that CRAG uses to predict
gravitational mass movements threatening critical infrastructure.

---

## What This Does

| Step | Module | Output |
|------|--------|--------|
| 1 | `load_dem` | Elevation array + georeferencing metadata |
| 2 | `compute_terrain_derivatives` | Slope, aspect, plan curvature (Horn 1981) |
| 3 | `compute_stress_proxy` | Gravitational dead-load stress proxy (σv, σh, τmax) |
| 4 | `compute_failure_probability` | Monte Carlo Pf map (n = 10,000, infinite-slope model) |
| 5 | Figures | 6 publication-quality PNG outputs |

---

## Study Area

**Hochkönig / Salzburg Alps region**
Bounding box: 13.05–13.20°E, 47.38–47.50°N
Resolution: ~18–22 m per pixel

> **Data note:** The DEM in `data/salzburg_dem.tif` is a high-fidelity
> synthetic terrain generated from published geomorphological parameters
> for the Northern Calcareous Alps (elevation 580–2941 m; glacially
> carved valley geometry; 20.1% of area with slopes >35°; limestone
> massif character consistent with Hochkönig).
> Any real SRTM-30m or Copernicus GLO-30 tile for this bounding box
> can replace it with zero code changes.

---

## Key Results

```
Elevation range   :  580 – 2,941 m
Mean slope        :  24.2°
Critical slopes (>35°) :  20.1% of area
Max failure probability Pf :  0.795
High-risk zones (Pf > 0.10) :  10.6% of area
Vertical stress proxy σv :  15 – 76 MPa
Max shear stress proxy τmax :  up to 25.5 MPa
```

---

## Methodological Notes

### Terrain derivatives
Slope, aspect and curvature are computed via the Horn (1981) 3×3
finite-difference gradient estimator — the same kernel used by
ArcGIS `Slope`, QGIS `Slope`, and `gdaldem slope`.

### Stress proxy
The topographic dead-load stress proxy here is a simplified 2-D
approximation:

```
σv  = ρ g H / 10⁶          (lithostatic vertical stress, MPa)
K0  = ν / (1 − ν)           (at-rest coefficient, elastic half-space)
σh  = K0 σv (1 + sin α)     (horizontal with slope amplification)
τmax = |σv − σh| / 2        (maximum shear stress)
```

This is **not** an implementation of the Finite Cell Method (FCM)
3-D stress computation developed by Haunsperger & Robl (2025, 2026),
nor of the Savage (1985) analytical solution. The FCM operates on full
volumetric meshes across entire massifs on HPC clusters and resolves
complete 3-D stress tensors. This proxy demonstrates the conceptual
pipeline only and produces the same first-order spatial pattern.

### Failure probability
The infinite-slope Monte Carlo treats cohesion c and friction angle φ
as normally distributed random variables:

```
c  ~ N(50, 15) kPa
φ  ~ N(35,  5) °
FS = (c + σn tan φ) / τ
Pf = P(FS < 1.0),  n = 10,000 realisations
```

This extends the probabilistic FEM methodology from the author's
M.Tech thesis (Rosenblueth PEM + Monte Carlo in RS2) to
DEM-scale alpine terrain.

---

## Connection to CRAG

The FCM-based 3-D stress modelling in CRAG
(primary supervisor: Assoc.-Prof. Jörg Robl, Univ. Salzburg;
co-supervisor: Prof. Barbara Schneider-Muntau, Univ. Innsbruck)
operates on the same conceptual chain — topography → stress field →
failure threshold — but resolved at full massif scale with volumetric
finite cell meshes on HPC clusters, coupled with real multiscale
monitoring data (InSAR, UAV, TLS, fissurometers) from Jan-Christoph
Otto's alpine field sites.

The identified gap this pipeline illustrates:
> *"None of the six foundational papers in CRAG's literature base
> computes a probabilistic failure output — Haunsperger & Robl produce
> deterministic stress tensors; Schneider-Muntau produces deterministic
> factors of safety. Propagating geological uncertainty through the
> stress and stability model to obtain Pf is the step this pipeline
> demonstrates."*

---

## Outputs

| File | Description |
|------|-------------|
| `00_overview_composite.png` | 4-panel summary: DEM · Slope · Stress · Pf |
| `01_dem_hillshade.png` | DEM overlaid on hillshade |
| `02_slope_map.png` | Slope angle + critical zones >35° |
| `03_stress_field.png` | σv and τmax proxy maps |
| `04_failure_probability.png` | Pf map + stress–stability scatter |
| `05_summary_statistics.png` | Statistical distributions of all products |

---

## Libraries

```
numpy       — numerical operations
scipy       — gradient computation, filtering (replaces richdem)
rasterio    — raster GeoTIFF I/O
matplotlib  — publication-quality figures
```

---

## Usage

```bash
# From project root
python3 src/terrain_analysis.py
```

Outputs appear in `outputs/`. To use a real DEM, replace
`data/salzburg_dem.tif` with any GeoTIFF at any resolution.

---

## Author

**Ashutosh Pratap Shastri**
Junior Research Fellow, Department of Mining Engineering, IIT (BHU) Varanasi
M.Tech Geotechnical Engineering, IIT Patna (2025)
GitHub: [github.com/Ashutosh-Pratap-Shastri](https://github.com/Ashutosh-Pratap-Shastri)
