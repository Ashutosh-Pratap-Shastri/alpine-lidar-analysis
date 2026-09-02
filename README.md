# 🏔️ Alpine Geohazard Analysis — Hochkönig / Salzburg Alps

\[!\[Streamlit App](https://static.streamlit.io/badges/streamlit\_badge\_black\_white.svg)](https://alpine-geohazard-analysis.stream
!\[Python](https://img.shields.io/badge/Python-3.12-blue)
!\[License](https://img.shields.io/badge/License-MIT-green)

**DEM-based slope stability and topographic stress assessment for high-alpine environments.**  
Preparatory portfolio work for the **CRAG PhD position** (*Critical Infrastructure Risk from Alpine Geohazards*) at the Paris Lodron University of Salzburg.

\---

## 🌐 Live Web App

> \*\*Anyone can run this analysis in their browser — no installation needed.\*\*  
> Upload your own GeoTIFF DEM or use the built-in Salzburg Alps terrain.  
> Adjust rock mass parameters, run Monte Carlo, download figures.

👉 👉 \*\*\[Launch the app on Streamlit Cloud](https://alpine-geohazard-analysis.streamlit.app)\*\*

\---

## 📊 Output Figures

### Overview (4-panel composite)

!\[Overview](outputs/00\_overview\_composite.png)

### DEM + Hillshade

!\[DEM](outputs/01\_dem\_hillshade.png)

### Slope Map — Critical Zones >35°

!\[Slope](outputs/02\_slope\_map.png)

### Topographic Stress Field

!\[Stress](outputs/03\_stress\_field.png)

### Failure Probability (Monte Carlo)

!\[Failure probability](outputs/04\_failure\_probability.png)

### Summary Statistics

!\[Stats](outputs/05\_summary\_statistics.png)

\---

## What This Does

|Step|Module|Output|
|-|-|-|
|1|`load\_dem`|Elevation array + georeferencing metadata|
|2|`compute\_terrain\_derivatives`|Slope, aspect, plan curvature (Horn 1981)|
|3|`compute\_stress\_proxy`|Gravitational dead-load stress proxy (σv, σh, τmax)|
|4|`compute\_failure\_probability`|Monte Carlo Pf map (n = 10,000, infinite-slope model)|
|5|Figures|6 publication-quality PNG outputs|

\---

## Key Results (Salzburg Alps synthetic DEM)

```
Elevation range        :  580 – 2,941 m
Mean slope             :  24.2°
Critical slopes (>35°) :  20.1% of area
Max failure prob Pf    :  0.795
High-risk (Pf > 0.10)  :  10.6% of area
Vertical stress σv     :  15 – 76 MPa
Max shear stress τmax  :  up to 25.5 MPa
```

\---

## Study Area

**Hochkönig / Salzburg Alps**  
Bounding box: 13.05–13.20°E, 47.38–47.50°N · Resolution: \~18–22 m/pixel

> \*\*Data note:\*\* The bundled DEM is a high-fidelity synthetic terrain built from
> published geomorphological parameters for the Northern Calcareous Alps
> (elevation 580–2941 m; glacially carved valleys; limestone massif character).
> Any real SRTM-30m or Copernicus GLO-30 tile for this bounding box replaces it
> with zero code changes.

\---

## Methodological Notes

### Terrain derivatives

Horn (1981) 3×3 finite-difference gradient estimator — same kernel as ArcGIS `Slope`, QGIS, `gdaldem`.

### Stress proxy

```
σv   = ρ g H / 10⁶          lithostatic vertical stress (MPa)
K0   = ν / (1 − ν)           at-rest coefficient
σh   = K0 σv (1 + sin α)     slope-amplified horizontal stress
τmax = |σv − σh| / 2         maximum shear stress
```

> ⚠️ This is a simplified 2-D dead-load proxy — \*\*not\*\* the FCM 3-D volumetric
> solution (Haunsperger \& Robl 2025, 2026). The FCM operates on full volumetric
> meshes across entire massifs on HPC clusters.

### Failure probability

```
c   \~ N(50, 15) kPa     cohesion (uncertain)
φ   \~ N(35,  5) °       friction angle (uncertain)
FS  = (c + σn tan φ) / τ
Pf  = P(FS < 1.0),  n = 10,000 Monte Carlo realisations
```

Extends the probabilistic FEM methodology from the author's M.Tech thesis
(Rosenblueth PEM + Monte Carlo in RS2) to DEM-scale alpine terrain.

\---

## Connection to CRAG

The FCM-based 3-D stress modelling in CRAG operates on the same conceptual chain —
**topography → stress field → failure threshold** — but resolved at full massif
scale with volumetric finite cell meshes on HPC clusters, coupled with real
multiscale monitoring data (InSAR, UAV, TLS, fissurometers).

> \*None of the six foundational CRAG papers computes a probabilistic Pf output.
> Haunsperger \& Robl produce deterministic stress tensors; Schneider-Muntau
> produces deterministic factors of safety. Propagating geological uncertainty
> through the model to obtain Pf is the step this pipeline demonstrates.\*

\---

## Usage

### Run locally

```bash
git clone https://github.com/Ashutosh-Pratap-Shastri/alpine-lidar-analysis
cd alpine-lidar-analysis
pip install -r requirements.txt
python3 src/terrain\_analysis.py      # batch pipeline → outputs/
streamlit run app.py                 # interactive web app
```

### Use your own DEM

Replace `data/salzburg\_dem.tif` with any GeoTIFF — or upload it directly
in the web app.

\---

## Libraries

|Library|Purpose|
|-|-|
|`numpy`|Numerical operations|
|`scipy`|Gradient computation, filtering|
|`rasterio`|GeoTIFF I/O|
|`matplotlib`|Publication-quality figures|
|`streamlit`|Interactive web app|

\---

## Author

**Ashutosh Pratap Shastri**  
Junior Research Fellow, Department of Mining Engineering, IIT (BHU) Varanasi  
M.Tech Geotechnical Engineering, IIT Patna (2025)  
[github.com/Ashutosh-Pratap-Shastri](https://github.com/Ashutosh-Pratap-Shastri)

