# Weighted Tweet Redistribution – GIS Project

This repository contains Python code implementing a weighted redistribution algorithm to visualise likely tweet locations based on administrative boundaries and population density. Built for GEOG3:71551 – Understanding GIS, the project extends Huck et al. (2015) and includes both efficient algorithm design and high-quality map visualisation outputs.

---

## Features

- Modularised functions: redistribution, seed generation, visualisation
- Support for both single-figure and 12-panel comparative outputs
- Population raster input (e.g. `100m_pop_2019.tif`)
- Optional reproducibility using fixed random seed
- Detailed comments and sectioned layout following rule-of-thirds coding style

---

## How to Use

### Setup
Install dependencies:
```bash
pip install geopandas rasterio shapely matplotlib scikit-image matplotlib-scalebar
```

### Required Input Files
- `gm-districts.shp` — Greater Manchester L3 districts
- `level3-tweets-subset.shp` — Geocoded tweet locations (collapsed to Level 3)
- `100m_pop_2019.tif` — UK population raster (British National Grid projection)

---

### Running the Code

#### Single Map
To produce a single standalone map with a scale bar, legend, and colour-coded likelihood surface:
```python
s = 0.1
sample_number = 25
output_surface = weighted_redistribution(admin_areas, point_data, raster, weight_data)
visualise_results(admin_areas, point_data, raster, output_surface)
```
This output is best suited for visualising the effects of a specific parameter configuration.

#### Multi-Panel Comparison Grid
To produce a 12-panel subplot comparing multiple `s` and `w` combinations:
```python
s_values = [0.05, 0.1, 0.2, 0.5]
w_values = [10, 25, 50]
visualise_multiple_distributions(admin_areas, point_data, raster, s_values, w_values)
```
Each panel will include a unique combination of `s` and `w`, with consistent styling and a shared legend.

---

## Algorithm Highlights

- Spatial indexing used to filter tweets efficiently by district
- Weighting surface values retrieved using raster row-column mapping
- Euclidean decay applied around highest-weight seed point
- Seed generation constrained within administrative geometry bounds
- Transparent colormap and perceptually uniform colour scale for clarity

---

## Analytical Considerations

- **Modifiable Areal Unit Problem (MAUP)**: Outputs are sensitive to how spatial units (districts) are defined and aggregated.
- **Parameter Subjectivity**: No objective criteria for selecting `s` or `w`; values must be interpreted contextually.
- **Ground Truth Unknown**: No way to verify redistributed tweet accuracy — results are plausible approximations.
- **False Positives / Dataset Completeness**: The dataset includes only tweets collapsed to L3 level, meaning many tweets are omitted. This may lead to underrepresentation in some areas.

---

## Author
Max Larner (10906064)  
University of Manchester – BSc Geography  
GEOG3:71551 – Understanding GIS
