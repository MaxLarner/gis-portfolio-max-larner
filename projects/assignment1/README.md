# Assignment 1 – Shortest International Border Finder

This project implements a Python-based GIS algorithm to identify the **shortest international land border** between two countries using a shapefile of global boundaries. The final output is a dynamic map highlighting the shortest border, based on ellipsoidal geodesic distance measurements.

---

## Features
- Reads and validates Natural Earth Admin 0 country shapefile
- Uses ellipsoidal measurements (not planar) for geodetic accuracy
- Implements a spatial index for efficient pairwise comparison
- Skips invalid border geometries (e.g. collapsed to points)
- Uses geodesic `Geod.inv()` distance calculation from `pyproj`
- Outputs a map with scale bar, north arrow, and a contextual inset map

---

## Requirements
Install dependencies via pip:
```bash
pip install geopandas shapely pyproj matplotlib matplotlib-scalebar
```

---

## Input Data
**Natural Earth Admin 0 Countries:**  
All data can be found from Natural Earth in /data/

---

## How to Run
From the project root:
```bash
python assignment1.py
```
This will:
- Parse the dataset
- Identify the shortest border
- Generate a visual map with all key cartographic elements

---

## Output
The final map is saved in the `/outputs/` directory and includes:
- Country A and Country B in contrasting colours
- The shortest border line highlighted
- A bottom-left inset showing regional context

---

## Skills Demonstrated
- Efficient spatial indexing with R-trees
- Handling of MultiLineString and LineString geometries
- Error handling with try-except for robust file loading
- Clear modular code and `if __name__ == '__main__'` structure
- Custom mapping with legends, annotation, and inset axes in `matplotlib`

---

## Author
Max Larner (10906064)  
BSc Geography – University of Manchester  
Module: GEOG3:71551 – Understanding GIS