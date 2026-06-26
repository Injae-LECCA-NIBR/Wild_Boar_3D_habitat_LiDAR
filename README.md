# Wild_Boar_3D_habitat_LiDAR

Data and code for: **"Drone-based LiDAR reveals behavior-specific 3D habitat attributes of wild boar (*Sus scrofa*) in urban-edge forests"** (manuscript ECOINF-D-26-01118, *Ecological Informatics*).

## Repository structure

```
wild-boar-3d-habitat-lidar/
├── data/
│   ├── Ethogram.csv                      # consolidated 10-day ethogram (behavior input)
│   ├── wildboar_environment_dataset.csv  # per-point habitat variables (analysis input)
│   ├── variable_description.csv          # data dictionary for the dataset above
│   └── DEM_study_site.tif                # 1-m LiDAR DEM (EPSG:32652)
├── scripts/
│   ├── ray_terrain_intersection.py       # flight log + ethogram -> boar coordinates
│   └── statistical_analysis.R            # track-level odds-ratio analysis
└── README.md
```

## Workflow

1. **Behavior classification** — UAV thermal/RGB observations are classified in
   the ethogram (`data/Ethogram.csv`) into traveling, daytime resting, and
   nighttime resting.
2. **Coordinate extraction** — `scripts/ray_terrain_intersection.py` projects a
   ray from each UAV pose (position, altitude, gimbal azimuth/pitch) and
   intersects it with the DEM to recover the ground coordinate of each observed
   boar.
3. **Habitat variables** — habitat rasters are overlaid on the coordinates in
   QGIS to produce `data/wildboar_environment_dataset.csv` (one row per point,
   one column per variable).
4. **Statistical analysis** — `scripts/statistical_analysis.R` runs the
   track-level univariate odds-ratio analysis with Benjamini-Hochberg
   correction.

## Reproducibility scope

The statistical analysis is fully reproducible from this repository:
`statistical_analysis.R` runs entirely on `wildboar_environment_dataset.csv` and
reproduces the reported results.

Coordinate extraction is not reproducible from the public files alone.
`ray_terrain_intersection.py` is provided for transparency, but the UAV flight
logs it reads and the individual boar coordinates it produces are restricted by
the data owner (see *Data availability*) and are available on request.

## Dataset notes

`wildboar_environment_dataset.csv` is the analysis-ready table: one row per
observation point, with habitat variables sampled at each point. Points without
a recovered coordinate were removed, and foraging was recoded to traveling
(anthropogenic-food-source foraging was identified by location in GIS and
removed beforehand). Individual boar coordinates are not included (see
*Data availability*).

`variable_description.csv` lists the 14 candidate analysis variables, each mapped
to its raw dataset column via `source_column`: `aspect` provides `northness` and
`eastness`, and `deck`/`trail` are stored as raw distances and
log(x + 1)-transformed in analysis. Twelve variables were retained after
multicollinearity screening (`gap_fraction` and `slope` excluded). Variable
derivations follow the methods cited in the manuscript.

## Coordinate reference systems

The DEM CRS is read dynamically from the raster (`dem.crs`); the DEM provided
here is **EPSG:32652** (WGS 84 / UTM zone 52N). Habitat rasters sampled in QGIS
are in **EPSG:5186** (Korea 2000 / Central Belt 2010). Extracted boar
coordinates are **WGS 84 geographic, EPSG:4326**.

## Running the extraction

```bash
pip install rasterio pyproj pandas numpy openpyxl
python scripts/ray_terrain_intersection.py
```

Run from the repository root. Flight logs are required and available on request
(see *Data availability*).

## Data availability

The consolidated ethogram, the analysis-ready habitat-variable dataset
(excluding individual coordinates), the digital elevation model, and all code
are openly available here and archived at Zenodo (DOI: 10.5281/zenodo.XXXXXXX).

Individual boar coordinates, UAV flight logs, and raw video are not released as
fully open data because, although wild boar are not an endangered species, the
data owner (Ministry of Climate, Energy and Environment, Republic of Korea)
restricts their full open-source release. These materials are available from the
corresponding author on reasonable request.

## License

Code: MIT.
Data: © Ministry of Climate, Energy and Environment. The Government of the Republic of Korea; available under the terms set by the data owner.
