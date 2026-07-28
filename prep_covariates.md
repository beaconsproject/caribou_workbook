# Covariate Selection Workplan
### Northern Mountain Caribou RSF / iSSF — 38 individuals, 2020–2026

---

## Phase 0 — Framework decisions (do this before touching data)

These choices determine *which* covariates you need and *at what resolution/extent*, so nail them down first.

| Decision | Why it matters | Suggested approach |
|---|---|---|
| **Order of selection** | 2nd-order (home range within study area) vs 3rd-order (point within home range) RSFs need different extents for "available" points and different covariate resolutions | Northern mountain caribou studies typically run 3rd-order RSFs (or iSSFs, which are inherently 3rd/4th-order) at the individual/seasonal level, often with a 2nd-order model layered on for range-level inference |
| **Seasonal stratification** | Northern mountain caribou show strong altitudinal migration (low-elevation old-growth in winter → alpine in summer/rut, distinct calving) — pooling seasons will wash out real selection signals | Run a **Net Squared Displacement (NSD) or k-means movement-phase analysis** on your own GPS data first (e.g., `migrateR`, `moveHMM`, or manual NSD fitting) to empirically define season breaks per animal/herd rather than relying purely on calendar dates from the literature |
| **iSSF step duration / fix rate** | Your covariates need to be extractable at the temporal resolution of your steps; irregular fix schedules need regularizing | Check fix-rate consistency across the 38 animals and across 2020–2026 (collar models/programming may have changed) — resample to a common interval with `amt::track_resample()` before building steps |
| **Spatial extent for "available"** | Determines the buffer around used points / step endpoints, and the clipping extent for every raster you download | For iSSF this is defined by the step-length distribution (empirical); for RSF, use MCP/KDE home ranges + buffer, or (for 2nd order) the herd range polygon |
| **Individual vs herd-level models** | With 38 animals over 6 years you likely span multiple herds/subpopulations with different disturbance regimes | Plan on mixed-effects RSF/iSSF (random slopes by animal or herd) from the start — this affects nothing about covariate *selection* but is worth deciding now since it affects how you structure the covariate extraction (long-format by animal-year)) |

---

## Phase 1 — Build an a priori, hypothesis-driven covariate list

Resist the urge to grab "everything available." Reviewers (and VIF) will punish a kitchen-sink approach. Build a table like this from the literature *before* you start downloading anything — it becomes your data-acquisition checklist.

Northern mountain caribou ecology centers on a well-established set of pressures: **predation risk avoidance** (spatial separation from wolves/moose/deer via elevation and cover), **forage** (terrestrial/arboreal lichen), **snow** (constrains both forage access and predator access), and **anthropogenic disturbance avoidance** (roads, cutblocks, seismic lines act as predator-facilitation corridors, not just direct avoidance). Good foundational literature to build hypotheses from: Johnson et al., DeCesare et al., Whittington et al., Dickie et al. (linear features & wolf movement), Apps & McLellan, and BC/Canada federal & provincial recovery strategy documents for the Northern Mountain population — these will also point you to which covariates were significant (or not) in nearby herds, which is useful prior information.

| Category | Example covariates | Ecological rationale | Season-dependence |
|---|---|---|---|
| **Topography** | Elevation, slope, aspect, terrain ruggedness index (TRI), terrain position index (TPI) | Elevation gradient drives the whole migratory strategy; ruggedness = predator-escape terrain | Elevation selection flips seasonally (low in winter, high in summer/calving) |
| **Land cover / vegetation** | Forest age/structure, canopy cover, open pine/spruce lichen woodland, alpine/subalpine parkland, wetland/bog | Direct habitat structure; old-growth lichen-bearing stands are the classic winter resource | Strongly season-specific |
| **Forage / productivity proxy** | NDVI/EVI (green-up timing & magnitude), lichen probability layers if available | Spring green-up and forage phenology; lichen abundance is rarely mapped directly except via provincial VRI attributes | Summer/calving mainly |
| **Snow** | Snow depth, snow water equivalent (SWE), snow-free date | Constrains both forage access and predator (wolf) movement efficiency | Winter/spring critical |
| **Fire history** | Time-since-fire, burn severity | Resets lichen recovery trajectory (decades-long); also affects moose/deer (apparent competition) habitat | Static-ish but must match GPS year |
| **Anthropogenic disturbance** | Distance to road, distance to cutblock edge, distance to seismic line/pipeline, cumulative human footprint index, well density | Predation-risk facilitation is the dominant mechanism in mountain/boreal caribou declines, not just direct disturbance | Mostly non-seasonal but avoidance strength often varies by season |
| **Predation-risk proxy** | Distance to low-elevation/early-seral habitat (moose/deer proxy), modeled wolf resource selection surfaces if available | Apparent competition is a central driver in this system | Often modeled as its own risk layer, not raw covariate |
| **Hydrology** | Distance to water, wetland density | Secondary influence on movement/foraging | Minor |
| **Protected status** | Inside/outside protected area or caribou habitat order | Useful for management-relevant models, not necessarily selection itself | Static |

**Deliverable for this phase:** a covariate table (category, variable, hypothesized direction, source, resolution needed, temporal resolution needed) that you can literally paste into your methods section later.

---

## Phase 2 — Data sources (R-first)

I checked current package status rather than relying on older documentation, because a couple of the "standard" MODIS packages in older tutorials are no longer maintained.

Layers that are national/circumpolar in coverage (terrain, MODIS/Sentinel, ERA5, global human footprint) apply identically to both jurisdictions — only the fine-scale, government-specific layers (vegetation/forest inventory, roads, disturbance, fire, hydrology, protected areas) need a separate BC vs. Yukon source. Yukon's open-data equivalent to the BC Data Catalogue is **GeoYukon** (`open.yukon.ca`), which exposes data two ways useful from R: a **CKAN API** (`open.yukon.ca/api/3/action/...`, readable with `httr`/`jsonlite`, or the generic `ckanr` package) and an **ArcGIS REST server** (`mapservices.gov.yk.ca/arcgis/rest/services/GeoYukon/...`), which is the more practical route for vector layers and works well with the `arcgislayers` package (or `esri2sf`-style REST queries) to pull data straight into `sf` objects — there's no `bcdata`-equivalent dedicated R package for Yukon, but the REST endpoints make this nearly as painless. One ecological note: Yukon herds see much less forestry/oil-and-gas disturbance than BC herds — the dominant anthropogenic layers there are **placer/quartz mining, the territorial highway network, and trails**, not cutblocks/seismic lines, so your disturbance covariate set may genuinely differ by jurisdiction rather than just the data source.

| Covariate | BC source | Yukon source | R acquisition path | Notes |
|---|---|---|---|---|
| Elevation / DEM | SRTM/CDEM | SRTM/CDEM (same, national coverage) | `elevatr::get_elev_raster()` | One source for both jurisdictions |
| Slope, aspect, TRI, TPI | — | — | `terra::terrain(dem, v = c("slope","aspect","TRI","TPI"))` | Derived, not downloaded |
| Land cover (national baseline) | Land Cover of Canada 2020 (30 m), or ESA WorldCover (10 m) | Same national/global products | `geodata::landcover()`; ESA WorldCover via direct URL + `terra::rast()`, or Earth Engine | Coarse for lichen-specific classes; use as fallback outside fine-scale inventory coverage |
| Land cover (fine-scale) | **BC Vegetation Resource Inventory (VRI)** — stand age, canopy cover, species composition | **Yukon Vegetation Inventory** (GeoYukon, 1:40,000 management-level; finer 1:10,000–1:5,000 coverage exists around some communities) — stand age/species/cover attributes | `bcdata::bcdc_query_geodata()` for BC; ArcGIS REST pull (`arcgislayers`) or CKAN download for Yukon | Yukon VI is management-level (coarser and patchier than VRI) — check coverage over your specific herd ranges before relying on it as the primary habitat-structure layer |
| Vegetation productivity / green-up (NDVI/EVI) | MODIS MOD13Q1/MYD13Q1, or Sentinel-2 | Same (national/global product) | **`rgee`** (Earth Engine via R, actively maintained) is the most reliable route; `luna` (rspatial) as an alternative for MODIS search/download. *Avoid `MODIStsp` — archived from CRAN Dec 2023.* | GEE route also gives per-pixel green-up timing (IRG) server-side — useful for iSSF forage-surfing covariates |
| Snow depth / SWE | ERA5-Land reanalysis (SNODAS has poor Canadian coverage) | Same — ERA5-Land | `ecmwfr` (free Copernicus CDS API key) or pull from the ERA5-Land collection via `rgee` | One source for both jurisdictions |
| Fire history | BC fire perimeters | **Yukon Fire History** (GeoYukon — landscape-level fire polygons back to 1946) + Fire Ignition Locations | `bcdata` for BC; CKAN API (`package_show?id=fire-history`) or ArcGIS REST for Yukon | Yukon's series goes back further (1946) than most BC products — useful if you want a long time-since-fire covariate, but check the pre-1997 minimum mapped fire size (200 ha) for bias toward large fires only |
| Roads | BC Digital Road Atlas | **Yukon Road Network** (GeoYukon) | `bcdata` for BC; ArcGIS REST/CKAN for Yukon | `osmdata` (OpenStreetMap) is a useful cross-check/supplement in both jurisdictions, especially for informal trails |
| Cutblocks / seismic lines / forestry disturbance | BC Consolidated Cutblocks, BC Oil & Gas Commission seismic lines | Generally minor in Yukon — check GeoYukon's forestry/land-use layers, but expect sparse coverage relative to BC | `bcdata` for BC; GeoYukon REST for Yukon if present | This is the layer most likely to need a genuinely different variable, not just a different source (see note above) |
| Mining disturbance | Present but secondary in BC | **Placer and quartz mining claims/workings** (GeoYukon — mineral tenure layers) | `bcdata` for BC mining layers; ArcGIS REST for Yukon mineral tenure | Worth adding as its own disturbance covariate for Yukon herds rather than folding into a generic "human footprint" layer, given how spatially concentrated placer mining is in Yukon watersheds |
| Human footprint index | Canadian Human Footprint (Watmough et al.) or Global Human Modification | Same national product | Direct raster download → `terra::rast()` | Same source for both; useful composite if you don't want to fit every disturbance type separately |
| Hydrology / wetlands | BC Freshwater Atlas | GeoYukon hydrology layers | `bcdata` for BC; ArcGIS REST for Yukon | |
| Protected areas / habitat protection | CPCAD, BC caribou habitat orders | CPCAD (national, covers Yukon too), + GeoYukon protected areas layers, Yukon caribou habitat/recovery planning boundaries if available | `bcdata` or direct download for BC-specific orders; CPCAD download or GeoYukon REST for Yukon | |
| Climate normals (optional) | ClimateNA, WorldClim | Same — both cover Yukon | ClimateNA is desktop software with a companion R wrapper, not a clean CRAN package; WorldClim via `geodata::worldclim_tile()` | Usually secondary to snow/NDVI for caribou work |

**R package cheat-sheet (currently maintained, worth installing up front):**
`terra`, `sf`, `elevatr`, `geodata`, `bcdata` (BC herds), `arcgislayers` or `httr`/`jsonlite` for GeoYukon's ArcGIS REST/CKAN endpoints (Yukon herds), `rgee` (+ Python/GEE setup), `ecmwfr`, `osmdata`, `amt` (movement/iSSF data structuring), `landscapemetrics` (patch-level metrics), `usdm` or `car` (VIF/collinearity screening).

---

## Phase 3 — Acquisition & preparation pipeline

1. **Define a master CRS and study-area polygon** (buffered herd ranges) — reproject everything to this once, at the start, so you're never guessing later.
2. **Download at the finest resolution you'll plausibly need**, then resample down — it's much easier to aggregate than to fabricate resolution you don't have.
3. **Match covariate year to relocation year.** This is the part people most often skip and it matters a lot given your 2020–2026 span: land cover, cutblocks, seismic lines, and fire perimeters all change over 6 years. Build a lookup so each animal-year of GPS data pulls the correct vintage of dynamic layers (VRI update cycle, annual cutblock/seismic layers, NDVI/snow by actual date).
4. **Build derived layers once, store them, don't recompute per-model:**
   - Distance-to-feature rasters (roads, cutblocks, seismic lines, water) via `terra::distance()`
   - Terrain derivatives via `terra::terrain()`
   - Time-since-fire via raster algebra on fire-year polygons rasterized per GPS-year
   - Consider **multi-scale versions** (e.g., % disturbed area in 250 m/500 m/1 km moving windows) since caribou avoidance of disturbance often operates at a scale broader than the raw pixel — this is a well-documented issue in this literature and worth testing a couple of buffer sizes early rather than assuming the native raster resolution is the right scale.
5. **Generate available points / random steps**
   - RSF: `sf`-based random point generation within the appropriate available-area polygon
   - iSSF: `amt::random_steps()` from your resampled tracks (this handles turn-angle/step-length distributions for you)
6. **Extract covariate values** to used points and available points/step-endpoints with `terra::extract()`, keeping the animal-year-season match established in step 3.
7. **Assemble one long-format analysis table**: `animal_id, herd, year, season, used(0/1), step_id (for iSSF), x, y, [covariates...]`.

---

## Phase 4 — Screening before modeling

- **Correlation/VIF check** on the full covariate set per season (`usdm::vifstep()` or `car::vif()` on a working GLM) — drop or combine before fitting the real models, not after convergence problems appear.
- **Univariate/exploratory screening** (simple GLMs, one covariate at a time) to catch coding errors and wildly implausible relationships before building multivariable models.
- **Scale-of-effect testing** for the disturbance and forage variables if you built multi-scale versions in Phase 3 — pick the buffer size with the best univariate fit per variable, then carry that version forward.
- **Standardize** continuous covariates (z-scores) before final model fitting, especially if you'll be comparing coefficient magnitudes across seasons/herds.
- **Sanity-map** a couple of finished covariate rasters against a known part of the landscape (e.g., check that your cutblock distance layer actually shows zero on a cutblock) before extracting to thousands of points.

---

## Rough sequencing / time budget

| Step | Approx. effort |
|---|---|
| Phase 0 (framework decisions + season definition from your own GPS data) | 1–2 weeks |
| Phase 1 (literature review + covariate table) | 1–2 weeks |
| Phase 2 (source identification, account/API setup for GEE, BC Data Catalogue, etc.) | 3–5 days |
| Phase 3 (download, reprocess, derive, extract) | 3–5 weeks — this is the long pole, mostly due to fire/cutblock/seismic year-matching and multi-scale derived layers |
| Phase 4 (screening) | 1 week |

---

## Handling the BC/Yukon split in practice
Since terrain, MODIS/Sentinel, ERA5, and national human-footprint layers are identical across both jurisdictions, the cleanest approach is usually a single processing pipeline with a `jurisdiction` flag per animal-year that switches only the fine-scale vegetation/disturbance/fire acquisition calls (`bcdata` vs. GeoYukon REST) — keep the derived-layer and extraction code (Phase 3, steps 4–7) identical for both, and build a lookup table mapping each herd/animal to its jurisdiction and study-area polygon at the very start of Phase 2 so you don't end up maintaining two parallel scripts.