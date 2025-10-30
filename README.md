[![2015 Stuttgart light rail network maps generated from GTFS data, with optimal line orderings, geographically correct (left), octilinear (middle), and orthoradial (right).](examples/render/stuttgart-example-small.png?raw=true)](examples/render/stuttgart-example.png?raw=true)
<sub>*2015 Stuttgart light rail network maps generated from GTFS data, with optimal line orderings, geographically correct (left), octilinear (middle), and orthoradial (right).* </sub>

# 🚊 Loom

[![Build](https://github.com/ad-freiburg/loom/actions/workflows/build.yml/badge.svg)](https://github.com/ad-freiburg/loom/actions/workflows/build.yml)

Create beautiful, geographically accurate or schematic transit maps from GTFS feeds and custom line graphs. Loom orchestrates data preparation, optimization, and rendering so you can focus on designing compelling transit experiences.

---

## 📋 Table of Contents

1. [Highlights](#-highlights)
2. [Research Background](#-research-background)
3. [Live Demos](#-live-demos)
4. [Requirements](#-requirements)
5. [Installing FreeType](#-installing-freetype)
6. [Building Loom](#-building-loom)
7. [Toolchain Overview](#-toolchain-overview)
8. [Quickstart](#-quickstart)
9. [Tips & Tricks](#-tips--tricks)
10. [Key `loom` Parameters](#-key-loom-parameters)
11. [Configuration](#-configuration)
12. [Tool Capabilities](#-tool-capabilities)
13. [Command-Line Reference](#-command-line-reference)
14. [Line Graph Extraction from GTFS](#-line-graph-extraction-from-gtfs)
15. [Landmarks & Location Markers](#-landmarks--location-markers)
16. [Docker Usage](#-docker-usage)

---

## ✨ Highlights

- **End-to-end pipeline** for transforming raw GTFS feeds into printable or digital transit maps.
- **Flexible layouts**: generate geographically faithful, octilinear, or orthoradial schematics with a single workflow.
- **Powerful optimization** engine with ILP, hybrid, stochastic, and greedy solvers tailored to complex transit networks.
- **Rich rendering controls** to fine-tune labels, padding, and exported layers for production-ready output.

## 📚 Research Background

Loom is based on the following publications:

- [Bast H., Brosi P., Storandt S., Efficient Generation of Geographically Accurate Transit Maps, SIGSPATIAL 2018](http://ad-publications.informatik.uni-freiburg.de/SIGSPATIAL_transitmaps_2018.pdf)
- [Bast H., Brosi P., Storandt S., Efficient Generation of Geographically Accurate Transit Maps (extended version), ACM TSAS, Vol. 5, No. 4, Article 25, 2019](http://ad-publications.informatik.uni-freiburg.de/ACM_efficient%20Generation%20of%20%20Geographically%20Accurate%20Transit%20Maps_extended%20version.pdf)
- [Bast H., Brosi P., Storandt S., Metro Maps on Octilinear Grid Graphs, EuroVis 2020](http://ad-publications.informatik.uni-freiburg.de/EuroVis%20octi-maps.pdf)
- [Bast H., Brosi P., Storandt S., Metro Maps on Flexible Base Grids, SSTD 2021](http://ad-publications.informatik.uni-freiburg.de/SSTD_Metro%20Maps%20on%20Flexible%20Base%20Grids.pdf)

A pipeline similar to Loom was also described by Anton Dubrau in this [Transit App blog post](https://blog.transitapp.com/how-we-built-the-worlds-prettiest-auto-generated-transit-maps-12d0c6fa502f).

## 🌐 Live Demos

- [Geographically accurate renderer](https://loom.cs.uni-freiburg.de/)
- [Global network explorer](https://loom.cs.uni-freiburg.de/global)
- [Octilinear schematic explorer](https://octi.cs.uni-freiburg.de)

## 🧰 Requirements

- `cmake`
- `gcc >= 5.0` (or `clang >= 3.9`)
- `libicu` development files (e.g., `libicu-dev`)
- Optional: `libglpk-dev`, `coinor-libcbc-dev`, `gurobi`, `libzip-dev`, `libprotobuf-dev`
- Optional: `freetype2` (libfreetype) for accurate label rendering

## 🧱 Installing FreeType

The map renderer uses the FreeType font engine to measure label text and the TT Norms Pro font bundled under `data/fonts` for rendering. Install the development package for your platform:

### Debian/Ubuntu

```bash
sudo apt-get install libfreetype6-dev
```

### macOS (Homebrew)

```bash
brew install freetype
```

### Windows (vcpkg)

```powershell
vcpkg install freetype
```

## 🏗️ Building Loom

Fetch this repository (including submodules):

```bash
git clone --recurse-submodules https://github.com/ad-freiburg/loom.git
```

If you cloned without `--recurse-submodules` or need to refresh dependencies later:

```bash
cd loom
git submodule update --init --recursive
```

Install the ICU development package before configuring the project (example for Debian/Ubuntu):

```bash
sudo apt-get install libicu-dev
```

Configure and build:

```bash
cd loom
mkdir build && cd build
cmake ..
make -j
```

Optionally install the binaries:

```bash
make install
```

> 💡 You can run the tools straight from `./build` without installing them system-wide.

## 🧰 Toolchain Overview

| Tool | Purpose |
| ---- | ------- |
| `gtfs2graph` | Convert GTFS feeds into GeoJSON line graphs. |
| `topo` | Clean and untangle graphs, remove overlapping segments, infer turn restrictions. |
| `loom` | Compute optimal line orderings using ILP, hybrid, or heuristic solvers. |
| `octi` | Produce schematic layouts on octilinear or orthoradial base grids. |
| `transitmap` | Render the final line graph to SVG or MVT vector tiles. |

## ▶️ Quickstart

The `examples/` folder contains several ready-made line graphs. To reproduce the Stuttgart map shown above:

```bash
cat examples/stuttgart.json | loom | transitmap > stuttgart.svg
```

Add labels:

```bash
cat examples/stuttgart.json | loom | transitmap -l > stuttgart-label.svg
```

Annotate only interchange segments (keep labels enabled but omit those on single-route edges):

```bash
cat examples/stuttgart.json | loom | transitmap -l --single-route-labels=false \
  > stuttgart-label-interchanges.svg
```

Generate an octilinear schematic:

```bash
cat examples/stuttgart.json | loom | octi | transitmap -l > stuttgart-octilin.svg
```

Switch `octi` to an orthoradial base graph:

```bash
cat examples/stuttgart.json | loom | octi -b orthoradial | transitmap -l > stuttgart-orthorad.svg
```

Adjust padding per side (unspecified sides default to `0`):

```bash
cat examples/stuttgart.json | loom | transitmap --padding-top 50 --padding-left 20 > stuttgart-pad.svg
```

## 💡 Tips & Tricks

- Run any tool with `-h` to inspect its full help text.
- Export individual layers by toggling renderer switches such as `--labels`, `--no-render-stations`, and `--no-render-node-connections`.
- Hide strokes entirely by setting `--line-width 0 --outline-width 0` before exporting a layer.

## 🎯 Key `loom` Parameters

Fine-tune preprocessing, optimization, and scoring with the following flag groups:

- **Preprocessing** – Disable untangling or pruning with `--no-untangle` and `--no-prune` if you must preserve the input topology.
- **Input format** – Parse DOT graphs by passing `-D/--from-dot` instead of GeoJSON.
- **Optimization strategy** – Select solvers via `-m/--optim-method`. Options include ILP variants (`ilp`, `ilp-naive`), hybrid (`comb`, `comb-no-ilp`), exhaustive search (`exhaust`), local search (`hillc`, `hillc-random`), simulated annealing (`anneal`, `anneal-random`), greedy strategies, and a `null` mode that keeps existing orderings. Combine with `--optim-runs <n>` to repeat stochastic solvers and keep the best result.
- **Penalty weights** – Adjust crossing and separation costs with `--same-seg-cross-pen`, `--diff-seg-cross-pen`, `--in-stat-cross-pen-same-seg`, `--in-stat-cross-pen-diff-seg`, `--sep-pen`, and `--in-stat-sep-pen`.
- **ILP solver configuration** – Configure solver, threads, and runtime via `--ilp-solver`, `--ilp-num-threads`, and `--ilp-time-limit`.
- **Diagnostics** – Enable textual or GeoJSON debug output using `--output-stats`, `--write-stats`, `--dbg-output-path`, and `--output-optgraph`.

> Start with the default `comb-no-ilp` strategy. Increase penalties for stubborn crossings or separations, and keep preprocessing enabled unless you need to inspect the raw input graph.

## 🛠️ Configuration

All tools consult optional `.loom.ini` files for default settings. Configuration is resolved in this order:

1. `$HOME/.loom.ini`
2. `.loom.ini` next to the executable
3. Command-line flags
4. `--config=<file>` overrides the above when supplied explicitly

Refer to [loom.ini](loom.ini) for available keys and defaults. Adjust the `log-level` key (or `--log-level`) from `0` (errors only) to `4` (verbose debug); the default is `2`.

### Terminus Route Labels

- `station-label` (default) – Anchor the route label stack to the positioned station name.
- `stop-footprint` – Anchor labels to the station polygon footprint.
- `node` – Use the raw node position when no station label or footprint exists.

### Station Label Candidates

- `station-label-angle-steps` – Number of evenly spaced orientations around each station. Must be positive and divisible by four.
- `station-label-angle-step-deg` – Angular distance between successive samples. Tune both settings to balance placement flexibility and runtime.

### Text Sizing

- Enable `textsize-constant` to interpret label sizes and spacing as typographic point sizes, converted to device-independent CSS pixels.
- Use `--textsize-range-pt <min,max>` to keep text resolution-dependent but clamp station and terminus label sizes to the provided range.

## 🧭 Tool Capabilities

- `gtfs2graph` – Generate GeoJSON line graphs from GTFS feeds.
- `topo` – Remove overlapping edges and infer turn restrictions.
- `loom` – Compute optimal line orderings for a line graph.
- `octi` – Convert a line graph into schematic layouts on base grids.
- `transitmap` – Render a line graph as SVG or MVT vector tiles.

## 🧾 Command-Line Reference

### `gtfs2graph`

- `-m`, `--mots <modes>` – MOTs to calculate shapes for; comma-separated list of mode names or GTFS codes (default `all`).
- `-p`, `--prune-threshold <0..1>` – Threshold for pruning seldomly occurring lines (default `0`).
- `-h`, `--help` – Show help message.
- `-v`, `--version` – Print version.

### `topo`

`topo` consumes the raw line graph created by `gtfs2graph` (or any equivalent tool) and cleans it up before layout. Key parameters include:

- `-d`, `--max-aggr-dist <meters>` – Maximum distance for merging parallel or overlapping segments (default `50`). Affects both the pre-pass and `MapConstructor::collapseShrdSegs` loop.
- `--sample-dist <length>` – Sampling step (in pseudometers) when constructing support points during aggregation (default `5`). Lower values follow the geometry more closely; higher values simplify it.
- `--smooth <factor>` – Optional smoothing factor applied after aggregation to soften sharp angles.
- `--max-comp-dist <meters>` – Maximum gap when grouping nodes into distance-based connected components (default `10000`). Controls how networks are partitioned.
- `--write-components`, `--write-components-path <dir>` – Annotate edges with component IDs and optionally export each component as its own GeoJSON file.
- `--infer-restr-max-dist <meters>` – Search radius for candidate turns during restriction inference (defaults to `--max-aggr-dist`).
- `--max-length-dev <meters>` – Maximum detour length accepted when validating inferred turn restrictions.
- `--turn-restr-full-turn-angle <angle>` / `--turn-restr-full-turn-pen <penalty>` – Penalize sharp “full turns.”
- `--no-infer-restrs` – Skip restriction inference entirely.
- `--random-colors` – Populate missing line colors before processing.
- `--write-stats`, `--aggr-stats` – Emit run statistics, optionally merging them with prior metadata when aggregation is requested.
- `-h`, `--help`; `-v`, `--version` – Standard CLI helpers.

### `loom`

- `--no-untangle` – Skip untangling rules.
- `--no-prune` – Skip pruning rules.
- `-m`, `--optim-method <method>` – Select optimization method (default `comb-no-ilp`).
- `--same-seg-cross-pen <weight>` – Penalty for same-segment crossings (default `4`).
- `--diff-seg-cross-pen <weight>` – Penalty for different-segment crossings (default `1`).
- `--in-stat-cross-pen-same-seg <weight>` – Penalty for same-segment crossings at stations (default `12`).
- `--in-stat-cross-pen-diff-seg <weight>` – Penalty for different-segment crossings at stations (default `3`).
- `--sep-pen <weight>` – Penalty for separations (default `3`).
- `--in-stat-sep-pen <weight>` – Penalty for separations at stations (default `9`).
- `--ilp-solver <solver>` – Preferred ILP solver (`glpk`, `cbc`, or `gurobi`; default `gurobi`).
- `--ilp-num-threads <n>` – Number of threads for the ILP solver (`0` for solver default).
- `--ilp-time-limit <sec>` – ILP solve time limit (`-1` for infinite runtime).
- `--output-stats`, `--write-stats` – Print or write statistics.
- `--dbg-output-path <path>` – Directory for debug output.
- `--output-optgraph` – Write the optimization graph to the debug path.
- `-h`, `--help`; `-v`, `--version` – Standard CLI helpers.

### `octi`

`octi` turns the topological line graph produced by `loom` into a schematic map that follows a configurable base grid. The flags below are grouped to mirror the CLI help.

**Optimization & input control**

- `-m`, `--optim-mode <heur|ilp>` – Choose between the fast heuristic placer (`heur`, default) and the exact ILP optimizer (`ilp`).
- `--obstacles <file>` – Provide a GeoJSON file with polygons that the layout must avoid when routing edges.
- `--edge-order <method>` – Pick how edges are ordered before the search (`num-lines`, `length`, `adj-nd-deg`, etc.; `all` tries several strategies).
- `--loc-search-max-iters <n>` – Limit how many refinement iterations the local improvement stage performs (default `100`).
- `--geo-pen <weight>` – Penalize deviation from the original line geometry so edges stay closer to their input shape (default `0`).

**Grid construction & base graph**

- `-g`, `--grid-size <len or %>` – Set the grid resolution as an absolute length or as a percentage of the average adjacent-station spacing (default `100%`).
- `-b`, `--base-graph <type>` – Select the base grid (`ortholinear`, `octilinear`, `orthoradial`, `quadtree`, or `octihanan`; default `octilinear`).
- `--hanan-iters <n>` – Number of refinement iterations when using the `octihanan` grid.
- `--max-grid-dist <n>` – Cap how many grid steps away from the original station position candidates may be generated (default `3`).

**Error handling**

- `--retry-on-error` – Retry the placement up to 30 times with an 85% grid size when the solver fails.
- `--skip-on-error` – Skip the current graph instead of aborting the pipeline when placement fails.

**ILP solver configuration**

- `--ilp-solver <solver>` – Choose the ILP backend (`glpk`, `cbc`, or `gurobi`; default `gurobi`).
- `--ilp-num-threads <n>` – Limit ILP solver threads (`0` uses the solver default).
- `--ilp-time-limit <sec>` – Cap solver runtime (`-1` removes the limit).
- `--ilp-cache-dir <dir>` – Directory used to cache subproblem solutions.
- `--ilp-cache-threshold <val>` – Minimum improvement required before caching an ILP result.

**Penalties & layout preferences**

- `--density-pen <weight>` – Discourage crowding by penalizing densely packed stations and edges (default `10`).
- `--vert-pen <weight>`, `--hori-pen <weight>`, `--diag-pen <weight>` – Weigh the use of vertical, horizontal, and diagonal grid directions.
- `--pen-180 <w>`, `--pen-135 <w>`, `--pen-90 <w>`, `--pen-45 <w>` – Add bend penalties for the corresponding angle changes.
- `--nd-move-pen <weight>` – Charge a cost for moving nodes away from their input location.

- `-h`, `--help`; `-v`, `--version` – Standard CLI helpers.

### `transitmap`

- `--line-width <px>` – Width of a single transit line (default `20`).
- `--line-spacing <px>` – Spacing between transit lines (default `10`).
- `--outline-width <px>` – Width of line outlines (default `1`).
- `--log-level <0..4>` – Logging verbosity, `0`=errors to `4`=very verbose (default `2`).
- `--render-dir-markers` – Render line direction markers (tails are always enabled when space allows).
- `--dir-marker-spacing <n>` – Edges between forced direction markers (default `1`).
- `--tail-ignore-sharp-angle` – Ignore the sharp-angle check when rendering marker tails (default off).
- `--bi-dir-marker` – Render markers for bidirectional edges (default off).
- `--crowded-line-thresh <n>` – Lines on edge to trigger direction marker (default `3`).
- `--sharp-turn-angle <rad>` – Turn angle in radians (0–π) to trigger direction markers (default `0.785398`). Values >π are treated as degrees.
- `-l`, `--labels` – Render labels.
- `-r`, `--route-labels` – Render route names at line termini.
- `--line-label-textsize <size>` – Text size for line labels (default `40`).
- `--line-label-bend-angle <rad>` – Max bend angle for line label candidates (default `0.349066`).
- `--line-label-length-ratio <ratio>` – Max length/straight distance ratio for line label candidates (default `1.1`).
- `--station-label-textsize <size>` – Text size for station labels (default `60`).
- `--textsize-constant` – Treat station labels, terminus route labels, and their spacing as typographic point values regardless of `--resolution`.
- `--textsize-range-pt <min,max>` – Clamp station and terminus label point sizes to a range while still scaling with `--resolution`.
- `--me-label-textsize <size>` – Text size for the "YOU ARE HERE" label (default `80`).
- `--font-svg-max <size>` – Max font size for station labels in SVG (`-1` for no limit; default `11`). Interpreted as typographic points when text-size options are active, otherwise SVG pixels.
- `--station-line-overlap-penalty <weight>` – Penalty multiplier for station-line overlaps (default `15`).
- `--station-line-overlap-per-line` – Count distinct transit lines when scoring station-line overlaps (default disabled).
- `--station-label-far-crowd-radius <px>` – Radius from the far end of a station label used to detect nearby features (default `0`, disables).
- `--station-label-far-crowd-penalty <weight>` – Penalty when the far label end crowds nearby features (default `25`).
- `--side-penalty-weight <weight>` – Weight for station label side preference penalties (default `2.5`).
- `--same-side-penalty <penalty>` – Penalty for station labels on opposite sides (default `100`).
- `--reposition-label <n>` – Perform additional passes after initial placement to relieve label crowding (default `0`).
- `--cluster-pen-scale <scale>` – Scale factor for station crowding penalties (default `1`).
- `--outside-penalty <weight>` – Penalty (positive) or bonus (negative) for labels outside the map bounds (default `-5`).
- `--orientation-penalties <p0,...,p7>` – Comma-separated penalties for eight label orientations (default `0,3,6,4,1,5,6,2`).
- `--terminus-angle-penalty <penalty>` – Penalty for non-axis-aligned terminus station labels (default `3`).
- `--route-label-gap <size>` – Gap between route label boxes (default `10`). Uses typographic points when text-size options are enabled, otherwise map units.
- `--route-label-terminus-gap <size>` – Gap between the terminus station label and the first route label box (default `80`).
- `--terminus-label-anchor <anchor>` – Anchor geometry for terminus route labels (`station-label`, `stop-footprint`, or `node`; default `station-label`).
- `--compact-terminal-label` – Arrange terminus route labels in multiple columns instead of a single row (default off).
- `--compact-route-label` – Stack edge route labels in multiple rows to avoid truncation (default off).
- `--highlight-terminal` – Highlight terminus stations (default off).
- `--terminus-highlight-fill <color>` – Fill color when highlighting terminus stations (default `black`).
- `--terminus-highlight-stroke <color>` – Stroke color for highlighted terminus stations (default `#BAB6B6`).
- `--no-deg2-labels` – Suppress labels for degree‑2 stations.
- `-D`, `--from-dot` – Input graph is in DOT format.
- `--resolution <res>` – Output resolution (default `0.1`).
- `--padding <padding>` – Padding around the map (`-1` for auto) applied to all sides when no side-specific padding is provided.
- `--padding-top <padding>` – Padding at the top of the map (`-1` for auto).
- `--padding-right <padding>` – Padding at the right side (`-1` for auto).
- `--padding-bottom <padding>` – Padding at the bottom (`-1` for auto).
- `--padding-left <padding>` – Padding at the left side (`-1` for auto).
- `--smoothing <factor>` – Input line smoothing (default `1`).
- `--ratio <value>` – Output width/height ratio (`width = height * ratio`).
- `--tl-ratio <value>` – Top-left anchored width/height ratio (default `-1`).
- `--random-colors` – Fill missing colors with random colors.
- `--tight-stations` – Do not expand node fronts for stations.
- `--no-render-stations` – Do not render stations.
- `--no-render-node-connections` – Do not render inner node connections.
- `--render-node-fronts` – Render node fronts.
- `--bg-map <file.geojson>` – Render additional GeoJSON geometry behind the network. Coordinates are expected in WGS84 latitude/longitude. Polygon and MultiPolygon features are supported; feature properties may specify `stroke`, `stroke-width`, `fill`, and `opacity`.
- `--bg-map-webmerc` – Treat `--bg-map` coordinates as already in Web Mercator.
- `--bg-map-opacity <value>` – Opacity for `--bg-map` geometry (default `1`).
- `--extend-with-bgmap` – Expand the output bounding box to include `--bg-map` geometry.
- `--geo-lock [bool]` – Ensure the map covers at least a default geographic bounding box (defaults to `true`).
- `--geo-lock-bbox <south,west,north,east>` – Custom bounding box for `--geo-lock` (latitude/longitude).
- `--landmark <spec>` – Add a landmark (`word:text,lat,lon[,fontSize[,color[,opacity]]]` or `iconPath,lat,lon[,size]`).
- `--landmarks <file>` – Read landmarks from a file, one per line.
- `--force-landmarks` – Keep landmarks even when they overlap existing geometry (default `true`). Pass `--force-landmarks=false` to displace or skip overlapping icons.
- `--landmark-search-radius <radius>` – Search radius for shifting overlapping landmark icons (default `10`).
- `--displacement-iterations <n>` – Maximum iterations for landmark displacement (default `100`).
- `--displacement-cooling <factor>` – Cooling factor for landmark displacement steps (default `0.9`).
- `--landmarks-webmerc` – Treat landmark and `--me` coordinates as already in Web Mercator.
- `--me <lat,lon>` – Mark the given coordinates with a red star (latitude/longitude by default).
- `--me-size <size>` – Star size (default `150`).
- `--me-label` – Add a "YOU ARE HERE" label.
- `--me-star[=<bool>]` – Force rendering of the "YOU ARE HERE" star when coordinates or a matching station are provided without other `--me` toggles (default `false`).
- `--me-station <name>` – Mark the current location by station label, preserving punctuation and casing for display while using a slugged identifier internally.
- `--me-with-bg[=<bool>]` – When `--me-station` is active, restyle the matched station label with a rounded badge; falls back to the standalone badge when the label cannot be restyled.
- `--me-bg-fill <color>` – Badge fill color for the `--me` background (default `#f5f5f5`).
- `--me-bg-stroke <color>` – Badge stroke color for the `--me` background (default `#d0d0d0`).
- `--me-label-color <color>` – Text color for the badge label (default `#3a3a3a`).
- `--me-station-fill <color>` – Fill color for the "me" marker (default `#f00`).
- `--me-station-border <color>` – Border color for the "me" marker (default none).
- `--print-stats` – Write statistics to stdout.
- `-h`, `--help`; `-v`, `--version` – Standard CLI helpers.


## 🗺️ Line Graph Extraction from GTFS

Use `gtfs2graph` to extract a line graph from a GTFS feed. For example, to build a Freiburg tram network graph:

```bash
gtfs2graph -m tram freiburg.zip > freiburg.json
```

The result will contain overlapping edges and stations. Clean it with `topo` before rendering:

```bash
gtfs2graph -m tram freiburg.zip | topo > freiburg.json
```

A full pipeline for creating an octilinear Freiburg map combines all tools:

```bash
gtfs2graph -m tram freiburg.zip | topo | loom | octi | transitmap -l > freiburg.svg
```

Use `--optim-method`, `--base-graph`, and renderer options to tailor the output to your network.

## 📍 Landmarks & Location Markers

Enhance maps with point-of-interest annotations using `transitmap` landmarks.

- `--landmark <spec>` – Add a single landmark inline. Use `word:text,lat,lon[,fontSize[,color[,opacity]]]` for text labels or `iconPath,lat,lon[,size]` for icons.
- `--landmarks <file>` – Load one landmark per line from a file. Relative icon paths are resolved against the landmark file.
- `--force-landmarks` – Keep landmarks even when they overlap existing geometry (default `true`). Pass `--force-landmarks=false` to allow the solver to displace or skip conflicting entries.
- `--landmark-search-radius <radius>` – Search radius used when nudging overlapping icon landmarks (default `10`).
- `--displacement-iterations <n>` / `--displacement-cooling <factor>` – Control the landmark displacement annealing process (defaults `100` / `0.9`).
- `--landmarks-webmerc` – Treat landmark coordinates as Web Mercator instead of WGS84.
- `--me`, `--me-station`, `--me-with-bg`, `--me-star`, and related `--me-*` flags – Highlight a specific location (for example "YOU ARE HERE") with a star and optional badge styling.

Example: render the sample landmarks provided in `examples/landmarks.txt`:

```bash
cat examples/stuttgart.json | loom | transitmap --landmarks examples/landmarks.txt > stuttgart-landmarks.svg
```

Add a single inline landmark from the command line:

```bash
cat examples/stuttgart.json | loom | transitmap --landmark word:CityHall,47.9210,106.9175,150,#ff0000 > stuttgart-cityhall.svg
```

The sample file uses locations around Ulaanbaatar, Mongolia; human-readable names are listed in `examples/ulaanbaatar_landmark_names.txt`.

## 🐳 Docker Usage

A Dockerfile is provided for reproducible builds.

### Build the image

```bash
docker build -t loom .
```

### Run a tool from the suite

```bash
docker run -i loom <TOOL>
```

Example: octilinearize the Freiburg demo:

```bash
cat examples/freiburg.json | sudo docker run -i loom octi
```

> **Note:** To use Gurobi inside the container, mount a directory containing a valid `gurobi.lic` to `/gurobi/`:

```bash
docker run -v /home/user/gurobi:/gurobi loom <TOOL>
```

Replace `<TOOL>` with the desired binary (`loom`, `octi`, etc.) and supply additional command-line flags as needed.

