[![2015 Stuttgart light rail network maps generated from GTFS data, with optimal line orderings, geographically correct (left), octilinear (middle), and  orthoradial (right).](examples/render/stuttgart-example-small.png?raw=true)](examples/render/stuttgart-example.png?raw=true)
*2015 Stuttgart light rail network maps generated from GTFS data, with optimal line orderings, geographically correct (left), octilinear (middle), and  orthoradial (right).*

[![Build](https://github.com/ad-freiburg/loom/actions/workflows/build.yml/badge.svg)](https://github.com/ad-freiburg/loom/actions/workflows/build.yml)

LOOM
====

Software suite for the automated generation of geographically correct or schematic transit maps.

Based on our work in the following papers:
[Bast H., Brosi P., Storandt S., Efficient Generation of Geographically Accurate Transit Maps, SIGSPATIAL 2018](http://ad-publications.informatik.uni-freiburg.de/SIGSPATIAL_transitmaps_2018.pdf)
[Bast H., Brosi P., Storandt S., Efficient Generation of Geographically Accurate Transit Maps (extended version), ACM TSAS, Vol. 5, No. 4, Article 25, 2019](http://ad-publications.informatik.uni-freiburg.de/ACM_efficient%20Generation%20of%20%20Geographically%20Accurate%20Transit%20Maps_extended%20version.pdf)
[Bast H., Brosi P., Storandt S., Metro Maps on Octilinear Grid Graphs, EuroVis 2020](http://ad-publications.informatik.uni-freiburg.de/EuroVis%20octi-maps.pdf)
[Bast H., Brosi P., Storandt S., Metro Maps on Flexible Base Grids, SSTD 2021](http://ad-publications.informatik.uni-freiburg.de/SSTD_Metro%20Maps%20on%20Flexible%20Base%20Grids.pdf).
A pipeline for generating geographically accurate transit maps which appears to be similar to ours was described by Anton Dubrau in a [blog post](https://blog.transitapp.com/how-we-built-the-worlds-prettiest-auto-generated-transit-maps-12d0c6fa502f).

Also see our web demos [here](https://loom.cs.uni-freiburg.de/), [here](https://loom.cs.uni-freiburg.de/global), and [here](https://octi.cs.uni-freiburg.de).

Requirements
------------

 * `cmake`
 * `gcc >= 5.0` (or `clang >= 3.9`)
 * `libicu` development files (e.g., `libicu-dev`)
 * Optional: `libglpk-dev`, `coinor-libcbc-dev`, `gurobi`, `libzip-dev`, `libprotobuf-dev`
 * Optional: `freetype2` (libfreetype) for accurate label rendering


Installing FreeType
-------------------

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


Building and Installation
-------------------------

Fetch this repository and init submodules:

```
git clone --recurse-submodules https://github.com/ad-freiburg/loom.git
```

If you cloned without the `--recurse-submodules` flag (or need to refresh the
dependencies later), make sure the bundled libraries are available before
building:

```
cd loom
git submodule update --init --recursive
```

Build and install:

Make sure the ICU development package is installed before configuring the
project, for example on Debian/Ubuntu:

```bash
sudo apt-get install libicu-dev
```

```
cd loom
mkdir build && cd build
cmake ..
make -j
```

To (optionally) install, type
```
make install
```

You can also use the binaries in `./build` directly.

Usage
=====

This suite consists of several tools:

* `gtfs2graph`, create a GeoJSON line graph from GTFS data
* `topo`, create an overlapping-free line graph from an arbitrary line graph
* `loom`, find optimal line orderings on a line graph
* `octi`, create a schematic version of a line graph
* `transitmap`, render a line graph into an SVG map

All tools output a graph, in the GeoJSON format, to `stdout`, and expect a GeoJSON graph at `stdin`. Exceptions are `gtfs2graph`, where the input is a GTFS feed, and `transitmap`, which writes SVG to `stdout`. Running a tool with `-h` will show a help message with all allowed options.

The `example` folder contains several overlapping-free line graphs.

To render the geographically correct Stuttgart map from above, use
```
cat examples/stuttgart.json | loom | transitmap > stuttgart.svg
```

To also render labels, use

```
cat examples/stuttgart.json | loom | transitmap -l > stuttgart-label.svg
```

To render an *octilinear* map, put the `octi` tool into the pipe:

```
cat examples/stuttgart.json | loom | octi | transitmap -l > stuttgart-octilin.svg
```

To render for example the orthoradial map from above, use a different base graph for `octi`:

```
cat examples/stuttgart.json | loom | octi -b orthoradial | transitmap -l > stuttgart-orthorad.svg
```

Padding can now be specified per side. When any side-specific padding is
provided, remaining sides default to `0`. For example, to add padding only to
the top and left of the map:

```
cat examples/stuttgart.json | loom | transitmap --padding-top 50 --padding-left 20 > stuttgart-pad.svg
```

Configuration
-------------

All tools consult optional `.loom.ini` files for default settings. First
`$HOME/.loom.ini` is loaded, then a `.loom.ini` located next to the running
binary, and finally any command-line flags. You can also specify an explicit
configuration file with `--config=<file>`. Later sources override earlier
ones. See the provided [loom.ini](loom.ini) for supported keys and default
values. The `log-level` key (or `--log-level` flag) controls logging verbosity
from `0` (errors only) to `4` (very verbose debug) and defaults to `2`.

Terminus route label placement can also be configured globally with the
`terminus-label-anchor` key:

* `station-label` (default) – anchor the route label stack to the positioned
  station name label, preserving legacy behavior.
* `stop-footprint` – derive the anchor from the combined polygon footprint of
  the station so the labels clear the rendered outline.
* `node` – fall back to the raw node position for cases without station label
  or footprint geometry.

Station label candidate generation can be tuned through two related options:

* `station-label-angle-steps` controls how many evenly spaced orientations are
  sampled around each station. The value must be positive and divisible by four
  so terminus labels can align horizontally and vertically.
* `station-label-angle-step-deg` sets the angular distance in degrees between
  successive samples. Together, the two settings let you trade off between
  placement flexibility and runtime.

Tool capabilities
-----------------

* `gtfs2graph` – generate a GeoJSON line graph from a GTFS feed.
* `topo` – remove overlapping edges and infer turn restrictions.
* `loom` – compute optimal line orderings for a line graph.
* `octi` – convert a line graph into a schematic layout on a base grid.
* `transitmap` – render a line graph as an SVG map or as MVT vector tiles.

Command-line parameters
-----------------------

### gtfs2graph

* `-m`, `--mots <modes>`: MOTs to calculate shapes for, comma‑separated list of mode names or GTFS codes (default `all`).
* `-p`, `--prune-threshold <0..1>`: threshold for pruning seldomly occurring lines (default `0`).
* `-h`, `--help`: show help message.
* `-v`, `--version`: print version.

### topo

* `-d`, `--max-aggr-dist <meters>`: maximum distance between segments (default `50`).
* `--infer-restr-max-dist <meters>`: maximum distance for turn-restriction checks (default uses `--max-aggr-dist`).
* `--max-comp-dist <meters>`: maximum distance between nodes in a component (default `10000`).
* `--sample-dist <length>`: sample length for map construction in pseudometers (default `5`).
* `--random-colors`: fill missing colors with random values.
* `--write-stats`: write statistics to the output file.
* `--no-infer-restrs`: don't infer turn restrictions.
* `-h`, `--help` and `-v`, `--version`.

### loom

* `--no-untangle`: don't apply untangling rules.
* `--no-prune`: don't apply pruning rules.
* `-m`, `--optim-method <method>`: optimization method (e.g., `ilp`, `comb-no-ilp`; default `comb-no-ilp`).
* `--same-seg-cross-pen <weight>`: penalty for same‑segment crossings (default `4`).
* `--diff-seg-cross-pen <weight>`: penalty for different‑segment crossings (default `1`).
* `--in-stat-cross-pen-same-seg <weight>`: penalty for same‑segment crossings at stations (default `12`).
* `--in-stat-cross-pen-diff-seg <weight>`: penalty for different‑segment crossings at stations (default `3`).
* `--sep-pen <weight>`: penalty for separations (default `3`).
* `--in-stat-sep-pen <weight>`: penalty for separations at stations (default `9`).
* `--ilp-solver <solver>`: preferred ILP solver (`glpk`, `cbc`, or `gurobi`; default `gurobi`).
* `--ilp-num-threads <n>`: number of threads for the ILP solver (`0` for solver default).
* `--ilp-time-limit <sec>`: ILP solve time limit (`-1` for infinite).
* `--output-stats` / `--write-stats`: print or write statistics.
* `--dbg-output-path <path>`: path used for debug output.
* `--output-optgraph`: write optimization graph to the debug path.
* `-h`, `--help` and `-v`, `--version`.

### octi

* `-m`, `--optim-mode <heur|ilp>`: optimization mode (default `heur`).
* `--obstacles <file>`: GeoJSON file containing obstacle polygons.
* `-g`, `--grid-size <len or %>`: grid cell length or percentage of adjacent station distance (default `100%`).
* `-b`, `--base-graph <type>`: base graph (`ortholinear`, `octilinear`, `orthoradial`, `quadtree`, or `octihanan`; default `octilinear`).
* `--retry-on-error`: retry at 85% grid size on error (30 attempts).
* `--skip-on-error`: skip graph on error.
* `--ilp-num-threads <n>` and `--ilp-time-limit <sec>`: ILP solver threads and time limit.
* `--ilp-cache-dir <dir>` and `--ilp-cache-threshold <val>`: ILP cache configuration.
* `--ilp-solver <solver>`: preferred ILP solver (`glpk`, `cbc`, or `gurobi`; default `gurobi`).
* `--hanan-iters <n>`: number of Hanan grid iterations.
* `--loc-search-max-iters <n>`: maximum local search iterations (default `100`).
* `--geo-pen <weight>`: enforce lines to follow input geometry (default `0`).
* `--max-grid-dist <n>`: maximum grid distance for station candidates (default `3`).
* `--edge-order <method>`: initial edge ordering method (e.g., `num-lines`, `length`, `adj-nd-deg`, etc.; default `all`).
* `--density-pen <weight>`: density penalty for local search (default `10`).
* `--vert-pen <weight>`, `--hori-pen <weight>`, `--diag-pen <weight>`: penalties for vertical, horizontal, and diagonal edges.
* `--pen-180 <w>`, `--pen-135 <w>`, `--pen-90 <w>`, `--pen-45 <w>`: penalties for bends.
* `--nd-move-pen <weight>`: penalty for node movement.
* `-h`, `--help` and `-v`, `--version`.

### transitmap

* `--line-width <px>`: width of a single transit line (default `20`).
* `--line-spacing <px>`: spacing between transit lines (default `10`).
* `--outline-width <px>`: width of line outlines (default `1`).
* `--log-level <0..4>`: logging verbosity, `0`=errors to `4`=very verbose (default `2`).
* `--render-dir-markers` and `--render-markers-tail`: render line direction markers and tails.
* `--dir-marker-spacing <n>`: edges between forced direction markers (default `1`).
* `--tail-ignore-sharp-angle`: ignore the sharp-angle check when rendering marker tails (default off).
* `--bi-dir-marker`: render markers for bidirectional edges (default off).
* `--crowded-line-thresh <n>`: lines on edge to trigger direction marker (default `3`).
* `--sharp-turn-angle <rad>`: turn angle in radians (0-π) to trigger direction marker (default `0.785398`). Values >π are treated as degrees.
* `-l`, `--labels`: render labels.
* `-r`, `--route-labels`: render route names at line termini.
* `--line-label-textsize <size>`: text size for line labels (default `40`).
* `--line-label-bend-angle <rad>`: max bend angle in radians for line label candidates (default `0.349066`).
* `--line-label-length-ratio <ratio>`: max length/straight distance ratio for line label candidates (default `1.1`).
* `--station-label-textsize <size>`: text size for station labels (default `60`).
* `--me-label-textsize <size>`: text size for "YOU ARE HERE" label (default `80`).
* `--font-svg-max <size>`: max font size for station labels in SVG, -1 for no limit (default `11`).
* `--station-line-overlap-penalty <weight>`: penalty multiplier for station-line overlaps (default `15`).
* `--station-line-overlap-per-line`: count distinct transit lines when scoring station-line overlaps (default disabled).
* `--station-label-far-crowd-radius <px>`: radius from the far end of a station label used to look for nearby edges, existing labels, station hulls, and landmark boxes (default `0`, disables).
* `--station-label-far-crowd-penalty <weight>`: penalty applied when the far label end crowds any of those nearby features (default `25`).
* `--side-penalty-weight <weight>`: weight for station label side preference penalties (default `2.5`).
* `--same-side-penalty <penalty>`: penalty for station labels on opposite sides (default `100`).
* `--reposition-label <n>`: perform `n` extra passes after initial placement to
  relieve label crowding (default `0`).
* `--cluster-pen-scale <scale>`: scale factor for station crowding penalties (default `1`).
* `--outside-penalty <weight>`: penalty (positive) or bonus (negative) for labels outside the map bounds (default `-5`).
* `--orientation-penalties <p0,...,p7>`: comma-separated penalties for eight label orientations (default `0,3,6,4,1,5,6,2`).
* `--route-label-gap <px>`: gap between route label boxes (default `10`).
* `--route-label-terminus-gap <px>`: gap between terminus station label and route labels (default `80`).
* `--terminus-label-anchor <anchor>`: anchor geometry for terminus route labels
  (`station-label`, `stop-footprint`, or `node`; default `station-label`).
* `--compact-terminal-label`: arrange terminus route labels in multiple columns instead of a single row (default off).
* `--compact-route-label`: stack edge route labels in multiple rows to avoid truncation (default off).
* `--highlight-terminal`: highlight terminus stations (default off).
* `--terminus-highlight-fill <color>`: fill color when highlighting terminus stations (default `black`).
* `--terminus-highlight-stroke <color>`: stroke color when highlighting terminus stations (default `#BAB6B6`).
* `--no-deg2-labels`: suppress labels for degree‑2 stations.
* `-D`, `--from-dot`: input graph is in DOT format.
* `--resolution <res>`: output resolution (default `0.1`).
* `--padding <padding>`: padding around the map (`-1` for auto); applied to all
  sides if no side-specific padding is provided.
* `--padding-top <padding>`: padding at the top of the map (`-1` for auto).
* `--padding-right <padding>`: padding at the right side (`-1` for auto).
* `--padding-bottom <padding>`: padding at the bottom (`-1` for auto).
* `--padding-left <padding>`: padding at the left side (`-1` for auto).

When any side-specific padding is provided, unspecified sides default to `0`.
* `--smoothing <factor>`: input line smoothing (default `1`).
* `--ratio <value>`: output width/height ratio (`width = height * ratio`).
* `--tl-ratio <value>`: top-left anchored width/height ratio (default `-1`).
* `--random-colors`: fill missing colors with random colors.
* `--tight-stations`: don't expand node fronts for stations.
* `--no-render-stations`: don't render stations.
* `--no-render-node-connections`: don't render inner node connections.
* `--render-node-fronts`: render node fronts.
* `--bg-map <file.geojson>`: render additional GeoJSON geometry behind the network. Coordinates
  are expected in latitude/longitude (WGS84). Polygon and MultiPolygon features are supported and
  feature `properties` may specify `stroke`, `stroke-width`, `fill`, and `opacity` for styling.
* `--bg-map-webmerc`: treat `--bg-map` coordinates as already in Web Mercator and skip conversion.
* `--bg-map-opacity <value>`: opacity for `--bg-map` geometry (default `1`).
* `--extend-with-bgmap`: expand the output bounding box to include `--bg-map` geometry.
* `--geo-lock [bool]`: ensure the map covers at least a default geographic bounding box (defaults to `true`).
* `--geo-lock-bbox <south,west,north,east>`: custom bounding box for `--geo-lock` (latitude/longitude).
* `--landmark <spec>`: add a landmark `word:text,lat,lon[,fontSize[,color[,opacity]]]` or
  `iconPath,lat,lon[,size]`.
* `--landmarks <file>`: read landmarks from a file, one per line.
* `--force-landmarks`: keep landmarks even when they overlap existing geometry (default `true`).
                        Passing `--force-landmarks=false` turns on the displacement search for icon
                        landmarks, nudging them within the configured radius and invoking the
                        optional skip path for entries that have nothing drawable; icons that still
                        collide after the search are rendered in the best position the solver
                        found, so overlaps can remain.
* `--landmark-search-radius <radius>`: search radius for shifting overlapping
  landmark icons (default `10`).
* `--displacement-iterations <n>`: maximum iterations for landmark
  displacement (default `100`).
* `--displacement-cooling <factor>`: cooling factor for displacement step
  (default `0.9`).
* `--landmarks-webmerc`: treat landmark and `--me` coordinates as already in
                          Web Mercator and skip conversion.
* `--me <lat,lon>`: mark the given coordinates with a red star (latitude and
                      longitude by default).
* `--me-size <size>`: star size (default `150`).
* `--me-label`: add a "YOU ARE HERE" label.
* `--me-station <name>`: mark current location by station label.
* `--me-with-bg[=<bool>]`: when `--me-station` is active, restyle the matched
  station label with a rounded horizontal badge and star; falls back to the
  standalone badge when the label cannot be restyled.
* `--me-bg-fill <color>`: badge fill color for the background behind the `--me`
  marker (default `#f5f5f5`).
* `--me-bg-stroke <color>`: badge stroke color for the `--me` background
  (default `#d0d0d0`).
* `--me-label-color <color>`: text color used for the badge label (default
  `#3a3a3a`).
* `--me-station-fill <color>`: fill color for "me" marker (default `#f00`).
* `--me-station-border <color>`: border color for "me" marker (default none).

Enabling `--me-with-bg` keeps the solver-selected station label in place when
labels are rendered and decorates it with the horizontal badge: the star sits to
the left of the name while the configured background fill and stroke frame the
text. The badge width adapts automatically to the station label so longer names
receive more space, and `--me-label-textsize` still controls the font size used
for the text. The badge now mirrors the asymmetric vertical padding of the
terminus route label boxes (28% above the content, 12% below) so the highlight
badge and terminus labels share the same slimmer silhouette. When no station
label is available (for example when `--labels` is omitted or the station could
not be matched), the standalone badge layout from earlier releases is drawn
instead so existing workflows continue to work unchanged. Solver collision
avoidance now reserves the entire badge footprint—including the star and badge
padding—so neighbouring labels and landmarks respect the combined space.
* `--print-stats`: write statistics to stdout.
* `-h`, `--help` and `-v`, `--version`.

Line graph extraction from GTFS
-------------------------------

To extract for example a line graph for streetcars from GTFS data, use `gtfs2graph` as follows:

```
gtfs2graph -m tram freiburg.zip > freiburg.json
```

This line graph will have many overlapping edges and stations. To create an overlapping-free line graph ready for rendering, add `topo`:

```
gtfs2graph -m tram freiburg.zip | topo > freiburg.json
```

A full pipeline for creating an octilinear map of the Freiburg tram network would look like this:
```
gtfs2graph -m tram freiburg.zip | topo | loom | octi | transitmap > freiburg-tram.svg
```

### GeoJSON output attributes

`gtfs2graph` exports each stop as a GeoJSON point with a `properties` map. In addition
to the station identifier and label, the `terminals` array now lists the GTFS routes
for which that node is a terminus (including their ID, short name, and color). Downstream
consumers can use this metadata to highlight terminal stops without re-deriving the
endpoints from the network topology.

A sample landmarks file with matching SVG icons is provided in
`examples/landmarks.txt`.

Each landmark line is either

* `word:<text>,lat,lon[,fontSize[,color[,opacity]]]` – render the given text at the
  latitude and longitude position. The optional `fontSize` (in pixels) defaults
  to `20`, `color` to `#000`, and `opacity` to `1`. If you want to specify only a
  color, omit the font size, e.g. `word:CityHall,47.92,106.91,#ff0000`.
* `iconPath,lat,lon[,size]` – place an SVG icon from `iconPath`. The optional
  `size` also defaults to `200`. Relative `iconPath` values are resolved
  relative to the landmarks file. If the icon file cannot be read, a warning is
  logged and the landmark is skipped.

Landmark coordinates are interpreted as WGS84 latitude/longitude and converted
to Web Mercator automatically. Use `--landmarks-webmerc` if the coordinates are
already given in Web Mercator.

Landmarks that would overlap with existing labels, map features, or previously
placed landmarks are repositioned using a small force-directed search. Use
`--force-landmarks=false` to skip those overlapping landmarks instead of moving
them.
To render the sample landmarks alongside a map, run

```
cat examples/stuttgart.json | loom | transitmap --landmarks examples/landmarks.txt > stuttgart-landmarks.svg
```

You can also add a single word landmark directly from the command line:

```
cat examples/stuttgart.json | loom | transitmap --landmark word:CityHall,47.9210,106.9175,150,#ff0000 > stuttgart-cityhall.svg
```

Landmarks for additional cities can be created in the same way. The sample
file uses locations around Ulaanbaatar, Mongolia; human-readable names are
listed in `examples/ulaanbaatar_landmark_names.txt`.

Usage via Docker
================

You can also use any tool in a Docker container via the provided Dockerfile.

To build the container:

```
docker build -t loom
```

To run a tool from the suite, use

```
docker run -i loom <TOOL>
```

For example, to octilinearize the Freiburg example, use

```
cat examples/freiburg.json | sudo docker run -i loom octi
```

*Note*: if you want to use gurobi for ILP optimization, you *must* mount a folder container a valid gurobi license file `gurobi.lic` to `/gurobi/` in the container. For example, if your `gurobi.lic` is in `/home/user/gurobi`:

```
docker run -v /home/user/gurobi:/gurobi loom <TOOL>
```
