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

Build and install:

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
* `transitmap`, render a line graph into an SVG map (`--render-engine=svg`) or into vector tiles (`--render-engine=mvt`)

All tools output a graph, in the GeoJSON format, to `stdout`, and expect a GeoJSON graph at `stdin`. Exceptions are `gtfs2graph`, where the input is a GTFS feed, and `transitmap`, which writes SVG to `stdout` or MVT vector tiles to a specified folder. Running a tool with `-h` will show a help message with all allowed options.

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

Padding can now be specified per side. For example, to add padding only to the
top and left of the map:

```
cat examples/stuttgart.json | loom | transitmap --padding 0 --padding-top 50 --padding-left 20 > stuttgart-pad.svg
```

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

* `--render-engine <svg|mvt>`: render engine (default `svg`).
* `--line-width <px>`: width of a single transit line (default `20`).
* `--line-spacing <px>`: spacing between transit lines (default `10`).
* `--outline-width <px>`: width of line outlines (default `1`).
* `--render-dir-markers` and `--render-markers-tail`: render line direction markers and tails.
* `-l`, `--labels`: render labels.
* `-r`, `--route-labels`: render route names at line termini.
* `--line-label-textsize <size>`: text size for line labels (default `40`).
* `--station-label-textsize <size>`: text size for station labels (default `60`).
* `--no-deg2-labels`: suppress labels for degree‑2 stations.
* `--zoom <levels>` and `--mvt-path <dir>`: zoom levels and output path for MVT tiles.
* `-D`, `--from-dot`: input graph is in DOT format.
* `--padding <padding>`: padding around the map (`-1` for auto).
* `--padding-top <padding>`: padding at the top of the map (`-1` for auto).
* `--padding-right <padding>`: padding at the right side (`-1` for auto).
* `--padding-bottom <padding>`: padding at the bottom (`-1` for auto).
* `--padding-left <padding>`: padding at the left side (`-1` for auto).
* `--smoothing <factor>`: input line smoothing (default `1`).
* `--ratio <value>`: output width/height ratio (`width = height * ratio`).
* `--random-colors`: fill missing colors with random colors.
* `--tight-stations`: don't expand node fronts for stations.
* `--no-render-stations`: don't render stations.
* `--no-render-node-connections`: don't render inner node connections.
* `--render-node-fronts`: render node fronts.
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

A sample landmarks file with matching SVG icons is provided in
`examples/landmarks.txt`. The landmark size defaults to 200 units when
not specified. To render these landmarks alongside a map, run

```
cat examples/stuttgart.json | loom | transitmap --landmarks examples/landmarks.txt > stuttgart-landmarks.svg
```

Landmarks for additional cities can be created in the same way. An example
dataset for Ulaanbaatar, Mongolia is provided in
`examples/ulaanbaatar_landmarks.txt` with accompanying human-readable names in
`examples/ulaanbaatar_landmark_names.txt`.

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
