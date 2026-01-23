LOOM Toolchain
==============

Automated pipeline for generating geographically correct or schematic transit maps from line graphs, with optional stop merging and SVG rendering. This README walks you through the full end-to-end pipeline without requiring code inspection or `-h` usage.

Overview
--------

LOOM is a suite of small CLI tools that stream GeoJSON line graphs through a map-production pipeline. The typical flow is:

```
input.json | topo | stopmerge | loom | octi | transitmap > out.svg
```

Quick Start
-----------

Minimal working pipeline on the included example:

```
cat examples/ub.json \
  | build/topo \
  | build/stopmerge --merge-stops always \
  | build/loom \
  | build/octi \
  | build/transitmap \
      --canvas-width 3000 \
  > /tmp/ub_map.svg
```

Toolchain at a Glance
---------------------

```
GeoJSON line graph
    |
  topo        (resolve overlaps / topology cleanup)
    |
  stopmerge   (optional nearby stop merging + metadata)
    |
  loom        (optimize line ordering)
    |
  octi        (schematize to base graph)
    |
  transitmap  (render SVG)
    |
  SVG output
```

Build / Install (if applicable)
-------------------------------

Requirements:

* `cmake`
* `gcc >= 5.1` (or `clang >= 3.9`)
* Optional: `libglpk-dev` (ILP solver), `coinor-libcbc-dev` (ILP solver), `gurobi` (ILP solver), `libzip-dev`

Build and install:

```
git clone --recurse-submodules https://github.com/ad-freiburg/loom.git
cd loom
mkdir build && cd build
cmake ..
make -j
```

Optional install:

```
make install
```

You can also use the binaries directly from `build/`.

Inputs and Outputs
------------------

All tools (except `gtfs2graph` and `transitmap`) read a GeoJSON line graph from `stdin` and write a GeoJSON line graph to `stdout`.

* `gtfs2graph`: input is a GTFS feed, output is a GeoJSON line graph.
* `transitmap`: input is a GeoJSON line graph, output is an SVG map.

At a high level, the initial JSON represents a transit line graph with stations and edges; downstream tools add ordering, schematization, and rendering data.

Pipeline Walkthrough
--------------------

Minimal pipeline (no background):

```
cat examples/ub.json \
  | build/topo \
  | build/stopmerge --merge-stops always \
  | build/loom \
  | build/octi \
  | build/transitmap \
      --canvas-width 3000 \
  > /tmp/ub_map.svg
```

With MBTiles background (requires `--mbtiles` and the example file):

```
cat examples/ub.json \
  | build/topo \
  | build/stopmerge --merge-stops always \
  | build/loom \
  | build/octi \
  | build/transitmap \
      --canvas-width 3000 \
      --mbtiles examples/ub.mbtiles \
  > /tmp/ub_map_mbtiles.svg
```

Intermediate JSON dumps:

```
cat examples/ub.json | build/topo | build/stopmerge --merge-stops always > /tmp/after_stopmerge.json
cat examples/ub.json | build/topo | build/stopmerge --merge-stops always | build/loom > /tmp/after_loom.json
cat examples/ub.json | build/topo | build/stopmerge --merge-stops always | build/loom | build/octi > /tmp/after_octi.json
```

Verify stopmerge metadata persistence:

```
rg -n "stopmerge_member_count|stopmerge_is_merged_rep" /tmp/after_stopmerge.json | head
rg -n "stopmerge_member_count|stopmerge_is_merged_rep" /tmp/after_loom.json | head
rg -n "stopmerge_member_count|stopmerge_is_merged_rep" /tmp/after_octi.json | head
```

Tool Reference
--------------

topo
~~~~

Purpose: clean up line graphs by resolving overlaps and producing a topology suitable for downstream optimization and rendering.

Typical usage:

```
cat input.json | build/topo > topo.json
```

Flags:

```
-v, --version
-h, --help
-d, --max-aggr-dist (=50)
--write-stats
--no-infer-restrs
--infer-restr-max-dist (=[-d])
--max-comp-dist (=10000)
--sample-dist (=5)
--max-length-dev (=500)
--turn-restr-full-turn-angle (=0)
--turn-restr-full-turn-pen (=0)
--random-colors
--write-components
--write-components-path
--smooth (=0)
--aggr-stats
```

stopmerge
~~~~~~~~~

Purpose: optionally merge nearby stops/stations and add merged-stop metadata to nodes.

Typical usage:

```
cat topo.json | build/stopmerge --merge-stops auto > merged.json
```

Flags:

```
--merge-stops (=never)                  never|auto|always
--merge-stops-parallel-corridors (=off) off|auto|always
--merge-stops-snap-dist-m (=40)
--merge-stops-radius-m (=40)
--merge-stops-chainage-m (=40)
--merge-stops-node-zone-m (=15)
--merge-stops-corridor-angle-deg (=25)
--merge-stops-max-cluster-size (=6)
--parallel-pair-max-dist-m (=25)
--parallel-pair-max-angle-deg (=20)
--parallel-pair-min-overlap-ratio (=0.35)
--parallel-pair-require-lines (=true)
--parallel-pair-min-group-length-m (=150)
--parallel-pair-max-median-dist-m (=30)
--hub-radius-m (=80)
--merge-stops-max-local-density (=8)
--hub-auto-require-base-name (=true)
--merge-stops-hub-fallback (=auto)      off|auto|always
--merge-stops-debug-csv
--merge-stops-debug-rejects-csv
--parallel-pair-debug-csv
--stop-merge-map-json
--merge-stops-overrides
-v, --version
-h, --help
```

loom
~~~~

Purpose: optimize line orderings to reduce crossings and improve readability.

Typical usage:

```
cat input.json | build/loom > ordered.json
```

Flags:

```
-v, --version
-h, --help
--no-untangle
--no-prune
-m, --optim-method (=comb-no-ilp)
--same-seg-cross-pen (=4)
--diff-seg-cross-pen (=1)
--in-stat-cross-pen-same-seg (=12)
--in-stat-cross-pen-diff-seg (=3)
--sep-pen (=3)
--in-stat-sep-pen (=9)
-D, --from-dot
--output-stats
--write-stats
--ilp-solver (=gurobi)
--ilp-num-threads (=0)
--ilp-time-limit (=-1)
--dbg-output-path (=.)
--output-optgraph
```

octi
~~~~

Purpose: schematize the graph onto a base grid (octilinear by default).

Typical usage:

```
cat input.json | build/octi > octi.json
```

Flags:

```
-v, --version
-h, --help
-m, --optim-mode (=heur)                heur|ilp
--obstacles
-g, --grid-size (=100%)
-b, -base-graph (=octilinear)           ortholinear|octilinear|orthoradial|quadtree|octihanan
--retry-on-error
--skip-on-error
--ilp-num-threads (=0)
--hanan-iters (=1)
--loc-search-max-iters (=100)
--ilp-cache-threshold (=inf)
--ilp-time-limit (=60)
--ilp-cache-dir (=.)
--ilp-solver (=gurobi)                  glpk|cbc|gurobi
--write-stats
-D, --from-dot
--no-deg2-heur
--geo-pen (=0)
--max-grid-dist (=3)
--restr-loc-search
--edge-order (=all)
--density-pen (=10)
--vert-pen (=0)
--hori-pen (=0)
--diag-pen (=.5)
--pen-180 (=0)
--pen-135 (=1)
--pen-90 (=1.5)
--pen-45 (=2)
--nd-move-pen (=.5)
```

transitmap
~~~~~~~~~~

Purpose: render a line graph into SVG.

Typical usage:

```
cat input.json | build/transitmap --canvas-width 3000 > map.svg
```

Flags:

```
-v, --version
-h, --help
--render-engine (=svg)
--line-width (=20)
--line-spacing (=10)
--outline-width (=1)
--render-dir-markers
-l, --labels
--line-label-textsize (=40)
--station-label-textsize (=60)
--no-deg2-labels
--show-merged-stop-members (=tooltip)   off|tooltip|inline|multiline
--merged-stop-view (=double_ring)       none|double_ring|merge_count|simple_cross
--station-radius (=-1)
--station-radius-mult (=1)
--merged-cross-scale (=1.1)
--merged-cross-stroke-mult (=1.4)
--debug-stop-labels
-D, --from-dot
--padding (=-1)
--padding-top (=-1)
--padding-right (=-1)
--padding-bottom (=-1)
--padding-left (=-1)
--smoothing (=1)
--random-colors
--tight-stations
--no-render-stations
--no-render-node-connections
--render-node-fronts
--print-stats
--mbtiles
--paper (=A4L)                          A4|A4L|A3|A3L
--canvas-width (=1000)
--canvas-unit (=px)                     px|mm
--zoom-levels
--max-tiles (=512)
--oversample (=1.25)
--background-pad-pct (=0.03)
--background-opacity (=1.0)
```

Examples Cookbook
-----------------

Render basic SVG without background:

```
cat examples/ub.json \
  | build/topo \
  | build/stopmerge --merge-stops always \
  | build/loom \
  | build/octi \
  | build/transitmap --canvas-width 3000 \
  > /tmp/ub_map.svg
```

Render with background MBTiles:

```
cat examples/ub.json \
  | build/topo \
  | build/stopmerge --merge-stops always \
  | build/loom \
  | build/octi \
  | build/transitmap --canvas-width 3000 --mbtiles examples/ub.mbtiles \
  > /tmp/ub_map_mbtiles.svg
```

Dump intermediate JSON:

```
cat examples/ub.json | build/topo | build/stopmerge --merge-stops always > /tmp/after_stopmerge.json
cat examples/ub.json | build/topo | build/stopmerge --merge-stops always | build/loom > /tmp/after_loom.json
cat examples/ub.json | build/topo | build/stopmerge --merge-stops always | build/loom | build/octi > /tmp/after_octi.json
```

Show merged stop visuals:

```
cat /tmp/after_octi.json \
  | build/transitmap \
      --merged-stop-view simple_cross \
      --canvas-width 3000 \
  > /tmp/ub_cross.svg
```

Show merged stop members inline:

```
cat /tmp/after_octi.json \
  | build/transitmap \
      --show-merged-stop-members inline \
      --canvas-width 3000 \
  > /tmp/ub_inline.svg
```

Tune station and merged-cross sizing:

```
cat /tmp/after_octi.json \
  | build/transitmap \
      --station-radius 30 \
      --station-radius-mult 1.2 \
      --merged-cross-scale 1.5 \
      --merged-cross-stroke-mult 2.0 \
      --canvas-width 3000 \
  > /tmp/ub_tuned.svg
```

Stopmerge debug outputs:

```
cat examples/ub.json \
  | build/topo \
  | build/stopmerge \
      --merge-stops always \
      --merge-stops-debug-csv /tmp/merge_debug.csv \
      --merge-stops-debug-rejects-csv /tmp/merge_rejects.csv \
      --parallel-pair-debug-csv /tmp/parallel_pairs.csv \
      --stop-merge-map-json /tmp/stop_merge_map.json \
  > /tmp/ub_after_stopmerge.json
```

Debugging and Troubleshooting
-----------------------------

* Pipeline is streaming; each tool reads `stdin` and writes `stdout`. Use intermediate dumps to isolate problems.
* If merged stop metadata is missing later in the pipeline, confirm it appears after `stopmerge` and persists through `loom` and `octi` with the `rg` checks shown above.
* If background rendering fails, confirm the MBTiles path and that `--mbtiles` is set.
* If output is unexpectedly sparse, confirm input GeoJSON is valid and contains stations/edges.

FAQ
---

Q: Where is `linegraph.json` produced?
A: It is produced by the previous stage in the pipeline; each tool writes its output graph to `stdout`. Capture with `> file.json` between stages.

Q: How do I confirm stopmerge metadata persists through loom and octi?
A: Use the `rg` checks in the Pipeline Walkthrough to verify `stopmerge_member_count` and `stopmerge_is_merged_rep` appear after `stopmerge`, `loom`, and `octi`.

Q: How do I render merged stop indicators?
A: Use `transitmap` with `--merged-stop-view` and `--show-merged-stop-members` (see Examples Cookbook).

Usage via Docker
================

You can also use any tool in a Docker container via the provided Dockerfile.

To build the container:

```
docker build -t loom .
```

To run a tool from the suite, use

```
docker run -i loom <TOOL>
```

For example, to octilinearize the Freiburg example, use

```
cat examples/freiburg.json | sudo docker run -i loom octi
```

*Note*: if you want to use gurobi for ILP optimization, you must mount a folder containing a valid `gurobi.lic` to `/gurobi/` in the container. For example, if your `gurobi.lic` is in `/home/user/gurobi`:

```
docker run -v /home/user/gurobi:/gurobi loom <TOOL>
```
