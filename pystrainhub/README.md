# pystrainhub

A Python implementation of the **StrainHub** phylogenetic network toolkit.  
Generates pathogen transmission networks from genomic / phylogenetic data and
computes graph-theoretic centrality metrics to rank the importance of each
location (or host) in the transmission network.

This library is a pure-Python equivalent of the
[StrainHub R package](../pkg) — it uses Python libraries directly, **not** an
R-to-Python wrapper.

---

## Installation

```bash
# Core (Newick + dataframe tree types)
pip install .

# With interactive visualization (pyvis + folium)
pip install ".[visualization]"

# With Bayesian tree support (dendropy)
pip install ".[bayesian]"

# Everything
pip install ".[all]"
```

Or install directly from the requirements file during development:

```bash
pip install -r requirements.txt
pip install -e .
```

---

## Quick start

```python
from Bio import Phylo
import pandas as pd
from pystrainhub import list_states, make_transnet, make_map

# --- Load data ---
tree     = Phylo.read("data/parsimonious/chikv/chikv_westernafrica.phy", "newick")
metadata = pd.read_csv("data/parsimonious/chikv/chikv_westernafrica_metadata.csv")
geodata  = pd.read_csv("data/parsimonious/chikv/chikv_geo.csv")

# --- Inspect available columns ---
list_states(tree, metadata, tree_type="parsimonious")
#    Index    Column
# 0      1  Accession
# 1      2      Host
# 2      3   Country
# ...

# --- Build transmission network ---
graph = make_transnet(
    tree,
    metadata,
    column_selection="Country",
    centrality_metric=6,          # Source Hub Ratio
    tree_type="parsimonious",
    metrics_output_file="StrainHub_metrics.csv",
)

# Centrality metrics DataFrame
print(graph["metrics"])

# Interactive network (requires pyvis)
graph["visualization"].show("network.html")

# --- Map (requires folium) ---
m = make_map(
    graph,
    geodata,
    column_selection="Country",
    show_labels=True,
    show_points=True,
)
m.save("map.html")
```

---

## Supported tree types

| `tree_type`     | Input `treedata`                                  | Extra required           |
|-----------------|--------------------------------------------------|--------------------------|
| `"parsimonious"`| Newick file path or `Bio.Phylo` tree object      | `metadata` DataFrame     |
| `"nj"`          | FASTA file path or `Bio.AlignIO` alignment object| `metadata` DataFrame     |
| `"bayesian"`    | BEAST nexus file path or `dendropy.Tree`         | –                        |
| `"dataframe"`   | `pd.DataFrame` with `from`/`to`[/`value`] cols  | –                        |

---

## Centrality metrics

| Value | Metric                  |
|-------|-------------------------|
| 0     | All metrics (CSV only)  |
| 1     | Indegree centrality     |
| 2     | Outdegree centrality    |
| 3     | Betweenness centrality  |
| 4     | Closeness centrality    |
| 5     | Degree centrality       |
| 6     | **Source Hub Ratio**    |

The **Source Hub Ratio** (outdegree / degree) ranges from ~0 (sink) to
0.5 (hub) to ~1 (source), and is the primary StrainHub metric for identifying
transmission sources.

---

## API reference

### `list_states(treedata, metadata, tree_type)`
Returns a `pd.DataFrame` listing all character-state columns available for
network construction.

### `make_transnet(treedata, metadata, column_selection, centrality_metric, ...)`
Builds a directed transmission network and returns a `dict` containing:

* `"nodes"` – `pd.DataFrame(id, label)`
* `"edges"` – `pd.DataFrame(from, to, value)`
* `"metrics"` – `pd.DataFrame` of six centrality values per node
* `"graph"` – `networkx.DiGraph`
* `"visualization"` – `pyvis.network.Network` (or `None`)

Full parameter reference: see docstring in `pystrainhub/strainhub.py`.

### `make_map(graph, geodata, column_selection, ...)`
Creates a *Folium* interactive Leaflet map of the transmission network with
coloured directed arrows between locations.

---

## Running tests

```bash
pip install ".[dev]"
pytest tests/ -v
```

---

## Citation

If you use StrainHub or pystrainhub in your research, please cite:

> Adriano de Bernardi Schneider, Colby T Ford, Reilly Hostager, John Williams,
> Michael Cioce, Ümit V Çatalyürek, Joel O Wertheim, Daniel Janies.
> **StrainHub: A phylogenetic tool to construct pathogen transmission networks.**
> *Bioinformatics*, btz646. <https://doi.org/10.1093/bioinformatics/btz646>
