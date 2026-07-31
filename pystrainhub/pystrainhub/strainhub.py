"""
pystrainhub – Python implementation of the StrainHub phylogenetic network toolkit.

Three public functions mirror the exported R package API:

* :func:`list_states`  – inspect available character-state columns.
* :func:`make_transnet` – build a transmission network and compute centrality
  metrics.
* :func:`make_map`      – render the network on an interactive Leaflet map
  (via *folium*).

Supported tree types
--------------------
``"parsimonious"``
    Newick tree + metadata CSV.  Ancestral states are inferred with Fitch
    parsimony (equivalent to ``castor::asr_max_parsimony``).
``"bayesian"``
    BEAST MCC nexus tree parsed with *dendropy*.  Node annotations supply
    state labels and probabilities.
``"nj"``
    Neighbor-joining tree built from a FASTA DNA alignment with *Biopython*.
    Ancestral states are then inferred with Fitch parsimony.
``"dataframe"``
    A :class:`pandas.DataFrame` with ``from`` / ``to`` columns (and an
    optional ``value`` column for pre-counted edge weights).
"""

from __future__ import annotations

import json
import math
import warnings
from typing import Optional, Union

import pandas as pd
import networkx as nx

from .parsimony import fitch_parsimony, extract_state_changes

# Optional heavy dependencies – imported lazily so the package still loads
# when they are absent.
try:
    from Bio import Phylo, AlignIO
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

try:
    import pyvis.network as pvnet

    HAS_PYVIS = True
except ImportError:
    HAS_PYVIS = False

try:
    import folium

    HAS_FOLIUM = True
except ImportError:
    HAS_FOLIUM = False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def list_states(
    treedata,
    metadata: Optional[pd.DataFrame] = None,
    tree_type: str = "parsimonious",
) -> pd.DataFrame:
    """
    Return a DataFrame of available character-state columns.

    Parameters
    ----------
    treedata:
        Tree object or file path.  For ``"bayesian"`` this must be either a
        path to a BEAST nexus file or a *dendropy* ``Tree`` object.  For other
        tree types the value is not inspected.
    metadata : pd.DataFrame, optional
        Required for ``"parsimonious"`` and ``"nj"`` tree types.  The first
        column should be ``"Accession"``; every column is listed.
    tree_type : str
        One of ``"parsimonious"``, ``"bayesian"``, ``"nj"``.

    Returns
    -------
    pd.DataFrame
        Columns ``"Index"`` (1-based) and ``"Column"`` (column name).
    """
    if tree_type in ("parsimonious", "nj"):
        if metadata is None:
            raise ValueError(
                f"metadata is required for tree_type='{tree_type}'"
            )
        cols = list(metadata.columns)
        return pd.DataFrame({"Index": range(1, len(cols) + 1), "Column": cols})

    if tree_type == "bayesian":
        try:
            import dendropy
            import re
        except ImportError as exc:
            raise ImportError(
                "dendropy is required for tree_type='bayesian'. "
                "Install with:  pip install dendropy"
            ) from exc

        tree = (
            dendropy.Tree.get(path=treedata, schema="nexus")
            if isinstance(treedata, str)
            else treedata
        )

        ann_names: set[str] = set()
        for node in tree.preorder_node_iter():
            for ann in node.annotations:
                ann_names.add(ann.name)

        # Mirror the R filter: drop HPD / median / range / prob / set / rate
        filtered = sorted(
            name
            for name in ann_names
            if not re.search(
                r"_0\.95_HPD|_median|_range|\.prob|\.set|\.rate", name
            )
        )
        return pd.DataFrame(
            {"Index": range(1, len(filtered) + 1), "Column": filtered}
        )

    raise ValueError(f"Unknown tree_type: '{tree_type}'")


def make_transnet(
    treedata,
    metadata: Optional[pd.DataFrame] = None,
    column_selection: Optional[str] = None,
    centrality_metric: int = 0,
    threshold: float = 0.9,
    threshold2: float = 0.9,
    bootstrap_value: Optional[int] = None,
    tree_type: str = "parsimonious",
    root_selection: Optional[str] = None,
    metrics_output_file: str = "StrainHub_metrics.csv",
    as_json: bool = False,
) -> Union[dict, str]:
    """
    Build a directed transmission network and compute centrality metrics.

    Parameters
    ----------
    treedata:
        * ``"parsimonious"`` / ``"nj"`` – :class:`Bio.Phylo.BaseTree.Tree`
          object **or** path to a Newick / FASTA file.
        * ``"bayesian"`` – path to a BEAST nexus file or a *dendropy* ``Tree``.
        * ``"dataframe"`` – :class:`pandas.DataFrame` with columns ``from``,
          ``to``, and optionally ``value``.
    metadata : pd.DataFrame, optional
        Required for ``"parsimonious"`` and ``"nj"``.  Must contain an
        ``"Accession"`` column.
    column_selection : str
        Metadata column (or BEAST annotation name) defining the character
        state used to build the network.
    centrality_metric : int
        | 0 – compute all metrics (no network visualization returned)
        | 1 – indegree centrality
        | 2 – outdegree centrality
        | 3 – betweenness centrality
        | 4 – closeness centrality
        | 5 – degree centrality
        | 6 – Source Hub Ratio
    threshold : float
        Minimum state-probability for Bayesian node-state acceptance (0.9).
    threshold2 : float
        Minimum posterior probability for Bayesian edge acceptance (0.9).
    bootstrap_value : int, optional
        For ``"nj"`` trees – collapse nodes with bootstrap support below this
        value (0–100).
    tree_type : str
        ``"parsimonious"``, ``"bayesian"``, ``"nj"``, or ``"dataframe"``.
    root_selection : str, optional
        Accession name to use as the outgroup root for ``"nj"`` trees.
    metrics_output_file : str
        Path for the output CSV of centrality metrics.  Pass ``""`` to skip.
    as_json : bool
        If ``True`` return a JSON string; otherwise return a ``dict``.

    Returns
    -------
    dict or str
        ``dict`` with keys:

        * ``"nodes"``         – :class:`pandas.DataFrame` (``id``, ``label``)
        * ``"edges"``         – :class:`pandas.DataFrame` (``from``, ``to``,
          ``value``)
        * ``"metrics"``       – :class:`pandas.DataFrame` of centrality values
        * ``"graph"``         – :class:`networkx.DiGraph`
        * ``"visualization"`` – :class:`pyvis.network.Network` (or ``None`` if
          *pyvis* is unavailable or ``centrality_metric == 0``)

        When *as_json* is ``True`` a JSON string is returned instead.
    """
    if tree_type == "parsimonious":
        edges_df, nodes_df = _process_parsimonious(
            treedata, metadata, column_selection
        )
    elif tree_type == "bayesian":
        edges_df, nodes_df = _process_bayesian(
            treedata, column_selection, threshold, threshold2
        )
    elif tree_type == "nj":
        edges_df, nodes_df = _process_nj(
            treedata, metadata, column_selection, root_selection, bootstrap_value
        )
    elif tree_type == "dataframe":
        edges_df, nodes_df = _process_dataframe(treedata)
    else:
        raise ValueError(f"Unknown tree_type: '{tree_type}'")

    # Build NetworkX directed graph for centrality calculations
    G: nx.DiGraph = nx.DiGraph()
    for _, row in nodes_df.iterrows():
        G.add_node(int(row["id"]))
    for _, row in edges_df.iterrows():
        G.add_edge(int(row["from"]), int(row["to"]), weight=int(row["value"]))

    # Compute centrality metrics
    metastates = list(nodes_df["label"])
    metrics = _compute_centrality_metrics(G, metastates)

    if metrics_output_file:
        metrics.to_csv(metrics_output_file, index=False)

    # Build optional visualization
    vis = (
        _build_pyvis_network(nodes_df, edges_df, centrality_metric, metrics)
        if centrality_metric in range(1, 7)
        else None
    )

    result: dict = {
        "nodes": nodes_df,
        "edges": edges_df,
        "metrics": metrics,
        "graph": G,
        "visualization": vis,
    }

    if as_json:
        return json.dumps(
            {
                "nodes": nodes_df.to_dict(orient="records"),
                "edges": edges_df.to_dict(orient="records"),
                "metrics": metrics.to_dict(orient="records"),
            }
        )

    return result


def make_map(
    graph: dict,
    geodata: pd.DataFrame,
    column_selection: str,
    basemap_layer: str = "CartoDB positron",
    hide_arrow_head: bool = False,
    arrow_filled: bool = True,
    show_labels: bool = False,
    label_color: str = "#000000",
    show_points: bool = False,
    point_color: str = "#000000",
    point_opacity: float = 0.5,
):
    """
    Create a *Folium* interactive map of the transmission network.

    Each directed edge is drawn as a coloured line with an arrow indicator.
    Location labels and circle markers can be added optionally.

    Parameters
    ----------
    graph : dict
        Output of :func:`make_transnet`.
    geodata : pd.DataFrame
        Geographic reference table.  Must contain a column matching
        *column_selection* plus ``"Latitude"`` and ``"Longitude"`` columns.
    column_selection : str
        Column name that matches the ``label`` column in *graph["nodes"]* and
        the location column in *geodata* (e.g. ``"Country"``).
    basemap_layer : str
        *Folium* tile-layer name (default ``"CartoDB positron"``).  Any tile
        name accepted by :func:`folium.Map` works, e.g. ``"OpenStreetMap"``.
    hide_arrow_head : bool
        Suppress directional arrow indicators (default ``False``).
    arrow_filled : bool
        Use a filled arrowhead character (``▶``); set to ``False`` for an
        outline arrow (``›``) (default ``True``).
    show_labels : bool
        Overlay location name labels on the map (default ``False``).
    label_color : str
        CSS colour for labels (default ``"#000000"``).
    show_points : bool
        Draw a filled circle at each location (default ``False``).
    point_color : str
        CSS colour for circles (default ``"#000000"``).
    point_opacity : float
        Fill opacity for circles, 0–1 (default ``0.5``).

    Returns
    -------
    folium.Map
    """
    if not HAS_FOLIUM:
        raise ImportError(
            "folium is required for make_map(). Install with:  pip install folium"
        )

    nodes_df = graph["nodes"]
    edges_df = graph["edges"]

    # Join nodes to geodata so each node gets lat/lon
    locations = nodes_df.merge(
        geodata, left_on="label", right_on=column_selection, how="inner"
    )
    if locations.empty:
        raise ValueError(
            "No matching rows found when joining network nodes to geodata. "
            f"Check that geodata has a column named '{column_selection}'."
        )

    lat_col, lon_col = "Latitude", "Longitude"

    # Build per-edge DataFrame with source and destination coordinates
    graph_df = (
        edges_df.copy()
        .merge(
            locations[["id", "label", lat_col, lon_col]],
            left_on="from",
            right_on="id",
            how="inner",
        )
        .rename(
            columns={
                "label": "label_from",
                lat_col: "Latitude_from",
                lon_col: "Longitude_from",
            }
        )
        .drop(columns=["id"])
        .merge(
            locations[["id", "label", lat_col, lon_col]],
            left_on="to",
            right_on="id",
            how="inner",
        )
        .rename(
            columns={
                "label": "label_to",
                lat_col: "Latitude_to",
                lon_col: "Longitude_to",
            }
        )
        .drop(columns=["id"])
    )

    if graph_df.empty:
        raise ValueError("No edges could be mapped to geographic coordinates.")

    center_lat = graph_df["Latitude_from"].mean()
    center_lon = graph_df["Longitude_from"].mean()

    m = folium.Map(location=[center_lat, center_lon], zoom_start=5, tiles=basemap_layer)

    # Distinct colours – cycle through a palette
    palette = [
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
        "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62",
    ]

    arrow_char = "▶" if arrow_filled else "›"

    for idx, row in graph_df.reset_index(drop=True).iterrows():
        from_lat = float(row["Latitude_from"])
        from_lon = float(row["Longitude_from"])
        to_lat = float(row["Latitude_to"])
        to_lon = float(row["Longitude_to"])
        color = palette[idx % len(palette)]
        tip = f"{row['label_from']} → {row['label_to']}"

        folium.PolyLine(
            locations=[[from_lat, from_lon], [to_lat, to_lon]],
            color=color,
            weight=2,
            opacity=0.85,
            tooltip=tip,
        ).add_to(m)

        if not hide_arrow_head:
            # Place a rotated arrow character at 75 % of the way along the edge
            mid_lat = from_lat + 0.75 * (to_lat - from_lat)
            mid_lon = from_lon + 0.75 * (to_lon - from_lon)
            angle = _bearing(from_lat, from_lon, to_lat, to_lon)
            html = (
                f'<div style="font-size:14px;color:{color};'
                f'transform:rotate({angle:.1f}deg);'
                f'transform-origin:center;">{arrow_char}</div>'
            )
            folium.Marker(
                location=[mid_lat, mid_lon],
                icon=folium.DivIcon(
                    html=html, icon_size=(16, 16), icon_anchor=(8, 8)
                ),
                tooltip=tip,
            ).add_to(m)

    # Collect all unique locations
    all_locs = pd.concat(
        [
            graph_df[["label_from", "Latitude_from", "Longitude_from"]].rename(
                columns={
                    "label_from": "location",
                    "Latitude_from": "latitude",
                    "Longitude_from": "longitude",
                }
            ),
            graph_df[["label_to", "Latitude_to", "Longitude_to"]].rename(
                columns={
                    "label_to": "location",
                    "Latitude_to": "latitude",
                    "Longitude_to": "longitude",
                }
            ),
        ]
    ).drop_duplicates(subset=["location"])

    if show_points:
        for _, row in all_locs.iterrows():
            folium.CircleMarker(
                location=[float(row["latitude"]), float(row["longitude"])],
                color=point_color,
                fill=True,
                fill_color=point_color,
                fill_opacity=point_opacity,
                radius=5,
                weight=1,
                tooltip=row["location"],
            ).add_to(m)

    if show_labels:
        for _, row in all_locs.iterrows():
            html = (
                f'<div style="font-size:13px;color:{label_color};'
                f'font-weight:bold;white-space:nowrap;">'
                f'{row["location"]}</div>'
            )
            folium.Marker(
                location=[float(row["latitude"]), float(row["longitude"])],
                icon=folium.DivIcon(
                    html=html,
                    icon_size=(120, 20),
                    icon_anchor=(0, 0),
                ),
                tooltip=row["location"],
            ).add_to(m)

    return m


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _bearing(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Return the compass bearing (degrees) from point 1 to point 2."""
    lat1r, lat2r = math.radians(lat1), math.radians(lat2)
    dlon = math.radians(lon2 - lon1)
    x = math.sin(dlon) * math.cos(lat2r)
    y = math.cos(lat1r) * math.sin(lat2r) - math.sin(lat1r) * math.cos(lat2r) * math.cos(dlon)
    return math.degrees(math.atan2(x, y))


def _compute_centrality_metrics(G: nx.DiGraph, metastates: list) -> pd.DataFrame:
    """Compute all six centrality metrics and return a metrics DataFrame."""
    sorted_nodes = sorted(G.nodes())

    indegree = dict(G.in_degree())
    outdegree = dict(G.out_degree())
    all_degree = dict(G.degree())

    # Betweenness: unnormalized to match igraph's default
    betweenness = nx.betweenness_centrality(G, normalized=False)

    # Closeness: computed on the undirected view (mode="all" in igraph)
    closeness = nx.closeness_centrality(G.to_undirected())

    rows = []
    for i, node_id in enumerate(sorted_nodes):
        label = metastates[i] if i < len(metastates) else str(node_id)
        od = outdegree.get(node_id, 0)
        ad = all_degree.get(node_id, 0)
        shr = od / ad if ad > 0 else 0.0
        rows.append(
            {
                "Metastates": label,
                "Degree Centrality": ad,
                "Indegree Centrality": indegree.get(node_id, 0),
                "Outdegree Centrality": od,
                "Betweenness Centrality": betweenness.get(node_id, 0.0),
                "Closeness Centrality": closeness.get(node_id, 0.0),
                "Source Hub Ratio": shr,
            }
        )
    return pd.DataFrame(rows)


def _build_pyvis_network(
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
    centrality_metric: int,
    metrics: pd.DataFrame,
) -> Optional["pvnet.Network"]:
    """Return a *pyvis* Network sized by the chosen centrality metric."""
    if not HAS_PYVIS:
        warnings.warn(
            "pyvis is not installed; no interactive visualization will be produced. "
            "Install with:  pip install pyvis",
            stacklevel=3,
        )
        return None

    metric_map = {
        1: ("Indegree Centrality", "Indegree Centrality"),
        2: ("Outdegree Centrality", "Outdegree Centrality"),
        3: ("Betweenness Centrality", "Betweenness Centrality"),
        4: ("Closeness Centrality", "Closeness Centrality"),
        5: ("Degree Centrality", "Degree Centrality"),
        6: ("Source Hub Ratio", "Source Hub Ratio: Sink ~0 / Hub = 0.5 / Source ~1"),
    }
    metric_col, title = metric_map.get(
        centrality_metric, ("Source Hub Ratio", "Transmission Network")
    )

    node_metric = dict(zip(metrics["Metastates"], metrics[metric_col]))

    net = pvnet.Network(
        directed=True,
        notebook=True,
        heading=title,
    )
    net.set_options(
        '{"physics": {"solver": "repulsion"},'
        ' "edges": {"arrows": {"to": {"enabled": true, "scaleFactor": 0.75}}}}'
    )

    for _, row in nodes_df.iterrows():
        label = str(row["label"])
        size = float(node_metric.get(label, 0))
        net.add_node(
            int(row["id"]),
            label=label,
            value=size,
            title=f"{label}: {size:.4f}",
        )

    for _, row in edges_df.iterrows():
        net.add_edge(
            int(row["from"]),
            int(row["to"]),
            value=int(row["value"]),
            title=str(int(row["value"])),
        )

    return net


# ---------------------------------------------------------------------------
# Tree-type processors
# ---------------------------------------------------------------------------


def _process_parsimonious(
    treedata, metadata: pd.DataFrame, column_selection: str
):
    """Process a Newick tree + metadata CSV (parsimonious ancestral states)."""
    if not HAS_BIOPYTHON:
        raise ImportError(
            "Biopython is required for tree_type='parsimonious'. "
            "Install with:  pip install biopython"
        )

    tree = (
        Phylo.read(treedata, "newick")
        if isinstance(treedata, str)
        else treedata
    )

    tip_labels = [c.name for c in tree.get_terminals()]

    sortingtable = pd.DataFrame(
        {"N_ID": range(1, len(tip_labels) + 1), "Accession": tip_labels}
    )
    data = (
        metadata.merge(sortingtable, on="Accession")
        .sort_values("N_ID")
        .reset_index(drop=True)
    )

    char_labels = sorted(data[column_selection].astype(str).unique().tolist())
    state_to_int = {label: i + 1 for i, label in enumerate(char_labels)}
    tip_states = {
        row["Accession"]: state_to_int[str(row[column_selection])]
        for _, row in data.iterrows()
    }

    node_states = fitch_parsimony(tree, tip_states)
    source_list, target_list = extract_state_changes(tree, node_states, char_labels)

    return _build_edges_nodes(source_list, target_list, char_labels)


def _process_bayesian(
    treedata, column_selection: str, threshold: float, threshold2: float
):
    """Process a BEAST MCC nexus tree (Bayesian phylogeography)."""
    try:
        import dendropy
    except ImportError as exc:
        raise ImportError(
            "dendropy is required for tree_type='bayesian'. "
            "Install with:  pip install dendropy"
        ) from exc

    tree = (
        dendropy.Tree.get(path=treedata, schema="nexus")
        if isinstance(treedata, str)
        else treedata
    )

    source_list: list[str] = []
    target_list: list[str] = []

    def _ann(node, name):
        val = node.annotations.get_value(name)
        return None if val is None else str(val).strip('"').strip("'")

    for edge in tree.preorder_edge_iter():
        if edge.tail_node is None or edge.head_node is None:
            continue

        parent = edge.tail_node
        child = edge.head_node

        p_state = _ann(parent, column_selection)
        c_state = _ann(child, column_selection)

        if p_state is None or c_state is None or p_state == c_state:
            continue

        try:
            p_prob = float(_ann(parent, f"{column_selection}.prob") or 1.0)
            c_prob = float(_ann(child, f"{column_selection}.prob") or 1.0)
            posterior = float(_ann(parent, "posterior") or 1.0)
        except (TypeError, ValueError):
            continue

        # Leaf nodes have probability 1 by definition
        if child.is_leaf():
            c_prob = 1.0
        if parent.parent_node is None:
            posterior = 1.0

        if p_prob >= threshold and c_prob >= threshold and posterior >= threshold2:
            source_list.append(p_state)
            target_list.append(c_state)

    all_states = sorted(set(source_list + target_list))
    return _build_edges_nodes(source_list, target_list, all_states)


def _process_nj(
    treedata,
    metadata: pd.DataFrame,
    column_selection: str,
    root_selection: Optional[str],
    bootstrap_value: Optional[int],
):
    """Build a neighbor-joining tree from a FASTA alignment then apply parsimony ASR."""
    if not HAS_BIOPYTHON:
        raise ImportError(
            "Biopython is required for tree_type='nj'. "
            "Install with:  pip install biopython"
        )

    if isinstance(treedata, str):
        alignment = AlignIO.read(treedata, "fasta")
    else:
        alignment = treedata

    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator)
    tree = constructor.nj(calculator.get_distance(alignment))

    if root_selection is not None:
        tree.root_with_outgroup({"name": root_selection})

    if bootstrap_value is not None:
        # Collapse nodes whose bootstrap support is below the threshold
        _collapse_low_support(tree, alignment, calculator, constructor, bootstrap_value)

    tip_labels = [c.name for c in tree.get_terminals()]

    sortingtable = pd.DataFrame(
        {"N_ID": range(1, len(tip_labels) + 1), "Accession": tip_labels}
    )
    data = (
        metadata.merge(sortingtable, on="Accession")
        .sort_values("N_ID")
        .reset_index(drop=True)
    )

    char_labels = sorted(data[column_selection].astype(str).unique().tolist())
    state_to_int = {label: i + 1 for i, label in enumerate(char_labels)}
    tip_states = {
        row["Accession"]: state_to_int[str(row[column_selection])]
        for _, row in data.iterrows()
    }

    node_states = fitch_parsimony(tree, tip_states)
    source_list, target_list = extract_state_changes(tree, node_states, char_labels)

    return _build_edges_nodes(source_list, target_list, char_labels)


def _process_dataframe(treedata: pd.DataFrame):
    """Process a pre-computed edge-list DataFrame."""
    all_nodes = sorted(
        set(treedata["from"].tolist() + treedata["to"].tolist())
    )
    node_to_id = {label: i + 1 for i, label in enumerate(all_nodes)}

    dat = pd.DataFrame(
        {
            "from": [node_to_id[f] for f in treedata["from"]],
            "to": [node_to_id[t] for t in treedata["to"]],
        }
    )
    edges = dat.groupby(["from", "to"]).size().reset_index(name="value")

    # If the input already has edge weights, use them
    if "value" in treedata.columns:
        weight_map = {
            (node_to_id[row["from"]], node_to_id[row["to"]]): row["value"]
            for _, row in treedata.iterrows()
        }
        edges["value"] = [
            weight_map.get((int(r["from"]), int(r["to"])), r["value"])
            for _, r in edges.iterrows()
        ]

    nodes_df = pd.DataFrame(
        {"id": range(1, len(all_nodes) + 1), "label": all_nodes}
    )
    return edges, nodes_df


def _build_edges_nodes(
    source_list: list, target_list: list, metastates: list
):
    """Convert source/target label lists into edges_df + nodes_df with integer IDs."""
    nodes_df = pd.DataFrame(
        {"id": range(1, len(metastates) + 1), "label": metastates}
    )
    label_to_id = {label: i + 1 for i, label in enumerate(metastates)}

    dat = pd.DataFrame({"from": source_list, "to": target_list})
    if dat.empty:
        edges_df = pd.DataFrame(columns=["from", "to", "value"])
    else:
        edges = dat.groupby(["from", "to"]).size().reset_index(name="value")
        edges_df = pd.DataFrame(
            {
                "from": [label_to_id[f] for f in edges["from"]],
                "to": [label_to_id[t] for t in edges["to"]],
                "value": edges["value"].values,
            }
        )
    return edges_df, nodes_df


def _collapse_low_support(tree, alignment, calculator, constructor, threshold: int):
    """
    Bootstrap the NJ tree and collapse branches with support < *threshold*.

    Modifies *tree* in-place by setting branch lengths of poorly supported
    internal branches to zero (equivalent to ``ape::di2multi``).
    """
    from Bio.Phylo.Consensus import bootstrap_trees

    n_bootstrap = 100
    bs_trees = list(bootstrap_trees(alignment, n_bootstrap, constructor.nj))

    # Count how often each internal clade appears in bootstrap trees
    terminal_names = {c.name for c in tree.get_terminals()}

    for clade in tree.get_nonterminals():
        clade_tips = frozenset(c.name for c in clade.get_terminals())
        if clade_tips == terminal_names:
            continue  # skip root
        count = sum(
            1
            for bt in bs_trees
            if any(
                frozenset(c.name for c in nc.get_terminals()) == clade_tips
                for nc in bt.get_nonterminals()
            )
        )
        support = (count / n_bootstrap) * 100
        if support < threshold and clade.branch_length is not None:
            clade.branch_length = 0.0
