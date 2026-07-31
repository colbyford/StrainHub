"""
Tests for pystrainhub.

Run with:  pytest pystrainhub/tests/
"""

from __future__ import annotations

import io
import textwrap

import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def simple_newick(tmp_path):
    """Write a tiny 4-tip Newick tree and return the file path."""
    content = "((A:1,B:1):1,(C:1,D:1):1);"
    p = tmp_path / "tree.phy"
    p.write_text(content)
    return str(p)


@pytest.fixture()
def simple_metadata():
    return pd.DataFrame(
        {
            "Accession": ["A", "B", "C", "D"],
            "Country": ["USA", "USA", "UK", "UK"],
            "Host": ["Human", "Human", "Human", "Mouse"],
        }
    )


@pytest.fixture()
def simple_geodata():
    return pd.DataFrame(
        {
            "Country": ["USA", "UK"],
            "Latitude": [38.9, 55.4],
            "Longitude": [-77.0, -3.4],
        }
    )


@pytest.fixture()
def dataframe_input():
    return pd.DataFrame(
        {
            "from": ["thingA", "thingB", "thingB"],
            "to": ["thingB", "thingD", "thingC"],
            "value": [5, 10, 1],
        }
    )


@pytest.fixture()
def simple_fasta(tmp_path):
    """Write a tiny 4-sequence FASTA alignment and return the file path."""
    content = (
        ">A\nATCGATCGATCG\n"
        ">B\nATCGATCGATCG\n"
        ">C\nTTCGATCGATCA\n"
        ">D\nTTCGATCGATCA\n"
    )
    p = tmp_path / "alignment.fasta"
    p.write_text(content)
    return str(p)


@pytest.fixture()
def beast_nexus(tmp_path):
    """Write a minimal BEAST-annotated NEXUS tree and return the file path."""
    content = """\
#NEXUS
Begin taxa;
    Dimensions ntax=4;
    Taxlabels A B C D;
End;

Begin trees;
    Translate
        1 A,
        2 B,
        3 C,
        4 D
    ;
    tree STATE_0 [&R] = (((1[&location="USA",location.prob=0.95,posterior=0.99]:0.1,\
2[&location="USA",location.prob=0.95,posterior=0.99]:0.1)\
[&location="USA",location.prob=0.95,posterior=0.99]:0.2,\
(3[&location="UK",location.prob=0.95,posterior=0.99]:0.1,\
4[&location="UK",location.prob=0.95,posterior=0.99]:0.1)\
[&location="UK",location.prob=0.95,posterior=0.99]:0.2)\
[&location="USA",location.prob=0.95,posterior=1.0]:0.0);
End;
"""
    p = tmp_path / "beast.nexus"
    p.write_text(content)
    return str(p)


# ---------------------------------------------------------------------------
# list_states
# ---------------------------------------------------------------------------


class TestListStates:
    def test_parsimonious_returns_columns(self, simple_newick, simple_metadata):
        from pystrainhub import list_states
        from Bio import Phylo

        tree = Phylo.read(simple_newick, "newick")
        result = list_states(tree, simple_metadata, tree_type="parsimonious")

        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ["Index", "Column"]
        assert list(result["Column"]) == list(simple_metadata.columns)
        assert list(result["Index"]) == list(range(1, len(simple_metadata.columns) + 1))

    def test_nj_returns_columns(self, simple_metadata):
        from pystrainhub import list_states

        result = list_states(None, simple_metadata, tree_type="nj")
        assert list(result["Column"]) == list(simple_metadata.columns)

    def test_missing_metadata_raises(self, simple_newick):
        from pystrainhub import list_states
        from Bio import Phylo

        tree = Phylo.read(simple_newick, "newick")
        with pytest.raises(ValueError, match="metadata is required"):
            list_states(tree, None, tree_type="parsimonious")

    def test_bayesian_returns_annotation_columns(self, beast_nexus):
        pytest.importorskip("dendropy")
        from pystrainhub import list_states

        result = list_states(beast_nexus, tree_type="bayesian")

        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ["Index", "Column"]
        # "location" should appear; "location.prob" should be filtered out
        assert "location" in result["Column"].tolist()
        assert not any(c.endswith(".prob") for c in result["Column"])

    def test_unknown_tree_type_raises(self, simple_metadata):
        from pystrainhub import list_states

        with pytest.raises(ValueError, match="Unknown tree_type"):
            list_states(None, simple_metadata, tree_type="bogus")


# ---------------------------------------------------------------------------
# make_transnet – parsimonious
# ---------------------------------------------------------------------------


class TestMakeTransnetParsimonious:
    def test_returns_dict_with_expected_keys(self, simple_newick, simple_metadata):
        from pystrainhub import make_transnet
        from Bio import Phylo

        tree = Phylo.read(simple_newick, "newick")
        result = make_transnet(
            tree,
            simple_metadata,
            column_selection="Country",
            centrality_metric=6,
            tree_type="parsimonious",
            metrics_output_file="",
        )

        assert isinstance(result, dict)
        for key in ("nodes", "edges", "metrics", "graph"):
            assert key in result, f"Missing key: {key}"

    def test_nodes_and_edges_are_dataframes(self, simple_newick, simple_metadata):
        from pystrainhub import make_transnet
        from Bio import Phylo

        tree = Phylo.read(simple_newick, "newick")
        result = make_transnet(
            tree,
            simple_metadata,
            column_selection="Country",
            centrality_metric=1,
            tree_type="parsimonious",
            metrics_output_file="",
        )
        assert isinstance(result["nodes"], pd.DataFrame)
        assert isinstance(result["edges"], pd.DataFrame)
        assert "id" in result["nodes"].columns
        assert "label" in result["nodes"].columns

    def test_metrics_contain_all_six_metrics(self, simple_newick, simple_metadata):
        from pystrainhub import make_transnet
        from Bio import Phylo

        tree = Phylo.read(simple_newick, "newick")
        result = make_transnet(
            tree,
            simple_metadata,
            column_selection="Country",
            centrality_metric=0,
            tree_type="parsimonious",
            metrics_output_file="",
        )
        expected_cols = {
            "Metastates",
            "Degree Centrality",
            "Indegree Centrality",
            "Outdegree Centrality",
            "Betweenness Centrality",
            "Closeness Centrality",
            "Source Hub Ratio",
        }
        assert expected_cols.issubset(set(result["metrics"].columns))

    def test_as_json_returns_string(self, simple_newick, simple_metadata):
        from pystrainhub import make_transnet
        from Bio import Phylo
        import json

        tree = Phylo.read(simple_newick, "newick")
        result = make_transnet(
            tree,
            simple_metadata,
            column_selection="Country",
            centrality_metric=1,
            tree_type="parsimonious",
            metrics_output_file="",
            as_json=True,
        )
        assert isinstance(result, str)
        parsed = json.loads(result)
        assert "nodes" in parsed
        assert "edges" in parsed
        assert "metrics" in parsed

    def test_metrics_written_to_file(self, simple_newick, simple_metadata, tmp_path):
        from pystrainhub import make_transnet
        from Bio import Phylo

        outfile = str(tmp_path / "metrics.csv")
        tree = Phylo.read(simple_newick, "newick")
        make_transnet(
            tree,
            simple_metadata,
            column_selection="Country",
            centrality_metric=1,
            tree_type="parsimonious",
            metrics_output_file=outfile,
        )
        written = pd.read_csv(outfile)
        assert "Metastates" in written.columns

    def test_state_labels_match_sorted_unique_values(
        self, simple_newick, simple_metadata
    ):
        from pystrainhub import make_transnet
        from Bio import Phylo

        tree = Phylo.read(simple_newick, "newick")
        result = make_transnet(
            tree,
            simple_metadata,
            column_selection="Country",
            centrality_metric=1,
            tree_type="parsimonious",
            metrics_output_file="",
        )
        labels = sorted(result["nodes"]["label"].tolist())
        expected = sorted(simple_metadata["Country"].unique().tolist())
        assert labels == expected


# ---------------------------------------------------------------------------
# make_transnet – dataframe
# ---------------------------------------------------------------------------


class TestMakeTransnetDataframe:
    def test_dataframe_type(self, dataframe_input):
        from pystrainhub import make_transnet

        result = make_transnet(
            dataframe_input,
            tree_type="dataframe",
            centrality_metric=6,
            metrics_output_file="",
        )
        assert isinstance(result, dict)
        assert not result["nodes"].empty
        assert not result["edges"].empty

    def test_dataframe_value_column_used(self, dataframe_input):
        from pystrainhub import make_transnet

        result = make_transnet(
            dataframe_input,
            tree_type="dataframe",
            centrality_metric=1,
            metrics_output_file="",
        )
        # The provided value column should appear in edges
        assert "value" in result["edges"].columns
        total_weight = result["edges"]["value"].sum()
        assert total_weight > 0

    def test_node_labels_sorted(self, dataframe_input):
        from pystrainhub import make_transnet

        result = make_transnet(
            dataframe_input,
            tree_type="dataframe",
            centrality_metric=1,
            metrics_output_file="",
        )
        labels = result["nodes"]["label"].tolist()
        assert labels == sorted(labels)


# ---------------------------------------------------------------------------
# make_transnet – nj
# ---------------------------------------------------------------------------


class TestMakeTransnetNJ:
    def test_returns_dict_with_expected_keys(self, simple_fasta, simple_metadata):
        from pystrainhub import make_transnet

        result = make_transnet(
            simple_fasta,
            simple_metadata,
            column_selection="Country",
            centrality_metric=6,
            tree_type="nj",
            metrics_output_file="",
        )

        assert isinstance(result, dict)
        for key in ("nodes", "edges", "metrics", "graph"):
            assert key in result, f"Missing key: {key}"

    def test_nodes_and_edges_are_dataframes(self, simple_fasta, simple_metadata):
        from pystrainhub import make_transnet

        result = make_transnet(
            simple_fasta,
            simple_metadata,
            column_selection="Country",
            centrality_metric=1,
            tree_type="nj",
            metrics_output_file="",
        )
        assert isinstance(result["nodes"], pd.DataFrame)
        assert isinstance(result["edges"], pd.DataFrame)
        assert "id" in result["nodes"].columns
        assert "label" in result["nodes"].columns

    def test_state_labels_from_metadata(self, simple_fasta, simple_metadata):
        from pystrainhub import make_transnet

        result = make_transnet(
            simple_fasta,
            simple_metadata,
            column_selection="Country",
            centrality_metric=1,
            tree_type="nj",
            metrics_output_file="",
        )
        labels = sorted(result["nodes"]["label"].tolist())
        expected = sorted(simple_metadata["Country"].unique().tolist())
        assert labels == expected

    def test_metrics_contain_all_six_metrics(self, simple_fasta, simple_metadata):
        from pystrainhub import make_transnet

        result = make_transnet(
            simple_fasta,
            simple_metadata,
            column_selection="Country",
            centrality_metric=0,
            tree_type="nj",
            metrics_output_file="",
        )
        expected_cols = {
            "Metastates",
            "Degree Centrality",
            "Indegree Centrality",
            "Outdegree Centrality",
            "Betweenness Centrality",
            "Closeness Centrality",
            "Source Hub Ratio",
        }
        assert expected_cols.issubset(set(result["metrics"].columns))

    def test_as_json_returns_string(self, simple_fasta, simple_metadata):
        import json
        from pystrainhub import make_transnet

        result = make_transnet(
            simple_fasta,
            simple_metadata,
            column_selection="Country",
            centrality_metric=1,
            tree_type="nj",
            metrics_output_file="",
            as_json=True,
        )
        assert isinstance(result, str)
        parsed = json.loads(result)
        assert "nodes" in parsed
        assert "edges" in parsed
        assert "metrics" in parsed

    def test_root_selection(self, simple_fasta, simple_metadata):
        from pystrainhub import make_transnet

        # Should not raise when root_selection is provided
        result = make_transnet(
            simple_fasta,
            simple_metadata,
            column_selection="Country",
            centrality_metric=1,
            tree_type="nj",
            root_selection="A",
            metrics_output_file="",
        )
        assert isinstance(result, dict)
        assert not result["nodes"].empty


# ---------------------------------------------------------------------------
# make_transnet – bayesian
# ---------------------------------------------------------------------------


class TestMakeTransnetBayesian:
    def test_returns_dict_with_expected_keys(self, beast_nexus):
        pytest.importorskip("dendropy")
        from pystrainhub import make_transnet

        result = make_transnet(
            beast_nexus,
            column_selection="location",
            centrality_metric=6,
            tree_type="bayesian",
            metrics_output_file="",
        )

        assert isinstance(result, dict)
        for key in ("nodes", "edges", "metrics", "graph"):
            assert key in result, f"Missing key: {key}"

    def test_nodes_and_edges_are_dataframes(self, beast_nexus):
        pytest.importorskip("dendropy")
        from pystrainhub import make_transnet

        result = make_transnet(
            beast_nexus,
            column_selection="location",
            centrality_metric=1,
            tree_type="bayesian",
            metrics_output_file="",
        )
        assert isinstance(result["nodes"], pd.DataFrame)
        assert isinstance(result["edges"], pd.DataFrame)
        assert "id" in result["nodes"].columns
        assert "label" in result["nodes"].columns

    def test_state_labels_from_annotations(self, beast_nexus):
        pytest.importorskip("dendropy")
        from pystrainhub import make_transnet

        result = make_transnet(
            beast_nexus,
            column_selection="location",
            centrality_metric=1,
            tree_type="bayesian",
            metrics_output_file="",
        )
        labels = sorted(result["nodes"]["label"].tolist())
        assert "UK" in labels
        assert "USA" in labels

    def test_metrics_contain_all_six_metrics(self, beast_nexus):
        pytest.importorskip("dendropy")
        from pystrainhub import make_transnet

        result = make_transnet(
            beast_nexus,
            column_selection="location",
            centrality_metric=0,
            tree_type="bayesian",
            metrics_output_file="",
        )
        expected_cols = {
            "Metastates",
            "Degree Centrality",
            "Indegree Centrality",
            "Outdegree Centrality",
            "Betweenness Centrality",
            "Closeness Centrality",
            "Source Hub Ratio",
        }
        assert expected_cols.issubset(set(result["metrics"].columns))

    def test_threshold_filters_low_probability(self, beast_nexus):
        pytest.importorskip("dendropy")
        from pystrainhub import make_transnet

        # With a very high threshold, low-confidence edges should be excluded
        result_strict = make_transnet(
            beast_nexus,
            column_selection="location",
            centrality_metric=1,
            tree_type="bayesian",
            threshold=0.99,
            threshold2=0.99,
            metrics_output_file="",
        )
        result_lenient = make_transnet(
            beast_nexus,
            column_selection="location",
            centrality_metric=1,
            tree_type="bayesian",
            threshold=0.5,
            threshold2=0.5,
            metrics_output_file="",
        )
        # Lenient threshold should produce at least as many edges as strict
        assert len(result_lenient["edges"]) >= len(result_strict["edges"])

    def test_as_json_returns_string(self, beast_nexus):
        import json
        pytest.importorskip("dendropy")
        from pystrainhub import make_transnet

        result = make_transnet(
            beast_nexus,
            column_selection="location",
            centrality_metric=1,
            tree_type="bayesian",
            metrics_output_file="",
            as_json=True,
        )
        assert isinstance(result, str)
        parsed = json.loads(result)
        assert "nodes" in parsed
        assert "edges" in parsed
        assert "metrics" in parsed


# ---------------------------------------------------------------------------
# Parsimony internals
# ---------------------------------------------------------------------------


class TestFitchParsimony:
    def _make_simple_tree(self):
        """Return a Bio.Phylo tree: ((A,B),(C,D))"""
        from Bio import Phylo

        newick = "((A:1,B:1):1,(C:1,D:1):1);"
        return Phylo.read(io.StringIO(newick), "newick")

    def test_no_state_change_same_states(self):
        from pystrainhub.parsimony import fitch_parsimony, extract_state_changes

        tree = self._make_simple_tree()
        tip_states = {"A": 1, "B": 1, "C": 1, "D": 1}
        node_states = fitch_parsimony(tree, tip_states)
        sources, targets = extract_state_changes(tree, node_states, ["X"])
        assert sources == []
        assert targets == []

    def test_state_change_detected(self):
        from pystrainhub.parsimony import fitch_parsimony, extract_state_changes

        tree = self._make_simple_tree()
        # A,B in state 1 ("USA"); C,D in state 2 ("UK")
        tip_states = {"A": 1, "B": 1, "C": 2, "D": 2}
        node_states = fitch_parsimony(tree, tip_states)
        sources, targets = extract_state_changes(tree, node_states, ["USA", "UK"])
        # There must be at least one transition between the two clades
        assert len(sources) > 0
        assert set(sources + targets) <= {"USA", "UK"}

    def test_all_leaves_assigned(self):
        from pystrainhub.parsimony import fitch_parsimony

        tree = self._make_simple_tree()
        tip_states = {"A": 1, "B": 2, "C": 1, "D": 2}
        node_states = fitch_parsimony(tree, tip_states)
        for clade in tree.get_terminals():
            assert clade._sh_id in node_states, f"Leaf {clade.name} not in node_states"

    def test_internal_nodes_assigned(self):
        from pystrainhub.parsimony import fitch_parsimony

        tree = self._make_simple_tree()
        tip_states = {"A": 1, "B": 1, "C": 2, "D": 2}
        node_states = fitch_parsimony(tree, tip_states)
        for clade in tree.get_nonterminals():
            assert clade._sh_id in node_states, "Internal node not in node_states"


# ---------------------------------------------------------------------------
# make_map
# ---------------------------------------------------------------------------


class TestMakeMap:
    def test_returns_folium_map(self, simple_newick, simple_metadata, simple_geodata):
        pytest.importorskip("folium")
        from pystrainhub import make_transnet, make_map
        from Bio import Phylo

        tree = Phylo.read(simple_newick, "newick")
        graph = make_transnet(
            tree,
            simple_metadata,
            column_selection="Country",
            centrality_metric=6,
            tree_type="parsimonious",
            metrics_output_file="",
        )
        m = make_map(graph, simple_geodata, column_selection="Country")
        import folium

        assert isinstance(m, folium.Map)

    def test_show_labels_and_points(
        self, simple_newick, simple_metadata, simple_geodata
    ):
        pytest.importorskip("folium")
        from pystrainhub import make_transnet, make_map
        from Bio import Phylo

        tree = Phylo.read(simple_newick, "newick")
        graph = make_transnet(
            tree,
            simple_metadata,
            column_selection="Country",
            centrality_metric=6,
            tree_type="parsimonious",
            metrics_output_file="",
        )
        # Should not raise
        make_map(
            graph,
            simple_geodata,
            column_selection="Country",
            show_labels=True,
            show_points=True,
        )

    def test_raises_without_folium(
        self, simple_newick, simple_metadata, simple_geodata, monkeypatch
    ):
        from pystrainhub import make_transnet
        from Bio import Phylo
        import pystrainhub.strainhub as sh_mod

        tree = Phylo.read(simple_newick, "newick")
        graph = make_transnet(
            tree,
            simple_metadata,
            column_selection="Country",
            centrality_metric=6,
            tree_type="parsimonious",
            metrics_output_file="",
        )

        monkeypatch.setattr(sh_mod, "HAS_FOLIUM", False)
        with pytest.raises(ImportError, match="folium"):
            from pystrainhub.strainhub import make_map as _make_map

            _make_map(graph, simple_geodata, column_selection="Country")
