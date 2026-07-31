"""
pystrainhub – Python implementation of the StrainHub phylogenetic network toolkit.

Quick-start
-----------
>>> import pandas as pd
>>> from Bio import Phylo
>>> from pystrainhub import list_states, make_transnet, make_map

>>> tree = Phylo.read("tree.phy", "newick")
>>> meta = pd.read_csv("metadata.csv")
>>> geo  = pd.read_csv("geodata.csv")

>>> list_states(tree, meta)
>>> graph = make_transnet(tree, meta, column_selection="Country",
...                       centrality_metric=6, tree_type="parsimonious")
>>> graph["metrics"]
>>> m = make_map(graph, geo, column_selection="Country")
>>> m.save("map.html")
"""

from .strainhub import list_states, make_transnet, make_map

__version__ = "1.0.0"
__author__ = "Colby T. Ford"
__all__ = ["list_states", "make_transnet", "make_map"]
