"""
Fitch parsimony ancestral state reconstruction.

Provides a pure-Python implementation of the Fitch parsimony algorithm
for inferring ancestral character states on a phylogenetic tree, equivalent
to the ``castor::asr_max_parsimony`` function used in the R strainhub package.
"""


def _assign_ids(clade, counter=None):
    """Recursively assign a unique integer ``_sh_id`` to every clade."""
    if counter is None:
        counter = [0]
    clade._sh_id = counter[0]
    counter[0] += 1
    for child in clade:
        _assign_ids(child, counter)


def fitch_parsimony(tree, tip_states):
    """
    Fitch parsimony ancestral state reconstruction.

    Assigns the most-parsimonious integer state to every internal node of
    *tree* given the observed integer states at the leaves (*tip_states*).
    The algorithm performs two passes over the tree:

    1. **Bottom-up**: computes a set of candidate states for each node.
       The candidate set of an internal node is the *intersection* of its
       children's candidate sets when non-empty, and the *union* otherwise
       (Fitch 1971).
    2. **Top-down**: picks a single state for each node.  If the parent's
       chosen state is in the node's candidate set, the parent's state is
       inherited; otherwise the minimum candidate state is chosen
       (deterministic tie-breaking).

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        A rooted phylogenetic tree (e.g. from ``Bio.Phylo.read``).
    tip_states : dict
        Mapping of leaf/tip name → integer state (1-indexed).  Tips absent
        from the dict are treated as missing data (all states equally likely).

    Returns
    -------
    dict
        Mapping of ``clade._sh_id`` → integer state (1-indexed) for every
        node in the tree (both leaves and internal nodes).
    """
    _assign_ids(tree.root)

    all_possible_states = frozenset(tip_states.values())
    state_sets = {}   # _sh_id -> frozenset of candidate states (bottom-up)
    node_states = {}  # _sh_id -> single chosen state (top-down)

    # ------------------------------------------------------------------
    # Bottom-up pass
    # ------------------------------------------------------------------
    def _bottom_up(clade):
        if clade.is_terminal():
            state = tip_states.get(clade.name)
            state_sets[clade._sh_id] = (
                frozenset({state}) if state is not None else all_possible_states
            )
            return state_sets[clade._sh_id]

        child_sets = [_bottom_up(child) for child in clade]

        intersection = child_sets[0]
        for cs in child_sets[1:]:
            intersection = intersection & cs

        if intersection:
            state_sets[clade._sh_id] = intersection
        else:
            union = frozenset()
            for cs in child_sets:
                union = union | cs
            state_sets[clade._sh_id] = union

        return state_sets[clade._sh_id]

    # ------------------------------------------------------------------
    # Top-down pass
    # ------------------------------------------------------------------
    def _top_down(clade, parent_state=None):
        candidate = state_sets.get(clade._sh_id, frozenset())

        if not candidate:
            return

        if clade.is_terminal():
            state = tip_states.get(clade.name)
            node_states[clade._sh_id] = (
                state if state is not None
                else (parent_state if parent_state in candidate else min(candidate))
            )
            return

        chosen = (
            parent_state
            if (parent_state is not None and parent_state in candidate)
            else min(candidate)
        )
        node_states[clade._sh_id] = chosen

        for child in clade:
            _top_down(child, chosen)

    _bottom_up(tree.root)
    _top_down(tree.root)

    return node_states


def extract_state_changes(tree, node_states, state_labels):
    """
    Walk all edges of *tree* and collect state changes.

    For each edge (parent → child) where the assigned states differ,
    the parent's label is appended to *source_list* and the child's label
    to *target_list*.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        The same tree used with :func:`fitch_parsimony`.
    node_states : dict
        Output of :func:`fitch_parsimony`.
    state_labels : list of str
        Sorted list of unique state labels.  State integer *k* (1-indexed)
        corresponds to ``state_labels[k - 1]``.

    Returns
    -------
    tuple[list[str], list[str]]
        ``(source_labels, target_labels)`` – parallel lists of state-change
        label strings.
    """
    source_list = []
    target_list = []

    def _walk(parent_clade):
        parent_state = node_states.get(parent_clade._sh_id)
        if parent_state is None:
            return
        for child in parent_clade:
            child_state = node_states.get(child._sh_id)
            if child_state is not None and parent_state != child_state:
                source_list.append(state_labels[parent_state - 1])
                target_list.append(state_labels[child_state - 1])
            _walk(child)

    _walk(tree.root)
    return source_list, target_list
