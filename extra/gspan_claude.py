"""=====================================================================================================================
gspan_claude.py

this code produced by claude sonnet 4.6 -- NOT TESTED

for mining frequent subgraphs from multiple graphs
results = gSpan(
    graphs,           # list of {'vertices': {id: label}, 'edges': [(u,v,elabel)]}
    min_support=3,
    max_vertices=10,  # strongly recommended for large graphs
)

Michael Gribskov 4/3/2026
====================================================================================================================="""
"""
gSpan: Graph-Based Substructure Pattern Mining
Optimized for graphs with UNLABELED vertices (edge labels only).
Supports graphs up to ~400 nodes.

Graph input format:
    {
        'vertices': [0, 1, 2, ...],          # just a list/set of vertex IDs
        'edges':    [(u, v, edge_label), ...]
    }
  OR simply omit 'vertices' and let edges define the graph.

Usage:
    from gspan_unlabeled import gSpan

    graphs = [
        {'edges': [(0,1,'a'), (1,2,'b'), (2,0,'c')]},
        {'edges': [(0,1,'a'), (1,2,'a')]},
        ...
    ]
    results = gSpan(graphs, min_support=2)
    for dfs_code, support in results:
        print(support, dfs_code)
"""

from __future__ import annotations
from collections import defaultdict
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

# ── Constants ────────────────────────────────────────────────────────────────
_VL = 1          # single vertex label used everywhere
_NONE = object() # sentinel

# DFS edge: (from_dfs, to_dfs, from_vlabel, edge_label, to_vlabel)
# Since vertex labels are all _VL, canonical comparisons reduce to edge labels only.
Edge = Tuple[int, int, int, Any, int]
DFSCode = Tuple[Edge, ...]


# ── Graph indexing ───────────────────────────────────────────────────────────

def _index(graph: dict) -> dict:
    """Build adjacency list; assign uniform vertex label."""
    adj: Dict[int, List[Tuple[int, Any]]] = defaultdict(list)
    vertices: Set[int] = set()
    for u, v, el in graph['edges']:
        adj[u].append((v, el))
        adj[v].append((u, el))
        vertices.add(u)
        vertices.add(v)
    # sort neighbour lists for determinism
    return {
        'adj': {v: sorted(nbrs) for v, nbrs in adj.items()},
        'vertices': vertices,
    }


# ── DFS code helpers ─────────────────────────────────────────────────────────

def _is_fwd(e: Edge) -> bool:
    return e[1] > e[0]


def _rightmost_path(dfs_code: DFSCode) -> List[int]:
    if not dfs_code:
        return []
    rm = max(e[1] for e in dfs_code if _is_fwd(e))
    path = [rm]
    cur = rm
    while cur != 0:
        for e in reversed(dfs_code):
            if _is_fwd(e) and e[1] == cur:
                cur = e[0]
                path.append(cur)
                break
    path.reverse()
    return path


# ── Edge sort key (min-DFS ordering) ─────────────────────────────────────────
# With uniform vertex labels the key simplifies to:
#   backward edges < forward edges
#   backward: sort by (from_dfs DESC, to_dfs ASC, el)   -> smaller to_dfs first
#   forward:  sort by (from_dfs ASC, el)

def _ekey(e: Edge):
    frm, to, _, el, _ = e
    if to < frm:           # backward
        return (0, -frm, to, el)
    else:                  # forward
        return (1, frm, el)


# ── Minimum DFS code check ───────────────────────────────────────────────────

def _is_min(dfs_code: DFSCode) -> bool:
    """
    Reconstruct the tiny pattern graph and verify dfs_code equals its min DFS code.
    With unlabeled vertices the reconstruction is simpler and faster.
    """
    if len(dfs_code) <= 1:
        return True  # single edge is always minimal

    # Build adj of the encoded pattern
    adj: Dict[int, List[Tuple[int, Any]]] = defaultdict(list)
    max_v = -1
    for frm, to, _, el, _ in dfs_code:
        adj[frm].append((to, el))
        adj[to].append((frm, el))
        max_v = max(max_v, frm, to)

    min_code: List[Edge] = []

    def grow(
        code: List[Edge],
        mapped: Dict[int, int],   # pattern_vid -> dfs_id
        rmpath: List[int],
        vis: FrozenSet,
    ) -> bool:
        """
        Returns False as soon as we know code > min_code prefix (prune).
        Populates min_code in place.
        """
        step = len(code)
        rm_dfs = rmpath[-1]
        rm_pat = mapped[rm_dfs]
        nxt = max(mapped.values()) + 1
        inv = {v: k for k, v in mapped.items()}

        candidates: List[Edge] = []

        # backward from rightmost
        for path_dfs in rmpath[:-1]:
            path_pat = mapped[path_dfs]
            for nbr_pat, el in adj[rm_pat]:
                if nbr_pat == path_pat:
                    eid = (min(rm_pat, path_pat), max(rm_pat, path_pat))
                    if eid not in vis:
                        candidates.append((rm_dfs, path_dfs, _VL, el, _VL))

        # forward from rightmost-path vertices
        for p_dfs in reversed(rmpath):
            p_pat = mapped[p_dfs]
            for nbr_pat, el in sorted(adj[p_pat]):
                if nbr_pat not in inv:
                    candidates.append((p_dfs, nxt, _VL, el, _VL))

        if not candidates:
            return True

        best_cand = min(candidates, key=_ekey)

        for cand in candidates:
            if cand != best_cand:
                continue

            frm_d, to_d, _, el, _ = cand
            is_fwd = to_d > frm_d

            # Compare against known min_code
            if step < len(min_code):
                if cand > min_code[step]:
                    return False
                if cand < min_code[step]:
                    min_code[step:] = []  # we found a better prefix
            elif step == len(min_code):
                min_code.append(cand)

            code.append(cand)

            if is_fwd:
                frm_pat = mapped[frm_d]
                for nbr_pat, el2 in adj[frm_pat]:
                    if nbr_pat not in inv and el2 == el:
                        new_mapped = dict(mapped)
                        new_mapped[to_d] = nbr_pat
                        eid = (min(frm_pat, nbr_pat), max(frm_pat, nbr_pat))
                        grow(code, new_mapped, rmpath + [to_d], vis | {eid})
                        break
            else:
                frm_pat = mapped[frm_d]
                to_pat  = mapped[to_d]
                eid = (min(frm_pat, to_pat), max(frm_pat, to_pat))
                grow(code, mapped, rmpath, vis | {eid})

            code.pop()

        return True

    # Try all starting vertices (all equivalent since unlabeled)
    for start in range(max_v + 1):
        mapped0 = {0: start}
        grow([], mapped0, [0], frozenset())

    return tuple(min_code) == dfs_code


# ── Projections ───────────────────────────────────────────────────────────────

# Each entry: (graph_id, mapping {dfs_vid: graph_vid}, visited_edges frozenset)
ProjEntry = Tuple[int, Dict[int, int], FrozenSet]
Projection = List[ProjEntry]


def _seed_projections(
    indexed: List[dict],
    min_sup: int,
) -> Dict[DFSCode, Projection]:
    """Enumerate all frequent single-edge patterns."""
    bucket: Dict[tuple, Projection] = defaultdict(list)

    for gid, ig in enumerate(indexed):
        seen: Set[tuple] = set()
        adj = ig['adj']
        for u in adj:
            for v, el in adj[u]:
                if u >= v:
                    continue
                # With no vertex labels the only canonical choice is edge label order.
                # A single edge always has code ((0,1,VL,el,VL),)
                code: DFSCode = ((0, 1, _VL, el, _VL),)
                key = (gid, code)
                if key not in seen:
                    seen.add(key)
                    bucket[code].append(
                        (gid, {0: u, 1: v}, frozenset([(u, v)]))
                    )

    return {k: v for k, v in bucket.items()
            if len({e[0] for e in v}) >= min_sup}


def _extend(
    dfs_code: DFSCode,
    proj: Projection,
    indexed: List[dict],
) -> Dict[Edge, Projection]:
    """Generate all one-edge extensions and their projections."""
    rmpath = _rightmost_path(dfs_code)
    rm_dfs = rmpath[-1]
    nxt_dfs = max(e[1] for e in dfs_code if _is_fwd(e)) + 1
    extensions: Dict[Edge, List[ProjEntry]] = defaultdict(list)

    for gid, mapping, vis in proj:
        adj = indexed[gid]['adj']
        inv: Dict[int, int] = {v: k for k, v in mapping.items()}

        rm_gv = mapping[rm_dfs]

        # Backward edges from rightmost vertex to rightmost-path vertices
        for path_dfs in rmpath[:-1]:
            path_gv = mapping[path_dfs]
            for nbr_gv, el in adj.get(rm_gv, []):
                if nbr_gv == path_gv:
                    eid = (min(rm_gv, path_gv), max(rm_gv, path_gv))
                    if eid not in vis:
                        e: Edge = (rm_dfs, path_dfs, _VL, el, _VL)
                        extensions[e].append((gid, mapping, vis | {eid}))

        # Forward edges from each rightmost-path vertex to new vertices
        for path_dfs in reversed(rmpath):
            path_gv = mapping[path_dfs]
            for nbr_gv, el in adj.get(path_gv, []):
                if nbr_gv not in inv:
                    eid = (min(path_gv, nbr_gv), max(path_gv, nbr_gv))
                    new_map = dict(mapping)
                    new_map[nxt_dfs] = nbr_gv
                    e_fwd: Edge = (path_dfs, nxt_dfs, _VL, el, _VL)
                    extensions[e_fwd].append((gid, new_map, vis | {eid}))

    return extensions


# ── Core recursive miner ──────────────────────────────────────────────────────

def _mine(
    dfs_code: DFSCode,
    proj: Projection,
    indexed: List[dict],
    min_sup: int,
    max_edges: Optional[int],
    results: list,
) -> None:
    support = len({p[0] for p in proj})
    if support < min_sup:
        return
    if not _is_min(dfs_code):
        return
    if dfs_code:
        results.append((dfs_code, support))

    if max_edges is not None and len(dfs_code) >= max_edges:
        return

    exts = _extend(dfs_code, proj, indexed)
    for new_edge in sorted(exts, key=_ekey):
        # deduplicate: keep one entry per graph
        by_graph: Dict[int, ProjEntry] = {}
        for entry in exts[new_edge]:
            gid = entry[0]
            if gid not in by_graph:
                by_graph[gid] = entry
        deduped = list(by_graph.values())
        if len(deduped) >= min_sup:
            _mine(dfs_code + (new_edge,), deduped, indexed,
                  min_sup, max_edges, results)


# ── Public API ────────────────────────────────────────────────────────────────

def gSpan(
    graphs: List[dict],
    min_support: int = 2,
    max_edges: Optional[int] = None,
) -> List[Tuple[DFSCode, int]]:
    """
    Frequent subgraph mining for graphs with unlabeled vertices.

    Parameters
    ----------
    graphs       : list of dicts with key 'edges': [(u, v, edge_label), ...]
    min_support  : minimum number of graphs a pattern must appear in
    max_edges    : optional cap on pattern size (number of edges)

    Returns
    -------
    List of (dfs_code, support) sorted by support descending.
    """
    indexed = [_index(g) for g in graphs]
    results: list = []

    seeds = _seed_projections(indexed, min_support)
    for code1, proj1 in sorted(seeds.items()):
        _mine(code1, proj1, indexed, min_support, max_edges, results)

    results.sort(key=lambda x: -x[1])
    return results


# ── Display ───────────────────────────────────────────────────────────────────

def print_results(results: List[Tuple[DFSCode, int]], limit: int = 20) -> None:
    for i, (code, sup) in enumerate(results[:limit]):
        n_vertices = max(max(e[0], e[1]) for e in code) + 1
        print(f"\nPattern #{i+1}  support={sup}  "
              f"vertices={n_vertices}  edges={len(code)}")
        for step, (frm, to, _, el, _) in enumerate(code):
            arrow = "→" if to > frm else "↩"
            print(f"  [{step}] v{frm} {arrow}[{el}] v{to}")


# ── Smoke test ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import time, random
    random.seed(0)

    LABELS = ['a', 'b', 'c']

    def rand_graph(n_v: int, n_e: int) -> dict:
        verts = list(range(n_v))
        edges_set: set = set()
        # ensure connectivity via a spanning path first
        for i in range(n_v - 1):
            edges_set.add((i, i+1, random.choice(LABELS)))
        while len(edges_set) < n_e:
            u, v = sorted(random.sample(verts, 2))
            edges_set.add((u, v, random.choice(LABELS)))
        return {'edges': list(edges_set)}

    # 50 graphs, up to 20 vertices each
    graphs = [rand_graph(random.randint(4, 20), random.randint(4, 25))
              for _ in range(50)]

    t0 = time.perf_counter()
    results = gSpan(graphs, min_support=8, max_edges=5)
    dt = time.perf_counter() - t0

    print(f"Found {len(results)} frequent subgraphs in {dt:.3f}s")
    print_results(results, limit=10)

# ======================================================================================================================
# canonical labelling
# ======================================================================================================================
"""
Canonical labeling for graphs with unlabeled vertices and edge labels.
Returns the minimum DFS code, which serves as a canonical certificate.
"""

from collections import defaultdict
from typing import Any, Dict, FrozenSet, List, Optional, Tuple

Edge = Tuple[int, int, Any]   # (u, v, edge_label)  — raw graph edge
DFSEdge = Tuple[int, int, Any]  # (from_dfs, to_dfs, edge_label)


def _ekey(e: DFSEdge):
    frm, to, el = e
    if to < frm:               # backward
        return (0, -frm, to, el)
    else:                      # forward
        return (1, frm, el)


def canonical_label(edges: List[Edge]) -> Tuple[DFSEdge, ...]:
    """
    Given a list of (u, v, edge_label) edges, return the minimum DFS code
    as a tuple of (from_dfs, to_dfs, edge_label) — the canonical certificate.

    Two graphs are isomorphic (up to edge-label preserving isomorphism)
    iff their canonical labels are equal.
    """
    # Build adjacency list
    adj: Dict[int, List[Tuple[int, Any]]] = defaultdict(list)
    vertices = set()
    for u, v, el in edges:
        adj[u].append((v, el))
        adj[v].append((u, el))
        vertices.add(u)
        vertices.add(v)

    best: List[DFSEdge] = []

    def grow(
        code: List[DFSEdge],
        mapped: Dict[int, int],   # orig_vid -> dfs_id
        rmpath: List[int],        # dfs_ids on rightmost path
        vis: FrozenSet,           # frozenset of (min,max) orig edge pairs
    ) -> None:
        step = len(code)
        rm_dfs = rmpath[-1]
        rm_orig = mapped[rm_dfs]
        inv = {v: k for k, v in mapped.items()}
        nxt = max(mapped.values()) + 1

        candidates: List[DFSEdge] = []

        # Backward edges from rightmost vertex
        for path_dfs in rmpath[:-1]:
            path_orig = mapped[path_dfs]
            for nbr, el in adj[rm_orig]:
                if nbr == path_orig:
                    eid = (min(rm_orig, path_orig), max(rm_orig, path_orig))
                    if eid not in vis:
                        candidates.append((rm_dfs, path_dfs, el))

        # Forward edges from rightmost-path vertices
        for p_dfs in reversed(rmpath):
            p_orig = mapped[p_dfs]
            for nbr, el in sorted(adj[p_orig]):
                if nbr not in inv:
                    candidates.append((p_dfs, nxt, el))

        if not candidates:
            return

        best_cand = min(candidates, key=_ekey)

        for cand in candidates:
            if cand != best_cand:
                continue

            frm_d, to_d, el = cand
            is_fwd = to_d > frm_d

            # Prune: if this prefix is already worse than best, stop
            if step < len(best):
                if cand > best[step]:
                    return
                if cand < best[step]:
                    del best[step:]
            if step >= len(best):
                best.append(cand)

            code.append(cand)

            if is_fwd:
                frm_orig = mapped[frm_d]
                for nbr, el2 in adj[frm_orig]:
                    if nbr not in inv and el2 == el:
                        new_mapped = dict(mapped)
                        new_mapped[to_d] = nbr
                        eid = (min(frm_orig, nbr), max(frm_orig, nbr))
                        grow(code, new_mapped, rmpath + [to_d], vis | {eid})
                        break
            else:
                frm_orig = mapped[frm_d]
                to_orig  = mapped[to_d]
                eid = (min(frm_orig, to_orig), max(frm_orig, to_orig))
                grow(code, mapped, rmpath, vis | {eid})

            code.pop()

    for start in vertices:
        grow([], {start: 0}, [0], frozenset())

    return tuple(best)


if __name__ == "__main__":
    # Two isomorphic triangles with different vertex numbering
    g1 = [(0, 1, 'a'), (1, 2, 'b'), (2, 0, 'c')]
    g2 = [(3, 7, 'b'), (7, 5, 'c'), (5, 3, 'a')]

    c1 = canonical_label(g1)
    c2 = canonical_label(g2)
    print("g1:", c1)
    print("g2:", c2)
    print("isomorphic:", c1 == c2)
