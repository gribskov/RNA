#!/usr/bin/env python3
"""
Single-file Python + C++ hybrid for canonical labeling of
DIRECTED edge-labeled, vertex-unlabeled graphs using a
gSpan-style canonical DFS code.

Optimized for sparse graphs up to ~400 vertices.
"""

import os
import subprocess
import ctypes
import numpy as np
import tempfile

# ============================================================
# Embedded C++ code (DIRECTED)
# ============================================================

CPP_CODE = r"""
#include <vector>
#include <algorithm>
#include <tuple>
#include <cstdint>
#include <limits>

extern "C" {

struct Edge {
    int to;
    int label;
    int id;
};

struct CodeEdge {
    int from;
    int to;
    int label;
    int dir; // 0 = backward, 1 = forward (DFS sense)
    bool operator<(const CodeEdge& o) const {
        return std::tie(from, to, label, dir)
             < std::tie(o.from, o.to, o.label, o.dir);
    }
};

int canonical_dfs(
    int n_vertices,
    int n_edges,
    const int32_t* edges,
    int32_t* out
) {
    // ---------- build directed adjacency ----------
    std::vector<std::vector<Edge>> adj(n_vertices);
    for (int i = 0; i < n_edges; ++i) {
        int u = edges[3*i];
        int v = edges[3*i+1];
        int l = edges[3*i+2];
        adj[u].push_back({v, l, i});
    }

    // ---------- find smallest edge label ----------
    int min_label = std::numeric_limits<int>::max();
    for (int u = 0; u < n_vertices; ++u)
        for (auto &e : adj[u])
            min_label = std::min(min_label, e.label);

    // ---------- DFS state ----------
    std::vector<int> g2dfs(n_vertices, -1);
    std::vector<int> dfs2g(n_vertices, -1);
    std::vector<char> used_edge(n_edges, 0);
    std::vector<int> rmp;   // right-most path
    std::vector<CodeEdge> code;

    // ---------- choose canonical root edge ----------
    int ru = -1, rv = -1, rid = -1;
    for (int u = 0; u < n_vertices; ++u) {
        for (auto &e : adj[u]) {
            if (e.label != min_label) continue;
            if (ru == -1 || u < ru || (u == ru && e.to < rv)) {
                ru = u;
                rv = e.to;
                rid = e.id;
            }
        }
    }
    if (ru == -1) return 0; // empty graph

    g2dfs[ru] = 0;
    g2dfs[rv] = 1;
    dfs2g[0] = ru;
    dfs2g[1] = rv;
    used_edge[rid] = 1;
    rmp = {0, 1};
    code.push_back({0, 1, min_label, 1});

    // ---------- canonical DFS construction ----------
    while (true) {
        bool extended = false;

        // ---- backward edges ----
        for (int i = (int)rmp.size() - 1; i >= 0 && !extended; --i) {
            int du = rmp[i];
            int u = dfs2g[du];
            for (auto &e : adj[u]) {
                if (used_edge[e.id]) continue;
                int dv = g2dfs[e.to];
                if (dv >= 0 && dv < du) {
                    used_edge[e.id] = 1;
                    code.push_back({du, dv, e.label, 0});
                    extended = true;
                    break;
                }
            }
        }
        if (extended) continue;

        // ---- forward edges ----
        for (int i = (int)rmp.size() - 1; i >= 0 && !extended; --i) {
            int du = rmp[i];
            int u = dfs2g[du];
            for (auto &e : adj[u]) {
                if (used_edge[e.id]) continue;
                if (g2dfs[e.to] == -1) {
                    int dv = (int)rmp.size();
                    g2dfs[e.to] = dv;
                    dfs2g[dv] = e.to;
                    used_edge[e.id] = 1;
                    rmp.push_back(dv);
                    code.push_back({du, dv, e.label, 1});
                    extended = true;
                    break;
                }
            }
        }

        if (!extended) break;
    }

    // ---------- serialize DFS code ----------
    int idx = 0;
    for (auto &e : code) {
        out[idx++] = e.from;
        out[idx++] = e.to;
        out[idx++] = e.label;
        out[idx++] = e.dir;
    }

    return idx;
}

} // extern "C"
"""

# ============================================================
# Build shared library
# ============================================================

def build_shared_library():
    tmpdir = tempfile.gettempdir()
    cpp_path = os.path.join(tmpdir, "canon_dfs_directed.cpp")
    so_path = os.path.join(tmpdir, "libcanon_dfs_directed.so")

    if not os.path.exists(so_path):
        with open(cpp_path, "w") as f:
            f.write(CPP_CODE)
        subprocess.check_call([
            "g++", "-O3", "-std=c++17", "-shared", "-fPIC",
            cpp_path, "-o", so_path
        ])
    return so_path

# ============================================================
# Python binding
# ============================================================

_lib = ctypes.CDLL(build_shared_library())
_lib.canonical_dfs.argtypes = [
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int32),
    ctypes.POINTER(ctypes.c_int32),
]
_lib.canonical_dfs.restype = ctypes.c_int


def canonical_label(n_vertices, edges):
    """
    edges: list of (u, v, label) for DIRECTED edges u -> v
    returns: canonical DFS code as tuple[int]
    """
    n_edges = len(edges)
    edge_arr = np.array(edges, dtype=np.int32).reshape(-1)
    out = np.zeros(4 * n_edges, dtype=np.int32)

    length = _lib.canonical_dfs(
        n_vertices,
        n_edges,
        edge_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
        out.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
    )
    return tuple(out[:length])


# ============================================================
# Example
# ============================================================

if __name__ == "__main__":
    # Directed 3-cycle
    n = 3
    edges = [
        (0, 1, 5),
        (1, 2, 3),
        (2, 0, 5),
    ]

    print("Canonical directed label:")
    print(canonical_label(n, edges))