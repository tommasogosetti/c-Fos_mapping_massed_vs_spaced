#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TOPOLOGY + SMALL-WORLDNESS: stratified Spaced-vs-Massed permutation analysis

Required input files in the SAME folder:
    ss net.csv
    ms net.csv
    sc net.csv
    mc net.csv

Optional:
    hc net.csv
If present, HC is summarized descriptively but is NOT included in the
Spaced-vs-Massed permutation test.

Network reconstruction:
    positive Pearson correlations
    r > 0.65
    conventional TWO-SIDED p < 0.04
    pairwise-complete observations

Observed validation expected from the manuscript networks:
    SS: 374 edges, giant component 121 nodes
    MS: 237 edges, giant component 114 nodes
    SC: 253 edges, giant component 121 nodes
    MC: 430 edges, giant component 113 nodes

Metrics (computed on the giant component):
    - global efficiency
    - characteristic path length
    - average clustering coefficient
    - transitivity
    - small-world coefficient sigma

Small-world coefficient:
    sigma = (C / C_rand) / (L / L_rand)

where C and L are the observed clustering coefficient and characteristic
path length, and C_rand / L_rand are the corresponding means from connected,
degree-preserving randomized networks.

Primary statistical test:
    subject labels are permuted WITHOUT replacement separately within:
       Spatial: SS + MS
       Cue:     SC + MC

For each permutation all four networks are rebuilt from subject-level data.
The primary statistic for each metric is:
    mean[(SS - MS), (SC - MC)]

Two-sided Monte Carlo p values use the +1 correction.

IMPORTANT:
This script does NOT bootstrap subjects with replacement.
"""

from __future__ import annotations

from pathlib import Path
import json
import math
import sys
import time
import traceback

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import t as student_t


# --------------------------------------------------------------------------
# SETTINGS
# --------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "topology_smallworld_results"

FILES = {
    "SS": BASE_DIR / "ss net.csv",
    "MS": BASE_DIR / "ms net.csv",
    "SC": BASE_DIR / "sc net.csv",
    "MC": BASE_DIR / "mc net.csv",
}
HC_FILE = BASE_DIR / "hc net.csv"

R_THRESHOLD = 0.65
P_THRESHOLD = 0.04

EXPECTED_EDGES = {"SS": 374, "MS": 237, "SC": 253, "MC": 430}
EXPECTED_GC_NODES = {"SS": 121, "MS": 114, "SC": 121, "MC": 113}

# --------------------------------------------------------------------------
# PYCHARM CONFIGURATION
# --------------------------------------------------------------------------
# Change ONLY this line before pressing the green Run button:
#   "OBSERVED" = reconstruct/validate the real networks only
#   "TEST"     = quick test with 100 permutations
#   "FINAL"    = final analysis with 5000 permutations
RUN_MODE = "FINAL"

SEED = 20260904
SWAP_MULTIPLIER = 10

MODE_SETTINGS = {
    "OBSERVED": {"permutations": 0,    "random_networks": 50},
    "TEST": {"permutations": 100, "random_networks": 20},
    "FINAL":    {"permutations": 5000, "random_networks": 20},
}


# --------------------------------------------------------------------------
# INPUT
# --------------------------------------------------------------------------

def read_group(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path.name}")

    # The current matrices are comma-separated, but tolerate semicolon files.
    df = pd.read_csv(path)
    if df.shape[1] == 1:
        df = pd.read_csv(path, sep=";")

    df = df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed")]
    df = df.dropna(axis=1, how="all")
    df.columns = [str(c).strip() for c in df.columns]

    subject_candidates = [
        c for c in df.columns
        if c.lower() in {"animale", "animal", "subject", "id"}
    ]
    if not subject_candidates:
        raise ValueError(
            f"{path.name}: no animal-ID column found. Expected 'animale'."
        )
    subject_col = subject_candidates[0]
    df = df.rename(columns={subject_col: "animal"})

    if df["animal"].duplicated().any():
        raise ValueError(f"{path.name}: duplicate animal IDs detected.")

    for c in df.columns:
        if c != "animal":
            df[c] = pd.to_numeric(df[c], errors="coerce")

    if df.columns.duplicated().any():
        raise ValueError(f"{path.name}: duplicate region names detected.")

    return df


def region_matrix(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop(columns=["animal"]).copy()


# --------------------------------------------------------------------------
# FAST PAIRWISE PEARSON, PAIRWISE COMPLETE
# --------------------------------------------------------------------------

def pairwise_pearson_r_p(X: pd.DataFrame):
    arr = X.to_numpy(dtype=float)
    valid = np.isfinite(arr)
    M = valid.astype(float)
    Z = np.where(valid, arr, 0.0)

    n_pair = M.T @ M
    sum_i = Z.T @ M
    sum_j = sum_i.T
    sum_i2 = (Z * Z).T @ M
    sum_j2 = sum_i2.T
    sum_ij = Z.T @ Z

    with np.errstate(divide="ignore", invalid="ignore"):
        cov_num = sum_ij - (sum_i * sum_j / n_pair)
        ss_i = sum_i2 - (sum_i * sum_i / n_pair)
        ss_j = sum_j2 - (sum_j * sum_j / n_pair)
        denom = np.sqrt(ss_i * ss_j)
        r = cov_num / denom

    r = np.clip(r, -1.0, 1.0)
    r[(n_pair < 3) | ~np.isfinite(r)] = np.nan
    np.fill_diagonal(r, 1.0)

    dfree = n_pair - 2.0
    with np.errstate(divide="ignore", invalid="ignore"):
        t_stat = r * np.sqrt(dfree / (1.0 - r * r))
        p_two = 2.0 * student_t.sf(np.abs(t_stat), dfree)

    np.fill_diagonal(p_two, 0.0)
    p_two[(n_pair < 3) | ~np.isfinite(p_two)] = np.nan

    return r, p_two, n_pair


# --------------------------------------------------------------------------
# NETWORK RECONSTRUCTION
# --------------------------------------------------------------------------

def build_network(df: pd.DataFrame, name: str) -> nx.Graph:
    X = region_matrix(df)
    regions = list(X.columns)
    r, p_two, n_pair = pairwise_pearson_r_p(X)

    G = nx.Graph(name=name)
    G.add_nodes_from(regions)

    n_regions = len(regions)
    for i in range(n_regions):
        for j in range(i + 1, n_regions):
            rij = r[i, j]
            pij = p_two[i, j]
            nij = int(n_pair[i, j]) if np.isfinite(n_pair[i, j]) else 0

            if (
                np.isfinite(rij)
                and np.isfinite(pij)
                and rij > R_THRESHOLD
                and pij < P_THRESHOLD
            ):
                G.add_edge(
                    regions[i], regions[j],
                    r=float(rij), p=float(pij), n=nij
                )
    return G


def giant_component(G: nx.Graph) -> nx.Graph:
    if G.number_of_nodes() == 0:
        raise ValueError("Empty graph.")
    comps = list(nx.connected_components(G))
    if not comps:
        raise ValueError("Graph has no connected components.")
    nodes = max(comps, key=len)
    return G.subgraph(nodes).copy()


# --------------------------------------------------------------------------
# TOPOLOGICAL METRICS
# --------------------------------------------------------------------------

def degree_preserving_random_reference(
    H: nx.Graph,
    n_random: int,
    swap_multiplier: int,
    rng: np.random.Generator,
):
    """
    Connected degree-preserving randomization via connected double-edge swaps.
    Returns mean clustering and mean characteristic path length.
    """
    if H.number_of_nodes() < 4 or H.number_of_edges() < 2:
        return np.nan, np.nan, 0

    cvals = []
    lvals = []

    nswap = max(100, swap_multiplier * H.number_of_edges())

    for _ in range(n_random):
        R = H.copy()
        seed = int(rng.integers(0, 2**32 - 1))

        try:
            nx.connected_double_edge_swap(
                R,
                nswap=nswap,
                seed=seed,
            )
            cvals.append(nx.average_clustering(R))
            lvals.append(nx.average_shortest_path_length(R))
        except Exception:
            # One retry with fewer swaps.
            try:
                R = H.copy()
                nx.connected_double_edge_swap(
                    R,
                    nswap=max(20, H.number_of_edges()),
                    seed=seed,
                )
                cvals.append(nx.average_clustering(R))
                lvals.append(nx.average_shortest_path_length(R))
            except Exception:
                continue

    if not cvals:
        return np.nan, np.nan, 0

    return float(np.mean(cvals)), float(np.mean(lvals)), len(cvals)


def calculate_metrics(
    df: pd.DataFrame,
    name: str,
    n_random: int,
    swap_multiplier: int,
    rng: np.random.Generator,
):
    G = build_network(df, name)
    H = giant_component(G)

    eff = nx.global_efficiency(H)
    path = nx.average_shortest_path_length(H)
    clustering = nx.average_clustering(H)
    transitivity = nx.transitivity(H)

    c_rand, l_rand, n_rand_valid = degree_preserving_random_reference(
        H, n_random, swap_multiplier, rng
    )

    if (
        np.isfinite(c_rand) and c_rand > 0
        and np.isfinite(l_rand) and l_rand > 0
        and path > 0
    ):
        sigma = (clustering / c_rand) / (path / l_rand)
    else:
        sigma = np.nan

    return {
        "group": name,
        "n_animals": int(len(df)),
        "n_regions": int(region_matrix(df).shape[1]),
        "missing_values": int(region_matrix(df).isna().sum().sum()),
        "edges": int(G.number_of_edges()),
        "gc_nodes": int(H.number_of_nodes()),
        "gc_edges": int(H.number_of_edges()),
        "density_gc": float(nx.density(H)),
        "global_efficiency": float(eff),
        "characteristic_path_length": float(path),
        "clustering_coefficient": float(clustering),
        "transitivity": float(transitivity),
        "C_random_mean": float(c_rand),
        "L_random_mean": float(l_rand),
        "small_world_sigma": float(sigma),
        "random_networks_valid": int(n_rand_valid),
    }


METRICS = [
    "global_efficiency",
    "characteristic_path_length",
    "clustering_coefficient",
    "transitivity",
    "small_world_sigma",
]


# --------------------------------------------------------------------------
# VALIDATION
# --------------------------------------------------------------------------

def validate_observed(observed: dict[str, dict], stop_on_fail: bool = True):
    rows = []
    exact = True

    for g in ["SS", "MS", "SC", "MC"]:
        obs_edges = observed[g]["edges"]
        obs_gc = observed[g]["gc_nodes"]
        edge_ok = obs_edges == EXPECTED_EDGES[g]
        gc_ok = obs_gc == EXPECTED_GC_NODES[g]

        rows.append({
            "group": g,
            "observed_edges": obs_edges,
            "expected_edges": EXPECTED_EDGES[g],
            "edge_match": edge_ok,
            "observed_gc_nodes": obs_gc,
            "expected_gc_nodes": EXPECTED_GC_NODES[g],
            "gc_match": gc_ok,
        })

        exact = exact and edge_ok and gc_ok

    pd.DataFrame(rows).to_csv(
        OUTPUT_DIR / "network_validation.csv", index=False
    )

    if stop_on_fail and not exact:
        raise RuntimeError(
            "STOP: reconstructed networks do not match the validated "
            "historical edge counts / giant-component sizes. "
            "Check the four input CSVs."
        )

    return exact


# --------------------------------------------------------------------------
# PERMUTATION TEST
# --------------------------------------------------------------------------

def compute_primary_statistics(group_metrics: dict[str, dict]):
    out = {}
    for metric in METRICS:
        ds = group_metrics["SS"][metric] - group_metrics["MS"][metric]
        dc = group_metrics["SC"][metric] - group_metrics["MC"][metric]
        omnibus = (ds + dc) / 2.0
        spaced_mean = (
            group_metrics["SS"][metric] + group_metrics["SC"][metric]
        ) / 2.0
        massed_mean = (
            group_metrics["MS"][metric] + group_metrics["MC"][metric]
        ) / 2.0

        out[metric] = {
            "SS_minus_MS": float(ds),
            "SC_minus_MC": float(dc),
            "mean_spaced": float(spaced_mean),
            "mean_massed": float(massed_mean),
            "mean_spaced_minus_massed": float(omnibus),
        }
    return out


def run_permutations(
    groups: dict[str, pd.DataFrame],
    n_perm: int,
    n_random: int,
    swap_multiplier: int,
    seed: int,
):
    spatial_pool = pd.concat(
        [groups["SS"], groups["MS"]], ignore_index=True
    )
    cue_pool = pd.concat(
        [groups["SC"], groups["MC"]], ignore_index=True
    )

    n_ss = len(groups["SS"])
    n_ms = len(groups["MS"])
    n_sc = len(groups["SC"])
    n_mc = len(groups["MC"])

    spatial_idx = np.arange(len(spatial_pool))
    cue_idx = np.arange(len(cue_pool))

    master_rng = np.random.default_rng(seed)

    records = []
    failures = []

    for b in range(n_perm):
        try:
            ps = master_rng.permutation(spatial_idx)
            pc = master_rng.permutation(cue_idx)

            perm_groups = {
                "SS": spatial_pool.iloc[ps[:n_ss]].reset_index(drop=True),
                "MS": spatial_pool.iloc[ps[n_ss:n_ss+n_ms]].reset_index(drop=True),
                "SC": cue_pool.iloc[pc[:n_sc]].reset_index(drop=True),
                "MC": cue_pool.iloc[pc[n_sc:n_sc+n_mc]].reset_index(drop=True),
            }

            # Deterministic but independent RNG stream for this permutation.
            perm_seed = int(master_rng.integers(0, 2**32 - 1))
            perm_rng = np.random.default_rng(perm_seed)

            gm = {}
            for g in ["SS", "MS", "SC", "MC"]:
                gm[g] = calculate_metrics(
                    perm_groups[g],
                    f"{g}_perm_{b+1}",
                    n_random=n_random,
                    swap_multiplier=swap_multiplier,
                    rng=perm_rng,
                )

            stats = compute_primary_statistics(gm)

            row = {"permutation": b + 1}
            for metric in METRICS:
                row[f"{metric}__SS_minus_MS"] = stats[metric]["SS_minus_MS"]
                row[f"{metric}__SC_minus_MC"] = stats[metric]["SC_minus_MC"]
                row[f"{metric}__mean_spaced_minus_massed"] = (
                    stats[metric]["mean_spaced_minus_massed"]
                )
            records.append(row)

        except Exception as exc:
            failures.append({
                "permutation": b + 1,
                "error": str(exc),
            })

        if (b + 1) % max(1, n_perm // 20) == 0:
            print(f"  permutations: {b+1}/{n_perm}")

    return pd.DataFrame(records), failures


def permutation_pvalues(observed_stats, dist: pd.DataFrame):
    rows = []

    for metric in METRICS:
        observed_delta = observed_stats[metric]["mean_spaced_minus_massed"]
        col = f"{metric}__mean_spaced_minus_massed"
        vals = pd.to_numeric(dist[col], errors="coerce").to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]

        if len(vals) == 0:
            p_two = np.nan
        else:
            p_two = (
                1 + np.sum(np.abs(vals) >= abs(observed_delta))
            ) / (
                1 + len(vals)
            )

        rows.append({
            "metric": metric,
            "SS": observed_metrics_global["SS"][metric],
            "MS": observed_metrics_global["MS"][metric],
            "SC": observed_metrics_global["SC"][metric],
            "MC": observed_metrics_global["MC"][metric],
            "mean_spaced": observed_stats[metric]["mean_spaced"],
            "mean_massed": observed_stats[metric]["mean_massed"],
            "SS_minus_MS": observed_stats[metric]["SS_minus_MS"],
            "SC_minus_MC": observed_stats[metric]["SC_minus_MC"],
            "observed_mean_spaced_minus_massed": observed_delta,
            "n_valid_permutations": int(len(vals)),
            "p_two_sided": float(p_two) if np.isfinite(p_two) else np.nan,
        })

    return pd.DataFrame(rows)


# Global holder used only when building the clean results table.
observed_metrics_global = {}


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------

def main():
    mode = RUN_MODE.strip().upper()
    if mode not in MODE_SETTINGS:
        raise ValueError(
            f"RUN_MODE must be one of {list(MODE_SETTINGS)}, not {RUN_MODE!r}"
        )

    n_perm = MODE_SETTINGS[mode]["permutations"]
    n_random = MODE_SETTINGS[mode]["random_networks"]
    swap_multiplier = SWAP_MULTIPLIER
    seed = SEED

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print(" TOPOLOGY + SMALL-WORLDNESS: SPACED VS MASSED PERMUTATION ANALYSIS")
    print("=" * 72)
    print(f"RUN MODE: {mode}")
    print(f"Script/data folder: {BASE_DIR}")
    print(f"r threshold: > {R_THRESHOLD}")
    print(f"two-sided p threshold: < {P_THRESHOLD}")
    print(f"Permutations: {n_perm}")
    print(f"Random degree-preserving networks per graph: {n_random}")
    print(f"Swap multiplier: {swap_multiplier}")
    print(f"Seed: {seed}")

    print("\nLooking for:")
    for g, p in FILES.items():
        print(f"  {g}: {p.name}")
    if HC_FILE.exists():
        print(f"  HC: {HC_FILE.name} [optional, found]")

    groups = {g: read_group(p) for g, p in FILES.items()}

    print("\nInput:")
    for g in ["SS", "MS", "SC", "MC"]:
        X = region_matrix(groups[g])
        print(
            f"  {g}: n={len(groups[g])}, "
            f"regions={X.shape[1]}, "
            f"missing={int(X.isna().sum().sum())}"
        )

    print("\nObserved network metrics:")
    obs_rng = np.random.default_rng(seed)

    global observed_metrics_global
    observed_metrics_global = {}

    for g in ["SS", "MS", "SC", "MC"]:
        observed_metrics_global[g] = calculate_metrics(
            groups[g],
            g,
            n_random=n_random,
            swap_multiplier=swap_multiplier,
            rng=obs_rng,
        )
        m = observed_metrics_global[g]
        print(
            f"  {g}: edges={m['edges']}, GC={m['gc_nodes']}, "
            f"Eglob={m['global_efficiency']:.6f}, "
            f"L={m['characteristic_path_length']:.6f}, "
            f"C={m['clustering_coefficient']:.6f}, "
            f"T={m['transitivity']:.6f}, "
            f"sigma={m['small_world_sigma']:.6f}"
        )

    if HC_FILE.exists():
        hc = read_group(HC_FILE)
        observed_metrics_global["HC"] = calculate_metrics(
            hc,
            "HC",
            n_random=n_random,
            swap_multiplier=swap_multiplier,
            rng=obs_rng,
        )
        print(
            f"  HC: edges={observed_metrics_global['HC']['edges']}, "
            f"GC={observed_metrics_global['HC']['gc_nodes']}, "
            f"sigma={observed_metrics_global['HC']['small_world_sigma']:.6f}"
        )

    pd.DataFrame(observed_metrics_global.values()).to_csv(
        OUTPUT_DIR / "observed_topology_metrics.csv", index=False
    )

    print("\nValidation against the already-validated manuscript networks:")
    exact = validate_observed(observed_metrics_global, stop_on_fail=True)
    print(f"  exact edge/GC validation = {exact}")

    observed_stats = compute_primary_statistics(observed_metrics_global)

    with open(
        OUTPUT_DIR / "observed_primary_statistics.json",
        "w", encoding="utf-8"
    ) as f:
        json.dump(observed_stats, f, indent=2)

    if mode == "OBSERVED":
        print("\nOBSERVED run completed successfully.")
        print(f"Results folder: {OUTPUT_DIR}")
        print('\nNext: change RUN_MODE = "TEST" and press Run again.')
        return

    print("\nRunning stratified subject-label permutations...")
    print("Spatial labels permuted within SS+MS; cue labels within SC+MC.")
    print("No subject bootstrap with replacement is used.")

    t0 = time.time()
    dist, failures = run_permutations(
        groups,
        n_perm=n_perm,
        n_random=n_random,
        swap_multiplier=swap_multiplier,
        seed=seed,
    )
    elapsed = time.time() - t0

    dist.to_csv(
        OUTPUT_DIR / "permutation_distributions.csv", index=False
    )

    if failures:
        pd.DataFrame(failures).to_csv(
            OUTPUT_DIR / "permutation_failures.csv", index=False
        )

    results = permutation_pvalues(observed_stats, dist)
    results.to_csv(
        OUTPUT_DIR / "spaced_vs_massed_permutation_results.csv",
        index=False
    )

    settings = {
        "run_mode": mode,
        "r_threshold": R_THRESHOLD,
        "p_threshold": P_THRESHOLD,
        "p_sidedness": "two-sided",
        "network_metrics_component": "giant component",
        "small_world_formula": "(C/C_rand)/(L/L_rand)",
        "random_reference": (
            "connected degree-preserving double-edge-swap random networks"
        ),
        "permutations_requested": n_perm,
        "permutations_valid": int(len(dist)),
        "permutations_failed": int(len(failures)),
        "random_networks_per_graph": n_random,
        "swap_multiplier": swap_multiplier,
        "seed": seed,
        "elapsed_seconds": elapsed,
    }
    with open(
        OUTPUT_DIR / "analysis_settings.json",
        "w", encoding="utf-8"
    ) as f:
        json.dump(settings, f, indent=2)

    print("\n" + "=" * 72)
    print("RESULTS")
    print("=" * 72)
    print(results.to_string(index=False))
    print(f"\nElapsed: {elapsed/60:.1f} min")
    print(f"Output folder: {OUTPUT_DIR}")

    if mode == "TEST":
        print('\nTEST completed. If everything is correct, change RUN_MODE = "FINAL".')
    else:
        print("\nFINAL analysis completed.")


if __name__ == "__main__":
    try:
        main()
    except ModuleNotFoundError as exc:
        print("\nMissing Python package:", exc)
        print(
            "\nIn the PyCharm terminal run:\n"
            "pip install numpy pandas scipy networkx"
        )
        raise
    except Exception:
        traceback.print_exc()
        raise
