#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PAIRED SUBJECT-LEVEL BOOTSTRAP:
Is betweenness-targeted disruption more destructive than
Collective Influence, Degree, Katz, and Random attack?

Experimental groups remain separate:
    SS, MS, SC, MC

For every bootstrap iteration and for every group:
1) resample that group's animals WITH replacement;
2) rebuild the c-Fos covariance network from the resampled animals;
3) apply ALL FIVE attack strategies to the SAME reconstructed network;
4) calculate normalized AUC of giant-component decay over q = 0..0.20;
5) calculate paired contrasts:
       AUC_betweenness - AUC_other

Interpretation:
    lower AUC = faster fragmentation / more destructive attack
    negative contrast = betweenness is MORE destructive than the comparator

The primary general contrast is the equal-weight mean across SS, MS, SC, MC.

This is a bootstrap analysis of sampling uncertainty, NOT a permutation of
attack labels. The main inferential output is the bootstrap confidence interval
of each paired AUC contrast. If the entire CI is below zero, that supports
greater destructiveness of betweenness.

Put this script in the same folder as:
    ss net.csv
    ms net.csv
    sc net.csv
    mc net.csv
    ssnet 0.65 0.04.graphml
    msnet 0.65.graphml
    scnet 0.65.graphml
    mcnet_0.65.graphml
"""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
import os
import traceback
import warnings

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import t as student_t


# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "attack_strategy_bootstrap_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

ZSCORE_FILES = {
    "SS": BASE_DIR / "ss net.csv",
    "MS": BASE_DIR / "ms net.csv",
    "SC": BASE_DIR / "sc net.csv",
    "MC": BASE_DIR / "mc net.csv",
}

REFERENCE_GRAPH_FILES = {
    "SS": BASE_DIR / "ssnet 0.65 0.04.graphml",
    "MS": BASE_DIR / "msnet 0.65.graphml",
    "SC": BASE_DIR / "scnet 0.65.graphml",
    "MC": BASE_DIR / "mcnet_0.65.graphml",
}

R_THRESHOLD = 0.65
P_THRESHOLD = 0.04

MIN_COMPONENT_SIZE = 5
N_ATTACK_MEASUREMENTS = 30

AUC_QMAX = 0.20
AUC_GRID_POINTS = 251

KATZ_ALPHA = 0.06
CI_RADIUS = 4

# Start with 1000. For a final paper run, 5000 is preferable if runtime allows.
N_BOOTSTRAPS = 1000

# Random is stochastic, so each bootstrap network gets the MEAN AUC of several
# independent random attacks. 100 is a good compromise.
N_RANDOM_ATTACKS_PER_BOOTSTRAP_NETWORK = 100

RANDOM_SEED = 20260903

# Windows-safe multiprocessing.
# Set N_JOBS = 1 if you prefer a single process or encounter multiprocessing
# issues. 4 is conservative for most desktop CPUs.
N_JOBS = min(4, max(1, (os.cpu_count() or 2) - 1))

# Confidence intervals.
CI_LEVEL = 0.95

# There are four planned primary comparisons:
# betweenness vs random, degree, katz, collective influence.
# We also report a conservative Bonferroni familywise interval.
N_PRIMARY_COMPARISONS = 4
FAMILYWISE_ALPHA = 0.05


# ============================================================================
# INPUT / CLEANING
# ============================================================================

def normalize_region_name(name):
    name = str(name).strip()
    return {
        "NACC SHELL": "NACC.SHELL",
        "NACC CORE": "NACC.CORE",
    }.get(name, name)


def read_zscore_matrix(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path.name}")

    df = pd.read_csv(path, na_values=["NA", "NaN", "nan", ""])
    df = df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed")]
    df = df.dropna(axis=1, how="all")
    df.columns = [normalize_region_name(c) for c in df.columns]

    subject_candidates = [
        c for c in df.columns
        if str(c).lower() in {"animale", "animal", "subject", "id"}
    ]
    if not subject_candidates:
        raise ValueError(f"{path.name}: subject column not found.")

    df = df.rename(columns={subject_candidates[0]: "animal"})

    for c in df.columns:
        if c != "animal":
            df[c] = pd.to_numeric(df[c], errors="coerce")

    return df


def region_matrix(df):
    return df.drop(columns=["animal"]).copy()


# ============================================================================
# EXACT NETWORK CONSTRUCTION
# ============================================================================

def pairwise_pearson_r_p(X):
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
        r = cov_num / np.sqrt(ss_i * ss_j)

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


def build_thresholded_network(df, name="network"):
    X = region_matrix(df)
    regions = list(X.columns)

    r, p_two, n_pair = pairwise_pearson_r_p(X)

    G = nx.Graph(name=name)
    G.add_nodes_from(regions)

    for i in range(len(regions)):
        for j in range(i + 1, len(regions)):
            rij = r[i, j]
            pij = p_two[i, j]
            if (
                np.isfinite(rij)
                and np.isfinite(pij)
                and rij > R_THRESHOLD
                and pij < P_THRESHOLD
            ):
                G.add_edge(
                    regions[i],
                    regions[j],
                    r=float(rij),
                    p=float(pij),
                    n=int(n_pair[i, j]),
                )

    return G


def load_original_graphml(path, name):
    if not path.exists():
        raise FileNotFoundError(f"Missing reference GraphML: {path.name}")

    G = nx.read_graphml(path).to_undirected()

    attrs = nx.get_node_attributes(G, "shared name")
    if attrs:
        nx.relabel_nodes(G, attrs, copy=False)

    mapping = {n: normalize_region_name(n) for n in G.nodes()}
    nx.relabel_nodes(G, mapping, copy=False)

    G.name = name
    return nx.Graph(G)


def edge_set(G):
    return {frozenset((u, v)) for u, v in G.edges() if u != v}


def validate_observed_networks(groups):
    rows = []

    print("\n--- Exact observed-network validation ---")
    for group in ["SS", "MS", "SC", "MC"]:
        recon = build_thresholded_network(groups[group], f"{group}_reconstructed")
        ref = load_original_graphml(REFERENCE_GRAPH_FILES[group], group)

        a = edge_set(ref)
        b = edge_set(recon)

        exact = a == b
        jaccard = len(a & b) / len(a | b) if (a | b) else 1.0

        rows.append({
            "group": group,
            "reference_edges": len(a),
            "reconstructed_edges": len(b),
            "common_edges": len(a & b),
            "edge_jaccard": jaccard,
            "exact_edge_match": exact,
        })

        print(
            f"{group}: ref={len(a)}, reconstructed={len(b)}, "
            f"Jaccard={jaccard:.6f}, EXACT={exact}"
        )

        if not exact:
            raise RuntimeError(
                f"STOP: {group} reconstructed network does not exactly match "
                "the historical GraphML."
            )

    pd.DataFrame(rows).to_csv(
        OUTPUT_DIR / "observed_network_validation.csv",
        index=False
    )


# ============================================================================
# COMMON ATTACK STARTING NETWORK
# ============================================================================

def largest_component_size(G):
    if G.number_of_nodes() == 0:
        return 0
    return max((len(c) for c in nx.connected_components(G)), default=0)


def common_initial_graph(G):
    """
    All attack strategies must begin from the same graph.

    We remove components with fewer than MIN_COMPONENT_SIZE nodes ONCE before
    the attack, matching the common strategy-comparison analysis.
    """
    G = nx.Graph(G.copy())

    comps = list(nx.connected_components(G))
    keep = set().union(*[c for c in comps if len(c) >= MIN_COMPONENT_SIZE]) if comps else set()

    if not keep:
        return nx.Graph()

    return nx.Graph(G.subgraph(keep).copy())


# ============================================================================
# CENTRALITY FUNCTIONS
# ============================================================================

def highest_degree_node(G):
    ranking = sorted(
        nx.degree(G),
        key=lambda x: x[1],
        reverse=True
    )
    return ranking[0][0], float(ranking[0][1])


def highest_betweenness_node(G):
    values = nx.betweenness_centrality(
        G,
        normalized=True,
        weight=None
    )
    ranking = sorted(
        values.items(),
        key=lambda x: x[1],
        reverse=True
    )
    return ranking[0][0], float(ranking[0][1])


def highest_katz_node(G):
    # Preserve the historical alpha used in the original analysis.
    # The direct numpy solver avoids iterative convergence failures.
    values = nx.katz_centrality_numpy(
        G,
        alpha=KATZ_ALPHA,
        beta=1.0,
        normalized=True,
        weight=None
    )
    ranking = sorted(
        values.items(),
        key=lambda x: x[1],
        reverse=True
    )
    return ranking[0][0], float(ranking[0][1])


def collective_influence_values(G, radius=CI_RADIUS):
    degree = dict(G.degree())
    values = {}

    for node in G.nodes():
        distances = nx.single_source_shortest_path_length(
            G,
            node,
            cutoff=radius
        )

        boundary_sum = sum(
            degree[other] - 1
            for other, distance in distances.items()
            if distance == radius
        )

        values[node] = (degree[node] - 1) * boundary_sum

    return values


def highest_ci_node(G):
    values = collective_influence_values(G, CI_RADIUS)

    # Stable sorting preserves graph node order for ties.
    ranking = sorted(
        values.items(),
        key=lambda x: x[1],
        reverse=True
    )
    return ranking[0][0], float(ranking[0][1])


# ============================================================================
# ATTACK CURVES
# ============================================================================

def targeted_attack_curve(G_input, strategy):
    G = common_initial_graph(G_input)

    if G.number_of_nodes() == 0:
        raise ValueError("No nodes remain after initial component filtering.")

    gc0 = largest_component_size(G)
    rows = []

    chooser = {
        "degree": highest_degree_node,
        "katz": highest_katz_node,
        "betweenness": highest_betweenness_node,
        "collective_influence": highest_ci_node,
    }[strategy]

    for step in range(N_ATTACK_MEASUREMENTS):
        gc = largest_component_size(G)

        rows.append({
            "step": step,
            "nodes_removed": step,
            "fraction_gc_nodes_removed": step / gc0,
            "gc_size": gc,
            "relative_gc": gc / gc0,
        })

        if step == N_ATTACK_MEASUREMENTS - 1 or G.number_of_nodes() == 0:
            break

        node, _ = chooser(G)
        G.remove_node(node)

    return pd.DataFrame(rows)


def random_attack_curve(G_input, rng):
    G = common_initial_graph(G_input)

    if G.number_of_nodes() == 0:
        raise ValueError("No nodes remain after initial component filtering.")

    gc0 = largest_component_size(G)
    rows = []

    # One random permutation is equivalent to choosing a random remaining node
    # at each step, and is faster/reproducible.
    removal_order = list(G.nodes())
    rng.shuffle(removal_order)

    for step in range(N_ATTACK_MEASUREMENTS):
        gc = largest_component_size(G)

        rows.append({
            "step": step,
            "nodes_removed": step,
            "fraction_gc_nodes_removed": step / gc0,
            "gc_size": gc,
            "relative_gc": gc / gc0,
        })

        if step == N_ATTACK_MEASUREMENTS - 1 or G.number_of_nodes() == 0:
            break

        G.remove_node(removal_order[step])

    return pd.DataFrame(rows)


def normalized_auc(curve):
    x = curve["fraction_gc_nodes_removed"].to_numpy(dtype=float)
    y = curve["relative_gc"].to_numpy(dtype=float)

    if len(x) < 2 or x[-1] < AUC_QMAX:
        raise ValueError(
            f"Attack curve reaches only q={x[-1] if len(x) else np.nan:.3f}; "
            f"requested q={AUC_QMAX:.3f}."
        )

    grid = np.linspace(0.0, AUC_QMAX, AUC_GRID_POINTS)
    yi = np.interp(grid, x, y)

    if hasattr(np, "trapezoid"):
        area = np.trapezoid(yi, grid)
    else:
        area = np.trapz(yi, grid)

    return float(area / AUC_QMAX)


def calculate_all_strategy_aucs(G, rng, n_random):
    results = {}

    for strategy in [
        "degree",
        "katz",
        "betweenness",
        "collective_influence",
    ]:
        curve = targeted_attack_curve(G, strategy)
        results[strategy] = normalized_auc(curve)

    random_aucs = []
    for _ in range(n_random):
        curve = random_attack_curve(G, rng)
        random_aucs.append(normalized_auc(curve))

    results["random"] = float(np.mean(random_aucs))
    results["random_sd_within_network"] = float(np.std(random_aucs, ddof=1))
    return results


# ============================================================================
# OBSERVED ANALYSIS
# ============================================================================

def observed_strategy_analysis(groups):
    rows = []

    print("\n--- Observed strategy AUCs ---")

    for gi, group in enumerate(["SS", "MS", "SC", "MC"]):
        G = build_thresholded_network(groups[group], f"{group}_observed")
        rng = np.random.default_rng(RANDOM_SEED + 10000 + gi)

        vals = calculate_all_strategy_aucs(
            G,
            rng,
            n_random=5000,  # stable observed random baseline
        )

        for strategy in [
            "random",
            "degree",
            "katz",
            "betweenness",
            "collective_influence",
        ]:
            rows.append({
                "group": group,
                "strategy": strategy,
                "auc": vals[strategy],
            })
            print(f"{group:2s} {strategy:22s} AUC = {vals[strategy]:.6f}")

    observed = pd.DataFrame(rows)
    observed.to_csv(
        OUTPUT_DIR / "observed_strategy_auc_by_group.csv",
        index=False
    )

    wide = observed.pivot(index="group", columns="strategy", values="auc")

    contrast_rows = []
    for comparator in [
        "random",
        "degree",
        "katz",
        "collective_influence",
    ]:
        deltas = wide["betweenness"] - wide[comparator]
        contrast_rows.append({
            "comparison": f"betweenness_minus_{comparator}",
            "observed_mean_delta_across_groups": float(deltas.mean()),
            "SS_delta": float(deltas.loc["SS"]),
            "MS_delta": float(deltas.loc["MS"]),
            "SC_delta": float(deltas.loc["SC"]),
            "MC_delta": float(deltas.loc["MC"]),
        })

    contrasts = pd.DataFrame(contrast_rows)
    contrasts.to_csv(
        OUTPUT_DIR / "observed_betweenness_contrasts.csv",
        index=False
    )

    return observed, contrasts


# ============================================================================
# BOOTSTRAP WORKER
# ============================================================================

WORKER_GROUPS = None


def worker_initializer():
    global WORKER_GROUPS
    WORKER_GROUPS = {
        g: read_zscore_matrix(ZSCORE_FILES[g])
        for g in ["SS", "MS", "SC", "MC"]
    }


def resample_subjects_with_replacement(df, rng):
    n = len(df)
    idx = rng.integers(0, n, size=n)
    return df.iloc[idx].reset_index(drop=True)


def one_bootstrap_iteration(iteration):
    """
    Returns one row per group containing the five AUCs.
    All strategies within a group are applied to the SAME bootstrap network.
    """
    global WORKER_GROUPS

    if WORKER_GROUPS is None:
        worker_initializer()

    seed = RANDOM_SEED + 1_000_003 * (iteration + 1)
    rng = np.random.default_rng(seed)

    rows = []

    for group_index, group in enumerate(["SS", "MS", "SC", "MC"]):
        df = WORKER_GROUPS[group]
        boot_df = resample_subjects_with_replacement(df, rng)

        G = build_thresholded_network(
            boot_df,
            f"{group}_bootstrap_{iteration}"
        )

        # Separate deterministic RNG stream for random attacks in this group.
        attack_rng = np.random.default_rng(
            seed + 10_000 * (group_index + 1)
        )

        vals = calculate_all_strategy_aucs(
            G,
            attack_rng,
            N_RANDOM_ATTACKS_PER_BOOTSTRAP_NETWORK,
        )

        rows.append({
            "bootstrap": iteration + 1,
            "group": group,
            "random": vals["random"],
            "degree": vals["degree"],
            "katz": vals["katz"],
            "betweenness": vals["betweenness"],
            "collective_influence": vals["collective_influence"],
            "random_sd_within_network": vals["random_sd_within_network"],
        })

    return rows


# ============================================================================
# BOOTSTRAP SUMMARIES
# ============================================================================

def percentile_ci(values, level):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    alpha = 1.0 - level
    lo = 100.0 * (alpha / 2.0)
    hi = 100.0 * (1.0 - alpha / 2.0)

    return (
        float(np.percentile(values, lo)),
        float(np.percentile(values, hi)),
    )


def bootstrap_tail_support(values):
    """
    Descriptive bootstrap sign support, NOT an exact permutation p-value.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    prop_negative = float(np.mean(values < 0))
    prop_positive = float(np.mean(values > 0))

    return prop_negative, prop_positive


def summarize_bootstrap(boot_df, observed_contrasts):
    comparison_names = [
        ("random", "betweenness_minus_random"),
        ("degree", "betweenness_minus_degree"),
        ("katz", "betweenness_minus_katz"),
        ("collective_influence", "betweenness_minus_collective_influence"),
    ]

    # Familywise conservative interval: alpha/4 for each comparison.
    per_comparison_alpha = FAMILYWISE_ALPHA / N_PRIMARY_COMPARISONS
    familywise_level = 1.0 - per_comparison_alpha

    general_rows = []
    group_rows = []

    for comparator, comparison in comparison_names:
        # Group-specific paired contrasts
        temp = boot_df.copy()
        temp["delta"] = temp["betweenness"] - temp[comparator]

        for group in ["SS", "MS", "SC", "MC"]:
            vals = temp.loc[temp["group"] == group, "delta"].to_numpy()
            ci95 = percentile_ci(vals, CI_LEVEL)
            prop_neg, prop_pos = bootstrap_tail_support(vals)

            group_rows.append({
                "comparison": comparison,
                "group": group,
                "bootstrap_mean_delta": float(np.nanmean(vals)),
                "bootstrap_sd_delta": float(np.nanstd(vals, ddof=1)),
                "ci95_low": ci95[0],
                "ci95_high": ci95[1],
                "proportion_bootstrap_delta_below_zero": prop_neg,
                "proportion_bootstrap_delta_above_zero": prop_pos,
                "n_valid": int(np.isfinite(vals).sum()),
            })

        # Equal-weight mean across the four experimental networks per bootstrap.
        wide_delta = temp.pivot(
            index="bootstrap",
            columns="group",
            values="delta"
        )

        general_delta = wide_delta[
            ["SS", "MS", "SC", "MC"]
        ].mean(axis=1)

        vals = general_delta.to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]

        ci95 = percentile_ci(vals, CI_LEVEL)
        ci_family = percentile_ci(vals, familywise_level)
        prop_neg, prop_pos = bootstrap_tail_support(vals)

        observed_value = float(
            observed_contrasts.loc[
                observed_contrasts["comparison"] == comparison,
                "observed_mean_delta_across_groups"
            ].iloc[0]
        )

        general_rows.append({
            "comparison": comparison,
            "interpretation": (
                "negative delta = betweenness more destructive"
            ),
            "observed_mean_delta_across_groups": observed_value,
            "bootstrap_mean_delta": float(np.mean(vals)),
            "bootstrap_median_delta": float(np.median(vals)),
            "bootstrap_sd_delta": float(np.std(vals, ddof=1)),
            "ci95_low": ci95[0],
            "ci95_high": ci95[1],
            "familywise_interval_level": familywise_level,
            "familywise_ci_low": ci_family[0],
            "familywise_ci_high": ci_family[1],
            "proportion_bootstrap_delta_below_zero": prop_neg,
            "proportion_bootstrap_delta_above_zero": prop_pos,
            "n_valid_bootstraps": int(len(vals)),
        })

        # Save the actual general bootstrap contrast distribution.
        pd.DataFrame({
            "bootstrap": wide_delta.index,
            "mean_delta_across_SS_MS_SC_MC": general_delta.values,
        }).to_csv(
            OUTPUT_DIR / f"{comparison}_bootstrap_distribution.csv",
            index=False
        )

    general = pd.DataFrame(general_rows)
    groups = pd.DataFrame(group_rows)

    general.to_csv(
        OUTPUT_DIR / "bootstrap_betweenness_general_contrast_summary.csv",
        index=False
    )
    groups.to_csv(
        OUTPUT_DIR / "bootstrap_betweenness_group_contrast_summary.csv",
        index=False
    )

    return general, groups


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("\n============================================================")
    print(" PAIRED SUBJECT-LEVEL BOOTSTRAP")
    print(" Is BETWEenness more destructive than other attack strategies?")
    print("============================================================")

    print(f"\nWorking folder: {BASE_DIR}")
    print(f"Bootstrap iterations: {N_BOOTSTRAPS}")
    print(
        "Random attacks per bootstrap network: "
        f"{N_RANDOM_ATTACKS_PER_BOOTSTRAP_NETWORK}"
    )
    print(f"Parallel workers: {N_JOBS}")
    print(f"Common AUC interval: q=0 to {AUC_QMAX}")
    print(f"Katz alpha: {KATZ_ALPHA}")
    print(f"CI radius: d={CI_RADIUS}")

    groups = {
        g: read_zscore_matrix(ZSCORE_FILES[g])
        for g in ["SS", "MS", "SC", "MC"]
    }

    for g in ["SS", "MS", "SC", "MC"]:
        X = region_matrix(groups[g])
        print(
            f"{g}: n={len(groups[g])}, regions={X.shape[1]}, "
            f"missing={int(X.isna().sum().sum())}"
        )

    # Sanity check: exact original network reconstruction.
    validate_observed_networks(groups)

    # Observed values.
    observed, observed_contrasts = observed_strategy_analysis(groups)

    print("\n--- Starting paired subject bootstrap ---")
    print(
        "Each bootstrap rebuilds SS/MS/SC/MC separately and applies all five "
        "strategies to the same bootstrap network within each group."
    )

    all_rows = []
    failed = []

    if N_JOBS == 1:
        worker_initializer()

        for b in range(N_BOOTSTRAPS):
            try:
                rows = one_bootstrap_iteration(b)
                all_rows.extend(rows)
            except Exception as exc:
                failed.append({
                    "bootstrap": b + 1,
                    "error": str(exc),
                })

            if (b + 1) % max(1, N_BOOTSTRAPS // 20) == 0:
                print(f"  {b+1}/{N_BOOTSTRAPS} bootstraps completed")

    else:
        with ProcessPoolExecutor(
            max_workers=N_JOBS,
            initializer=worker_initializer,
        ) as pool:

            futures = {
                pool.submit(one_bootstrap_iteration, b): b
                for b in range(N_BOOTSTRAPS)
            }

            completed = 0

            for future in as_completed(futures):
                b = futures[future]
                completed += 1

                try:
                    rows = future.result()
                    all_rows.extend(rows)
                except Exception as exc:
                    failed.append({
                        "bootstrap": b + 1,
                        "error": str(exc),
                    })

                if completed % max(1, N_BOOTSTRAPS // 20) == 0:
                    print(
                        f"  {completed}/{N_BOOTSTRAPS} bootstraps completed"
                    )

    if not all_rows:
        raise RuntimeError("No valid bootstrap iterations were produced.")

    boot_df = pd.DataFrame(all_rows).sort_values(
        ["bootstrap", "group"]
    )

    boot_df.to_csv(
        OUTPUT_DIR / "bootstrap_strategy_auc_by_group.csv",
        index=False
    )

    if failed:
        pd.DataFrame(failed).to_csv(
            OUTPUT_DIR / "bootstrap_failures.csv",
            index=False
        )

    valid_bootstraps = boot_df["bootstrap"].nunique()

    print(
        f"\nValid bootstrap iterations represented: "
        f"{valid_bootstraps}/{N_BOOTSTRAPS}"
    )
    print(f"Failed iterations: {len(failed)}")

    general, group_summary = summarize_bootstrap(
        boot_df,
        observed_contrasts
    )

    print("\n============================================================")
    print(" GENERAL BETWEEnNESS CONTRASTS")
    print(" Negative delta = BETWEEnNESS MORE DESTRUCTIVE")
    print("============================================================")

    for _, row in general.iterrows():
        print(f"\n{row['comparison']}")
        print(
            f"  observed mean delta = "
            f"{row['observed_mean_delta_across_groups']:.6f}"
        )
        print(
            f"  bootstrap mean delta = "
            f"{row['bootstrap_mean_delta']:.6f}"
        )
        print(
            f"  95% bootstrap CI = "
            f"[{row['ci95_low']:.6f}, {row['ci95_high']:.6f}]"
        )
        print(
            f"  Bonferroni familywise CI "
            f"({100*row['familywise_interval_level']:.2f}%) = "
            f"[{row['familywise_ci_low']:.6f}, "
            f"{row['familywise_ci_high']:.6f}]"
        )
        print(
            f"  bootstrap proportion delta < 0 = "
            f"{row['proportion_bootstrap_delta_below_zero']:.4f}"
        )

    metadata = {
        "n_bootstraps_requested": N_BOOTSTRAPS,
        "n_valid_bootstraps": int(valid_bootstraps),
        "n_failed_bootstraps": len(failed),
        "random_attacks_per_bootstrap_network":
            N_RANDOM_ATTACKS_PER_BOOTSTRAP_NETWORK,
        "auc_qmax": AUC_QMAX,
        "katz_alpha": KATZ_ALPHA,
        "ci_radius": CI_RADIUS,
        "r_threshold": R_THRESHOLD,
        "p_threshold_two_sided": P_THRESHOLD,
        "primary_general_statistic":
            "equal-weight mean of paired AUC contrasts across SS, MS, SC, MC",
        "contrast_definition":
            "AUC_betweenness - AUC_comparator",
        "negative_contrast_interpretation":
            "betweenness is more destructive",
        "ci_level": CI_LEVEL,
        "familywise_correction":
            f"Bonferroni across {N_PRIMARY_COMPARISONS} planned comparisons",
    }

    with open(
        OUTPUT_DIR / "bootstrap_analysis_metadata.json",
        "w",
        encoding="utf-8"
    ) as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    print("\nResults saved in:")
    print(OUTPUT_DIR.resolve())

    print("\nMain files:")
    print("  observed_strategy_auc_by_group.csv")
    print("  observed_betweenness_contrasts.csv")
    print("  bootstrap_strategy_auc_by_group.csv")
    print("  bootstrap_betweenness_general_contrast_summary.csv")
    print("  bootstrap_betweenness_group_contrast_summary.csv")

    print("\nANALISI COMPLETATA.")


def run():
    try:
        main()
    except Exception as exc:
        print("\n" + "=" * 78)
        print("ERRORE DURANTE L'ANALISI")
        print("=" * 78)
        print(str(exc))
        print("\nDettaglio tecnico:")
        traceback.print_exc()
    finally:
        try:
            input("\nPremi INVIO per chiudere...")
        except Exception:
            pass


if __name__ == "__main__":
    run()
