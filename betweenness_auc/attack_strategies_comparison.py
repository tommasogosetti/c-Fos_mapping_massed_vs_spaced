#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ATTACK STRATEGIES COMPARISON
Random vs Degree vs Katz vs Betweenness vs Collective Influence (CI)
for the four validated c-Fos covariance networks: SS, MS, SC, MC.

GOAL
----
This is a DESCRIPTIVE head-to-head comparison of how efficiently five node-removal
strategies fragment the SAME network. It is NOT a Spaced-vs-Massed hypothesis test.

For a fair comparison, every strategy:
  1) starts from exactly the same thresholded network for a given group;
  2) uses the same initial component filtering (retain components >= 5 nodes);
  3) tracks the same outcome: largest connected component / original GC;
  4) is summarized by normalized AUC over q = 0..0.20 of original GC nodes removed.

Interpretation:
    LOWER AUC = faster fragmentation = more destructive attack
    HIGHER AUC = greater network resilience to that attack

Random attack is stochastic, so it is repeated many times per network. The targeted
attacks are adaptive: centrality is recalculated after each removal.

IMPORTANT
---------
The script first reconstructs each network from the final already-z-scored matrices
and requires an EXACT edge match to the historical GraphML before continuing.

Collective Influence is computed adaptively with radius d=4, using:
    CI_l(i) = (k_i - 1) * sum_{j: dist(i,j)=l} (k_j - 1)
For the cross-strategy comparison, CI starts from the same component-filtered graph
as all other strategies so that AUC values are directly comparable.
"""

from pathlib import Path
import json
import traceback
import warnings

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import t as student_t


# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "attack_strategies_comparison_results"
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
P_THRESHOLD = 0.04  # two-sided, as validated against the historical GraphML
MIN_COMPONENT_SIZE = 5
N_ATTACK_MEASUREMENTS = 30  # initial state + 29 removals
AUC_QMAX = 0.20
AUC_GRID_POINTS = 251

KATZ_ALPHA = 0.06
CI_RADIUS = 4

# Random baseline. 5000 is the final recommended value.
N_RANDOM_ATTACKS = 5000
RANDOM_SEED = 20260903

GROUPS = ["SS", "MS", "SC", "MC"]
TARGETED_STRATEGIES = ["degree", "katz", "betweenness", "collective_influence"]
ALL_STRATEGIES = ["random"] + TARGETED_STRATEGIES


# =============================================================================
# INPUT / NETWORK RECONSTRUCTION
# =============================================================================

def normalize_region_name(name):
    name = str(name).strip()
    return {
        "NACC SHELL": "NACC.SHELL",
        "NACC CORE": "NACC.CORE",
    }.get(name, name)


def read_zscore_matrix(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path.name}")

    df = pd.read_csv(path, sep=",", na_values=["NA", "NaN", "nan", ""])
    df = df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed")]
    df = df.dropna(axis=1, how="all")
    df.columns = [normalize_region_name(c) for c in df.columns]

    candidates = [
        c for c in df.columns
        if str(c).lower() in {"animale", "animal", "subject", "id"}
    ]
    if not candidates:
        raise ValueError(f"{path.name}: subject-ID column not found.")

    df = df.rename(columns={candidates[0]: "animal"})

    for c in df.columns:
        if c != "animal":
            df[c] = pd.to_numeric(df[c], errors="coerce")

    return df


def region_matrix(df):
    return df.drop(columns=["animal"]).copy()


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


def build_thresholded_network(df, name):
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

    nx.relabel_nodes(
        G,
        {n: normalize_region_name(n) for n in G.nodes()},
        copy=False,
    )
    G.name = name
    return nx.Graph(G)


def edge_set(G):
    return {frozenset((u, v)) for u, v in G.edges() if u != v}


def validate_exact_networks(groups):
    rows = []
    reconstructed = {}

    print("\n--- Exact covariance-network validation ---")

    for group in GROUPS:
        recon = build_thresholded_network(groups[group], f"{group}_reconstructed")
        ref = load_original_graphml(REFERENCE_GRAPH_FILES[group], group)

        er = edge_set(ref)
        en = edge_set(recon)
        exact = er == en
        union = er | en
        jaccard = len(er & en) / len(union) if union else 1.0

        print(
            f"{group}: ref={len(er)}, reconstructed={len(en)}, "
            f"Jaccard={jaccard:.6f}, EXACT={exact}"
        )

        rows.append({
            "group": group,
            "reference_edges": len(er),
            "reconstructed_edges": len(en),
            "common_edges": len(er & en),
            "edge_jaccard": jaccard,
            "exact_edge_match": exact,
        })

        if not exact:
            raise RuntimeError(
                f"STOP: {group} does not exactly reproduce the historical GraphML."
            )

        reconstructed[group] = recon

    pd.DataFrame(rows).to_csv(
        OUTPUT_DIR / "network_validation.csv",
        index=False,
    )

    return reconstructed


# =============================================================================
# COMMON ATTACK BASE GRAPH
# =============================================================================

def purge_components_below_five(G):
    """Retain every connected component with at least MIN_COMPONENT_SIZE nodes."""
    G = nx.Graph(G.copy())
    if G.number_of_nodes() == 0:
        return G

    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    keep_comps = [c for c in comps if len(c) >= MIN_COMPONENT_SIZE]
    if not keep_comps:
        return nx.Graph()

    keep = set().union(*keep_comps)
    return nx.Graph(G.subgraph(keep).copy())


def largest_component_size(G):
    if G.number_of_nodes() == 0:
        return 0
    return max((len(c) for c in nx.connected_components(G)), default=0)


# =============================================================================
# CENTRALITY / ATTACK RULES
# =============================================================================

def highest_degree_node(G):
    ranking = sorted(nx.degree(G), key=lambda x: x[1], reverse=True)
    return ranking[0][0], float(ranking[0][1])


def highest_katz_node(G):
    values = nx.katz_centrality_numpy(
        G,
        alpha=KATZ_ALPHA,
        beta=1.0,
        normalized=True,
        weight=None,
    )
    ranking = sorted(values.items(), key=lambda x: x[1], reverse=True)
    return ranking[0][0], float(ranking[0][1])


def highest_betweenness_node(G):
    values = nx.betweenness_centrality(G, normalized=True, weight=None)
    ranking = sorted(values.items(), key=lambda x: x[1], reverse=True)
    return ranking[0][0], float(ranking[0][1])


def collective_influence_values(G, radius=CI_RADIUS):
    degree = dict(G.degree())
    values = {}

    for node in G.nodes():
        distances = nx.single_source_shortest_path_length(
            G,
            node,
            cutoff=radius,
        )
        boundary_sum = sum(
            degree[other] - 1
            for other, d in distances.items()
            if d == radius
        )
        values[node] = (degree[node] - 1) * boundary_sum

    return values


def highest_ci_node(G):
    values = collective_influence_values(G, CI_RADIUS)
    ranking = sorted(values.items(), key=lambda x: x[1], reverse=True)
    return ranking[0][0], float(ranking[0][1])


def choose_target_node(G, strategy):
    if strategy == "degree":
        return highest_degree_node(G)
    if strategy == "katz":
        return highest_katz_node(G)
    if strategy == "betweenness":
        return highest_betweenness_node(G)
    if strategy == "collective_influence":
        return highest_ci_node(G)
    raise ValueError(f"Unknown targeted strategy: {strategy}")


# =============================================================================
# CURVES AND AUC
# =============================================================================

def attack_curve_targeted(G_input, strategy):
    G = purge_components_below_five(G_input)
    if G.number_of_nodes() == 0:
        raise ValueError("No component survives the common pre-attack filter.")

    gc0 = largest_component_size(G)
    rows = []

    for step in range(N_ATTACK_MEASUREMENTS):
        gc = largest_component_size(G)
        rows.append({
            "step": step,
            "nodes_removed": step,
            "fraction_gc_nodes_removed": step / gc0,
            "gc_size": gc,
            "relative_gc": gc / gc0,
            "node_removed_next": None,
            "score_removed_next": np.nan,
        })

        if step == N_ATTACK_MEASUREMENTS - 1 or G.number_of_nodes() == 0:
            break

        node, score = choose_target_node(G, strategy)
        rows[-1]["node_removed_next"] = node
        rows[-1]["score_removed_next"] = score
        G.remove_node(node)

    return pd.DataFrame(rows)


def attack_curve_random(G_input, rng):
    G = purge_components_below_five(G_input)
    if G.number_of_nodes() == 0:
        raise ValueError("No component survives the common pre-attack filter.")

    gc0 = largest_component_size(G)
    rows = []

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

        nodes = list(G.nodes())
        node = nodes[int(rng.integers(0, len(nodes)))]
        G.remove_node(node)

    return pd.DataFrame(rows)


def normalized_auc(curve):
    x = curve["fraction_gc_nodes_removed"].to_numpy(dtype=float)
    y = curve["relative_gc"].to_numpy(dtype=float)

    if x[-1] < AUC_QMAX:
        raise ValueError(
            f"Curve reaches only q={x[-1]:.3f}; requested q={AUC_QMAX:.3f}."
        )

    grid = np.linspace(0.0, AUC_QMAX, AUC_GRID_POINTS)
    y_grid = np.interp(grid, x, y)

    if hasattr(np, "trapezoid"):
        auc = np.trapezoid(y_grid, grid)
    else:
        auc = np.trapz(y_grid, grid)

    return float(auc / AUC_QMAX)


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def main():
    print("\n============================================================")
    print(" ATTACK STRATEGIES COMPARISON")
    print(" Random vs Degree vs Katz vs Betweenness vs CI")
    print("============================================================")
    print(f"\nWorking folder: {BASE_DIR}")
    print(f"Random attacks per group: {N_RANDOM_ATTACKS}")
    print(f"Common AUC interval: q=0 to {AUC_QMAX}")
    print(f"Katz alpha: {KATZ_ALPHA}")
    print(f"CI radius: d={CI_RADIUS}")
    print(
        "IMPORTANT: every strategy starts from the SAME component-filtered "
        "network within each group."
    )

    groups = {g: read_zscore_matrix(ZSCORE_FILES[g]) for g in GROUPS}
    for g in GROUPS:
        X = region_matrix(groups[g])
        print(
            f"{g}: n={len(groups[g])}, regions={X.shape[1]}, "
            f"missing={int(X.isna().sum().sum())}"
        )

    networks = validate_exact_networks(groups)

    targeted_rows = []
    random_auc_rows = []
    random_curve_rows = []

    # -------------------------------------------------------------------------
    # Targeted attacks
    # -------------------------------------------------------------------------
    print("\n--- Targeted attacks ---")

    for group in GROUPS:
        G = networks[group]

        for strategy in TARGETED_STRATEGIES:
            curve = attack_curve_targeted(G, strategy)
            auc = normalized_auc(curve)

            targeted_rows.append({
                "group": group,
                "strategy": strategy,
                "auc": auc,
            })

            curve.to_csv(
                OUTPUT_DIR / f"{group}_{strategy}_curve.csv",
                index=False,
            )

            print(f"{group} {strategy:22s} AUC = {auc:.6f}")

    targeted_df = pd.DataFrame(targeted_rows)
    targeted_df.to_csv(
        OUTPUT_DIR / "targeted_attack_auc_by_group.csv",
        index=False,
    )

    # -------------------------------------------------------------------------
    # Random attacks
    # -------------------------------------------------------------------------
    print("\n--- Random attacks ---")

    for gi, group in enumerate(GROUPS):
        G = networks[group]
        rng = np.random.default_rng(RANDOM_SEED + gi)
        curves = []

        for rep in range(N_RANDOM_ATTACKS):
            curve = attack_curve_random(G, rng)
            auc = normalized_auc(curve)

            random_auc_rows.append({
                "group": group,
                "random_attack": rep + 1,
                "auc": auc,
            })

            # Save the trajectory for summary statistics; only 30 values each.
            tmp = curve[["step", "fraction_gc_nodes_removed", "relative_gc"]].copy()
            tmp["random_attack"] = rep + 1
            curves.append(tmp)

            if (rep + 1) % max(1, N_RANDOM_ATTACKS // 10) == 0:
                print(f"  {group}: {rep+1}/{N_RANDOM_ATTACKS} random attacks")

        group_random = pd.DataFrame(
            [r for r in random_auc_rows if r["group"] == group]
        )

        q025, q975 = np.quantile(group_random["auc"], [0.025, 0.975])
        print(
            f"{group} random AUC: mean={group_random['auc'].mean():.6f}, "
            f"median={group_random['auc'].median():.6f}, "
            f"95% empirical interval=[{q025:.6f}, {q975:.6f}]"
        )

        long_curves = pd.concat(curves, ignore_index=True)
        curve_summary = (
            long_curves
            .groupby(["step", "fraction_gc_nodes_removed"], as_index=False)
            ["relative_gc"]
            .agg(
                mean_relative_gc="mean",
                median_relative_gc="median",
                lower_2_5=lambda s: np.quantile(s, 0.025),
                upper_97_5=lambda s: np.quantile(s, 0.975),
            )
        )
        curve_summary.insert(0, "group", group)
        random_curve_rows.append(curve_summary)

    random_auc_df = pd.DataFrame(random_auc_rows)
    random_auc_df.to_csv(
        OUTPUT_DIR / "random_attack_auc_distribution.csv",
        index=False,
    )

    pd.concat(random_curve_rows, ignore_index=True).to_csv(
        OUTPUT_DIR / "random_attack_curve_summary.csv",
        index=False,
    )

    # -------------------------------------------------------------------------
    # Build one directly comparable table
    # -------------------------------------------------------------------------
    random_summary = (
        random_auc_df.groupby("group")["auc"]
        .agg(random_mean_auc="mean", random_median_auc="median")
        .reset_index()
    )

    random_ci = (
        random_auc_df.groupby("group")["auc"]
        .quantile([0.025, 0.975])
        .unstack()
        .reset_index()
        .rename(columns={0.025: "random_auc_q025", 0.975: "random_auc_q975"})
    )
    random_summary = random_summary.merge(random_ci, on="group")

    comparison_rows = []

    for group in GROUPS:
        rs = random_summary.loc[random_summary["group"] == group].iloc[0]
        comparison_rows.append({
            "group": group,
            "strategy": "random",
            "auc": float(rs["random_mean_auc"]),
            "auc_lower_95": float(rs["random_auc_q025"]),
            "auc_upper_95": float(rs["random_auc_q975"]),
            "note": "mean of random-attack AUC distribution",
        })

        for strategy in TARGETED_STRATEGIES:
            auc = float(
                targeted_df.loc[
                    (targeted_df["group"] == group)
                    & (targeted_df["strategy"] == strategy),
                    "auc",
                ].iloc[0]
            )
            comparison_rows.append({
                "group": group,
                "strategy": strategy,
                "auc": auc,
                "auc_lower_95": np.nan,
                "auc_upper_95": np.nan,
                "note": "deterministic adaptive targeted attack",
            })

    comparison_df = pd.DataFrame(comparison_rows)
    comparison_df.to_csv(
        OUTPUT_DIR / "attack_strategy_auc_by_group.csv",
        index=False,
    )

    # -------------------------------------------------------------------------
    # General descriptive summary across the four group networks
    # -------------------------------------------------------------------------
    general = (
        comparison_df.groupby("strategy", as_index=False)["auc"]
        .agg(
            mean_auc_across_groups="mean",
            sd_auc_across_groups="std",
            min_auc="min",
            max_auc="max",
        )
    )

    # Rank lower AUC as more destructive.
    general["destructiveness_rank_lower_auc_is_stronger"] = (
        general["mean_auc_across_groups"]
        .rank(method="min", ascending=True)
        .astype(int)
    )
    general = general.sort_values("mean_auc_across_groups").reset_index(drop=True)

    general.to_csv(
        OUTPUT_DIR / "attack_strategy_general_summary.csv",
        index=False,
    )

    # -------------------------------------------------------------------------
    # Direct descriptive contrast with betweenness
    # Positive value here means OTHER attack has higher AUC than betweenness,
    # i.e. betweenness is more destructive.
    # -------------------------------------------------------------------------
    bet_rows = []

    for group in GROUPS:
        auc_bet = float(
            comparison_df.loc[
                (comparison_df["group"] == group)
                & (comparison_df["strategy"] == "betweenness"),
                "auc",
            ].iloc[0]
        )

        for other in ["random", "degree", "katz", "collective_influence"]:
            auc_other = float(
                comparison_df.loc[
                    (comparison_df["group"] == group)
                    & (comparison_df["strategy"] == other),
                    "auc",
                ].iloc[0]
            )

            bet_rows.append({
                "group": group,
                "comparison": f"{other}_minus_betweenness",
                "auc_other": auc_other,
                "auc_betweenness": auc_bet,
                "delta_auc_other_minus_betweenness": auc_other - auc_bet,
            })

    bet_df = pd.DataFrame(bet_rows)
    bet_df.to_csv(
        OUTPUT_DIR / "betweenness_vs_other_attacks_by_group.csv",
        index=False,
    )

    bet_general = (
        bet_df.groupby("comparison", as_index=False)["delta_auc_other_minus_betweenness"]
        .agg(
            mean_delta_across_groups="mean",
            sd_delta_across_groups="std",
            min_delta="min",
            max_delta="max",
        )
    )
    bet_general.to_csv(
        OUTPUT_DIR / "betweenness_vs_other_attacks_general_summary.csv",
        index=False,
    )

    # -------------------------------------------------------------------------
    # Console summary
    # -------------------------------------------------------------------------
    print("\n============================================================")
    print(" GENERAL DESCRIPTIVE RANKING ACROSS SS, MS, SC, MC")
    print(" LOWER AUC = MORE DESTRUCTIVE ATTACK")
    print("============================================================")

    for _, row in general.iterrows():
        print(
            f"#{int(row['destructiveness_rank_lower_auc_is_stronger'])} "
            f"{row['strategy']:22s} mean AUC = {row['mean_auc_across_groups']:.6f}"
        )

    print("\nBetweenness vs other attacks (mean AUC difference across groups):")
    print("Positive delta = the other attack has higher AUC, so betweenness fragments more.")
    for _, row in bet_general.iterrows():
        print(
            f"  {row['comparison']:36s} "
            f"mean delta = {row['mean_delta_across_groups']:.6f}"
        )

    metadata = {
        "analysis": "descriptive attack-strategy comparison",
        "groups": GROUPS,
        "strategies": ALL_STRATEGIES,
        "network_threshold": {"r": R_THRESHOLD, "p_two_sided": P_THRESHOLD},
        "minimum_component_size": MIN_COMPONENT_SIZE,
        "attack_measurements": N_ATTACK_MEASUREMENTS,
        "auc_qmax": AUC_QMAX,
        "katz_alpha": KATZ_ALPHA,
        "ci_radius": CI_RADIUS,
        "random_attacks_per_group": N_RANDOM_ATTACKS,
        "interpretation": "lower AUC = lower resilience / more destructive attack",
        "important_note": (
            "Means across SS/MS/SC/MC are descriptive summaries of four group-level networks, "
            "not an inferential n=4 test. A subject-level paired bootstrap is required for "
            "formal uncertainty on differences between targeted strategies."
        ),
    }
    with open(OUTPUT_DIR / "analysis_metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    print(f"\nResults saved in:\n{OUTPUT_DIR.resolve()}")
    print("\nMain files:")
    print("  attack_strategy_auc_by_group.csv")
    print("  attack_strategy_general_summary.csv")
    print("  random_attack_auc_distribution.csv")
    print("  betweenness_vs_other_attacks_general_summary.csv")


def run():
    try:
        main()
        print("\nANALISI COMPLETATA SENZA ERRORI.")
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
