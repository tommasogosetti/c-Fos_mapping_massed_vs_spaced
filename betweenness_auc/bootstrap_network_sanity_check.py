#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BOOTSTRAP NETWORK SANITY CHECK

Diagnose whether subject-level bootstrap resampling substantially changes the
topology of the c-Fos covariance networks relative to the observed SS, MS, SC,
and MC networks.

Uses the same:
- subject-level z-scored matrices
- pairwise-complete Pearson correlations
- r > 0.65
- two-sided p < 0.04
- 1000 bootstrap iterations
- random seed

For each bootstrap/group records:
- unique animals represented
- edge count
- density
- number of components
- isolated nodes
- giant-component size
- topology after removing components <5 nodes
- edge Jaccard similarity to the observed network
- edge/density/GC ratios relative to the observed network

If the previous bootstrap AUC file exists, topology diagnostics are merged with
the previous attack AUC results and Spearman correlations are calculated.
"""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import os
import traceback

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import t as student_t, spearmanr


BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "bootstrap_network_sanity_results"
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

PREVIOUS_BOOTSTRAP_AUC_FILE = (
    BASE_DIR / "attack_strategy_bootstrap_results" /
    "bootstrap_strategy_auc_by_group.csv"
)

R_THRESHOLD = 0.65
P_THRESHOLD = 0.04
MIN_COMPONENT_SIZE = 5
N_BOOTSTRAPS = 1000
RANDOM_SEED = 20260903
N_JOBS = min(4, max(1, (os.cpu_count() or 2) - 1))


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
                    regions[i], regions[j],
                    r=float(rij),
                    p=float(pij),
                    n=int(n_pair[i, j]),
                )
    return G


def load_reference_graph(path, name):
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


def common_initial_graph(G):
    G = nx.Graph(G.copy())
    comps = list(nx.connected_components(G))
    retained = [c for c in comps if len(c) >= MIN_COMPONENT_SIZE]

    if not retained:
        return nx.Graph()

    keep = set().union(*retained)
    return nx.Graph(G.subgraph(keep).copy())


def graph_diagnostics(G, observed_edge_set=None, observed=None):
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    density = nx.density(G) if n_nodes > 1 else 0.0

    comps = list(nx.connected_components(G))
    n_components = len(comps)
    isolated = sum(1 for n in G.nodes() if G.degree(n) == 0)
    gc_size = max((len(c) for c in comps), default=0)

    filtered = common_initial_graph(G)
    filt_nodes = filtered.number_of_nodes()
    filt_edges = filtered.number_of_edges()
    filt_density = nx.density(filtered) if filt_nodes > 1 else 0.0
    filt_comps = list(nx.connected_components(filtered))
    filt_gc = max((len(c) for c in filt_comps), default=0)
    filt_n_components = len(filt_comps)

    current_edges = edge_set(G)

    if observed_edge_set is not None:
        union = current_edges | observed_edge_set
        edge_jaccard = (
            len(current_edges & observed_edge_set) / len(union)
            if union else 1.0
        )
        edges_in_common = len(current_edges & observed_edge_set)
    else:
        edge_jaccard = np.nan
        edges_in_common = np.nan

    out = {
        "nodes": n_nodes,
        "edges": n_edges,
        "density": density,
        "n_components": n_components,
        "isolated_nodes": isolated,
        "giant_component_size": gc_size,
        "filtered_nodes": filt_nodes,
        "filtered_edges": filt_edges,
        "filtered_density": filt_density,
        "filtered_n_components": filt_n_components,
        "filtered_giant_component_size": filt_gc,
        "edges_in_common_with_observed": edges_in_common,
        "edge_jaccard_to_observed": edge_jaccard,
    }

    if observed is not None:
        out["edge_count_ratio_to_observed"] = (
            n_edges / observed["edges"] if observed["edges"] > 0 else np.nan
        )
        out["density_ratio_to_observed"] = (
            density / observed["density"] if observed["density"] > 0 else np.nan
        )
        out["gc_ratio_to_observed"] = (
            gc_size / observed["giant_component_size"]
            if observed["giant_component_size"] > 0 else np.nan
        )
        out["filtered_gc_ratio_to_observed"] = (
            filt_gc / observed["filtered_giant_component_size"]
            if observed["filtered_giant_component_size"] > 0 else np.nan
        )
    else:
        out["edge_count_ratio_to_observed"] = np.nan
        out["density_ratio_to_observed"] = np.nan
        out["gc_ratio_to_observed"] = np.nan
        out["filtered_gc_ratio_to_observed"] = np.nan

    return out


def prepare_observed(groups):
    observed_rows = []
    observed_info = {}

    print("\n--- Observed networks and exact validation ---")

    for group in ["SS", "MS", "SC", "MC"]:
        recon = build_thresholded_network(
            groups[group],
            f"{group}_observed_reconstructed"
        )
        ref = load_reference_graph(REFERENCE_GRAPH_FILES[group], group)

        exact = edge_set(recon) == edge_set(ref)
        union = edge_set(recon) | edge_set(ref)
        jaccard = (
            len(edge_set(recon) & edge_set(ref)) / len(union)
            if union else 1.0
        )

        if not exact:
            raise RuntimeError(
                f"STOP: {group} reconstructed network does not exactly match "
                "the historical GraphML."
            )

        diag = graph_diagnostics(recon)
        observed_info[group] = {**diag, "edge_set": edge_set(recon)}

        observed_rows.append({
            "group": group,
            "exact_graphml_match": exact,
            "graphml_edge_jaccard": jaccard,
            **diag,
        })

        print(
            f"{group}: EXACT={exact}, edges={diag['edges']}, "
            f"density={diag['density']:.4f}, "
            f"GC={diag['giant_component_size']}, "
            f"filtered GC={diag['filtered_giant_component_size']}"
        )

    observed_df = pd.DataFrame(observed_rows)
    observed_df.to_csv(
        OUTPUT_DIR / "observed_network_topology.csv",
        index=False
    )
    return observed_info, observed_df


WORKER_GROUPS = None
WORKER_OBSERVED = None


def worker_initializer():
    global WORKER_GROUPS, WORKER_OBSERVED

    WORKER_GROUPS = {
        g: read_zscore_matrix(ZSCORE_FILES[g])
        for g in ["SS", "MS", "SC", "MC"]
    }

    WORKER_OBSERVED = {}

    for g in ["SS", "MS", "SC", "MC"]:
        G = build_thresholded_network(
            WORKER_GROUPS[g],
            f"{g}_observed_worker"
        )
        d = graph_diagnostics(G)
        d["edge_set"] = edge_set(G)
        WORKER_OBSERVED[g] = d


def one_bootstrap_iteration(iteration):
    global WORKER_GROUPS, WORKER_OBSERVED

    if WORKER_GROUPS is None or WORKER_OBSERVED is None:
        worker_initializer()

    # Same seed formula as the previous paired attack-strategy bootstrap.
    seed = RANDOM_SEED + 1_000_003 * (iteration + 1)
    rng = np.random.default_rng(seed)

    rows = []

    for group in ["SS", "MS", "SC", "MC"]:
        df = WORKER_GROUPS[group]
        n = len(df)

        idx = rng.integers(0, n, size=n)
        boot_df = df.iloc[idx].reset_index(drop=True)

        n_unique_rows = int(len(np.unique(idx)))
        sampled_labels = df.iloc[idx]["animal"].astype(str).tolist()
        n_unique_animal_labels = int(len(set(sampled_labels)))

        G = build_thresholded_network(
            boot_df,
            f"{group}_bootstrap_{iteration+1}"
        )

        obs = WORKER_OBSERVED[group]
        diag = graph_diagnostics(
            G,
            observed_edge_set=obs["edge_set"],
            observed=obs,
        )

        rows.append({
            "bootstrap": iteration + 1,
            "group": group,
            "n_resampled": n,
            "n_unique_source_rows": n_unique_rows,
            "n_unique_animal_labels": n_unique_animal_labels,
            **diag,
        })

    return rows


SUMMARY_METRICS = [
    "n_unique_source_rows",
    "edges",
    "density",
    "n_components",
    "isolated_nodes",
    "giant_component_size",
    "filtered_nodes",
    "filtered_edges",
    "filtered_density",
    "filtered_n_components",
    "filtered_giant_component_size",
    "edge_jaccard_to_observed",
    "edge_count_ratio_to_observed",
    "density_ratio_to_observed",
    "gc_ratio_to_observed",
    "filtered_gc_ratio_to_observed",
]


def distribution_summary(values):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]

    return {
        "mean": float(np.mean(x)),
        "sd": float(np.std(x, ddof=1)) if len(x) > 1 else np.nan,
        "median": float(np.median(x)),
        "p2_5": float(np.percentile(x, 2.5)),
        "p25": float(np.percentile(x, 25)),
        "p75": float(np.percentile(x, 75)),
        "p97_5": float(np.percentile(x, 97.5)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
        "n": int(len(x)),
    }


def summarize_bootstraps(boot_df, observed_df):
    rows = []

    for group in ["SS", "MS", "SC", "MC"]:
        sub = boot_df[boot_df["group"] == group]
        obsrow = observed_df[observed_df["group"] == group].iloc[0]

        for metric in SUMMARY_METRICS:
            s = distribution_summary(sub[metric].to_numpy())

            if metric == "n_unique_source_rows":
                observed_value = len(sub) * np.nan  # no observed analogue
                observed_value = np.nan
            elif metric in obsrow.index and pd.notna(obsrow[metric]):
                observed_value = float(obsrow[metric])
            else:
                observed_value = np.nan

            rows.append({
                "group": group,
                "metric": metric,
                "observed_value": observed_value,
                **s,
            })

    summary = pd.DataFrame(rows)
    summary.to_csv(
        OUTPUT_DIR / "bootstrap_network_topology_summary.csv",
        index=False
    )
    return summary


def summarize_unique_subject_counts(boot_df):
    counts = (
        boot_df
        .groupby(["group", "n_unique_source_rows"])
        .size()
        .reset_index(name="bootstrap_count")
    )
    totals = counts.groupby("group")["bootstrap_count"].transform("sum")
    counts["proportion"] = counts["bootstrap_count"] / totals

    counts.to_csv(
        OUTPUT_DIR / "bootstrap_unique_subject_count_distribution.csv",
        index=False
    )
    return counts


def auc_topology_relationships(boot_df):
    if not PREVIOUS_BOOTSTRAP_AUC_FILE.exists():
        print(
            "\nPrevious bootstrap AUC file not found; "
            "skipping topology-vs-AUC correlations."
        )
        return None

    print(
        "\nFound previous attack-strategy bootstrap AUC file. "
        "Calculating topology-vs-AUC relationships..."
    )

    auc = pd.read_csv(PREVIOUS_BOOTSTRAP_AUC_FILE)

    required = {
        "bootstrap", "group", "random", "degree", "katz",
        "betweenness", "collective_influence"
    }
    missing = required - set(auc.columns)

    if missing:
        print("AUC file missing columns: " + ", ".join(sorted(missing)))
        return None

    merged = boot_df.merge(
        auc,
        on=["bootstrap", "group"],
        how="inner",
        validate="one_to_one"
    )

    merged["betweenness_minus_random"] = (
        merged["betweenness"] - merged["random"]
    )
    merged["betweenness_minus_degree"] = (
        merged["betweenness"] - merged["degree"]
    )
    merged["betweenness_minus_katz"] = (
        merged["betweenness"] - merged["katz"]
    )
    merged["betweenness_minus_collective_influence"] = (
        merged["betweenness"] - merged["collective_influence"]
    )

    merged.to_csv(
        OUTPUT_DIR / "bootstrap_topology_plus_attack_auc.csv",
        index=False
    )

    topology_vars = [
        "n_unique_source_rows",
        "edges",
        "density",
        "giant_component_size",
        "filtered_giant_component_size",
        "edge_jaccard_to_observed",
        "edge_count_ratio_to_observed",
    ]

    outcome_vars = [
        "random",
        "degree",
        "katz",
        "betweenness",
        "collective_influence",
        "betweenness_minus_random",
        "betweenness_minus_degree",
        "betweenness_minus_katz",
        "betweenness_minus_collective_influence",
    ]

    corr_rows = []

    for group in ["SS", "MS", "SC", "MC"]:
        sub = merged[merged["group"] == group]

        for xvar in topology_vars:
            for yvar in outcome_vars:
                tmp = sub[[xvar, yvar]].dropna()

                if len(tmp) < 3 or tmp[xvar].nunique() < 2 or tmp[yvar].nunique() < 2:
                    rho, p = np.nan, np.nan
                else:
                    rho, p = spearmanr(
                        tmp[xvar].to_numpy(),
                        tmp[yvar].to_numpy()
                    )

                corr_rows.append({
                    "group": group,
                    "topology_variable": xvar,
                    "auc_or_contrast_variable": yvar,
                    "spearman_rho": rho,
                    "p_uncorrected_descriptive": p,
                    "n": len(tmp),
                })

    pd.DataFrame(corr_rows).to_csv(
        OUTPUT_DIR / "bootstrap_topology_auc_spearman_correlations.csv",
        index=False
    )
    return merged


def print_key_diagnostics(summary):
    print("\n============================================================")
    print(" KEY SANITY-CHECK RESULTS")
    print("============================================================")

    for group in ["SS", "MS", "SC", "MC"]:
        print(f"\n--- {group} ---")

        for metric, label in [
            ("edges", "Edges"),
            ("density", "Density"),
            ("giant_component_size", "GC size"),
            ("filtered_giant_component_size", "Filtered GC size"),
            ("edge_jaccard_to_observed", "Edge Jaccard vs observed"),
            ("edge_count_ratio_to_observed", "Edge ratio vs observed"),
            ("n_unique_source_rows", "Unique animals/rows"),
        ]:
            row = summary[
                (summary["group"] == group) &
                (summary["metric"] == metric)
            ].iloc[0]

            obs = row["observed_value"]
            obs_text = f"{obs:.6f}" if np.isfinite(obs) else "n/a"

            print(
                f"{label:26s}: observed={obs_text:>10s} | "
                f"bootstrap median={row['median']:.6f} | "
                f"95% interval=[{row['p2_5']:.6f}, {row['p97_5']:.6f}]"
            )

    print("\nInterpretation guide:")
    print("  edge ratio ~1.0  -> bootstrap edge count near observed")
    print("  ratio >>1.0      -> bootstrap networks much denser")
    print("  ratio <<1.0      -> bootstrap networks much sparser")
    print("  Jaccard near 1   -> bootstrap edges resemble observed")
    print("  Jaccard near 0   -> bootstrap topology differs strongly")


def main():
    print("\n============================================================")
    print(" BOOTSTRAP NETWORK SANITY CHECK")
    print(" Topological stability of subject-level bootstrap networks")
    print("============================================================")

    print(f"\nWorking folder: {BASE_DIR}")
    print(f"Bootstrap iterations: {N_BOOTSTRAPS}")
    print(f"Parallel workers: {N_JOBS}")
    print(
        f"Network criterion: r > {R_THRESHOLD}, "
        f"two-sided p < {P_THRESHOLD}"
    )

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

    _, observed_df = prepare_observed(groups)

    print("\n--- Running topology-only bootstrap diagnostics ---")

    all_rows = []
    failures = []

    if N_JOBS == 1:
        worker_initializer()
        for b in range(N_BOOTSTRAPS):
            try:
                all_rows.extend(one_bootstrap_iteration(b))
            except Exception as exc:
                failures.append({"bootstrap": b + 1, "error": str(exc)})

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
                    all_rows.extend(future.result())
                except Exception as exc:
                    failures.append({"bootstrap": b + 1, "error": str(exc)})

                if completed % max(1, N_BOOTSTRAPS // 20) == 0:
                    print(f"  {completed}/{N_BOOTSTRAPS} bootstraps completed")

    if not all_rows:
        raise RuntimeError("No valid bootstrap diagnostics were produced.")

    boot_df = pd.DataFrame(all_rows).sort_values(
        ["bootstrap", "group"]
    )

    boot_df.to_csv(
        OUTPUT_DIR / "bootstrap_network_diagnostics_by_iteration.csv",
        index=False
    )

    if failures:
        pd.DataFrame(failures).to_csv(
            OUTPUT_DIR / "bootstrap_network_diagnostic_failures.csv",
            index=False
        )

    summary = summarize_bootstraps(boot_df, observed_df)
    summarize_unique_subject_counts(boot_df)
    auc_topology_relationships(boot_df)
    print_key_diagnostics(summary)

    metadata = {
        "n_bootstraps": N_BOOTSTRAPS,
        "random_seed": RANDOM_SEED,
        "r_threshold": R_THRESHOLD,
        "p_threshold_two_sided": P_THRESHOLD,
        "min_component_size": MIN_COMPONENT_SIZE,
        "purpose": (
            "Assess whether subject-level bootstrap resampling changes network "
            "topology enough to affect interpretation of paired attack-strategy "
            "bootstrap confidence intervals."
        ),
    }

    with open(
        OUTPUT_DIR / "bootstrap_network_sanity_metadata.json",
        "w",
        encoding="utf-8"
    ) as f:
        json.dump(metadata, f, indent=2)

    print("\n============================================================")
    print(" ANALYSIS COMPLETE")
    print("============================================================")
    print(f"\nResults saved in:\n{OUTPUT_DIR.resolve()}")
    print("\nMain files:")
    print("  observed_network_topology.csv")
    print("  bootstrap_network_diagnostics_by_iteration.csv")
    print("  bootstrap_network_topology_summary.csv")
    print("  bootstrap_unique_subject_count_distribution.csv")

    if PREVIOUS_BOOTSTRAP_AUC_FILE.exists():
        print("  bootstrap_topology_plus_attack_auc.csv")
        print("  bootstrap_topology_auc_spearman_correlations.csv")


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
