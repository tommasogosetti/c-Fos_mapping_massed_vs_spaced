#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BETWEENNESS ROBUSTNESS ANALYSIS
Updated version for the FINAL, ALREADY Z-SCORED subject x region matrices.

Put this script in the SAME folder as:

    ss net.csv
    ms net.csv
    sc net.csv
    mc net.csv

and the original reference networks:

    ssnet 0.65 0.04.graphml
    msnet 0.65.graphml
    scnet 0.65.graphml
    mcnet_0.65.graphml

IMPORTANT
---------
The four CSV files supplied to this version are ALREADY standardized within subject.
The script DOES NOT z-score them again.

The script:
1. checks that each supplied row is approximately mean=0, SD=1;
2. reconstructs each positive covariance network with Pearson r > 0.65 AND
   two-sided p < 0.04;
3. validates the reconstructed networks against the original GraphML files;
4. reproduces the original adaptive, UNWEIGHTED betweenness attack;
5. calculates q50 and normalized AUC of the giant-component decay curve;
6. can run a subject-label permutation test of the AUC difference.

The permutation test pools the ALREADY Z-SCORED subjects, permutes subject labels,
reconstructs both networks from scratch, reruns the adaptive betweenness attack,
and compares Delta AUC between spaced and massed groups.
"""

from pathlib import Path
import json
import sys
import traceback
import warnings

try:
    import numpy as np
    import pandas as pd
    import networkx as nx
    from scipy.stats import t as student_t
except ImportError as exc:
    print("\nERRORE: manca un pacchetto Python necessario.")
    print(f"Dettaglio: {exc}")
    print("\nInstalla con:")
    print("  python -m pip install numpy pandas scipy networkx")
    input("\nPremi INVIO per chiudere...")
    raise SystemExit(1)


# ============================================================================
# 1. CONFIGURATION
# ============================================================================

BASE_DIR = Path(__file__).resolve().parent

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

OUTPUT_DIR = BASE_DIR / "betweenness_robustness_results"

# Exact network-construction rule that reproduces the supplied GraphML networks
R_THRESHOLD = 0.65
P_THRESHOLD = 0.04
P_SIDEDNESS = "two-sided"

# Original disruption settings
MIN_COMPONENT_SIZE = 5
N_ATTACK_MEASUREMENTS = 30   # original code: 30 GC measurements, therefore 29 removals

# AUC is compared on the same x-domain in every network
AUC_QMAX = 0.20
AUC_GRID_POINTS = 251

# First run: validation only.
# After validation reports exact reference matching, change this to True.
RUN_PERMUTATIONS = True

# Start with e.g. 100 only as a smoke test.
# For the final inferential analysis use 5000.
N_PERMUTATIONS = 5000
RANDOM_SEED = 20260902

# Strict gate: permutation testing is allowed only if the observed networks
# reconstructed from the z-score matrices exactly match the original GraphML.
REQUIRE_EXACT_REFERENCE_MATCH = True


# ============================================================================
# 2. INPUT AND LABEL NORMALIZATION
# ============================================================================

def normalize_region_name(name: str) -> str:
    """
    Match historical GraphML node labels.

    The supplied final z-score CSVs use:
       NACC SHELL, NACC CORE
    whereas the historical GraphML files use:
       NACC.SHELL, NACC.CORE
    """
    name = str(name).strip()
    mapping = {
        "NACC SHELL": "NACC.SHELL",
        "NACC CORE": "NACC.CORE",
    }
    return mapping.get(name, name)


def read_zscore_matrix(path: Path) -> pd.DataFrame:
    """Read one final already-z-scored animal x region matrix."""
    if not path.exists():
        raise FileNotFoundError(f"File non trovato: {path.name}")

    # These final matrices are comma-separated.
    df = pd.read_csv(path, sep=",", na_values=["NA", "NaN", "nan", ""])

    df = df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed")]
    df = df.dropna(axis=1, how="all")

    # Clean/normalize column names
    df.columns = [normalize_region_name(c) for c in df.columns]

    subject_candidates = [
        c for c in df.columns
        if str(c).lower() in {"animale", "animal", "subject", "id"}
    ]
    if not subject_candidates:
        raise ValueError(
            f"{path.name}: non trovo la colonna identificativa 'animale'."
        )

    subject_col = subject_candidates[0]
    df = df.rename(columns={subject_col: "animal"})

    if df["animal"].duplicated().any():
        raise ValueError(f"{path.name}: ID animale duplicati.")

    for c in df.columns:
        if c != "animal":
            df[c] = pd.to_numeric(df[c], errors="coerce")

    if df.columns.duplicated().any():
        dup = df.columns[df.columns.duplicated()].tolist()
        raise ValueError(f"{path.name}: colonne duplicate dopo normalizzazione: {dup}")

    return df


def region_matrix(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop(columns=["animal"]).copy()


def verify_already_zscored(df: pd.DataFrame, group: str) -> pd.DataFrame:
    """Check row-wise mean and SD without altering the data."""
    X = region_matrix(df)

    means = X.mean(axis=1, skipna=True)
    sds = X.std(axis=1, skipna=True, ddof=1)

    check = pd.DataFrame({
        "group": group,
        "animal": df["animal"].astype(str),
        "row_mean": means,
        "row_sd": sds,
        "abs_mean": means.abs(),
        "abs_sd_minus_1": (sds - 1.0).abs(),
    })

    # These files should be extremely close to mean=0 and SD=1.
    if (check["abs_mean"] > 1e-5).any() or (check["abs_sd_minus_1"] > 1e-5).any():
        warnings.warn(
            f"{group}: alcune righe non risultano esattamente z-standardizzate "
            "entro la tolleranza 1e-5. I dati NON vengono modificati."
        )

    return check


# ============================================================================
# 3. PAIRWISE PEARSON r AND p
# ============================================================================

def pairwise_pearson_r_p(X: pd.DataFrame):
    """
    Pairwise-complete Pearson correlations.

    Returns:
      r      : Pearson correlation matrix
      p_two  : conventional two-sided p-value matrix
      n_pair : number of available subjects for each region pair
    """
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


# ============================================================================
# 4. NETWORK CONSTRUCTION
# ============================================================================

def build_thresholded_network(
    df: pd.DataFrame,
    r_threshold: float = R_THRESHOLD,
    p_threshold: float = P_THRESHOLD,
    name: str = "network",
) -> nx.Graph:

    # DO NOT z-score here: supplied matrices are already final z-scores.
    X = region_matrix(df)
    regions = list(X.columns)

    r, p_two, n_pair = pairwise_pearson_r_p(X)

    G = nx.Graph(name=name)
    G.add_nodes_from(regions)

    p_regions = len(regions)

    for i in range(p_regions):
        for j in range(i + 1, p_regions):
            rij = r[i, j]
            pij = p_two[i, j]
            nij = int(n_pair[i, j]) if np.isfinite(n_pair[i, j]) else 0

            if (
                np.isfinite(rij)
                and np.isfinite(pij)
                and rij > r_threshold
                and pij < p_threshold
            ):
                # r/p retained only for auditing.
                # Betweenness disruption itself is UNWEIGHTED.
                G.add_edge(
                    regions[i],
                    regions[j],
                    r=float(rij),
                    p=float(pij),
                    n=nij,
                )

    return G


# ============================================================================
# 5. REFERENCE GRAPHML
# ============================================================================

def load_original_graphml(path: Path, name: str) -> nx.Graph:
    if not path.exists():
        raise FileNotFoundError(f"GraphML non trovato: {path.name}")

    G = nx.read_graphml(path).to_undirected()
    G.name = name

    attrs = nx.get_node_attributes(G, "shared name")
    if attrs:
        nx.relabel_nodes(G, attrs, copy=False)
    else:
        name_attrs = nx.get_node_attributes(G, "name")
        if name_attrs:
            nx.relabel_nodes(G, name_attrs, copy=False)

    # Normalize labels on the historical graph too, defensively.
    mapping = {n: normalize_region_name(n) for n in G.nodes()}
    nx.relabel_nodes(G, mapping, copy=False)

    return nx.Graph(G)


def edge_set(G: nx.Graph):
    return {frozenset((u, v)) for u, v in G.edges() if u != v}


def compare_graphs(reference: nx.Graph, reconstructed: nx.Graph) -> dict:
    e_ref = edge_set(reference)
    e_new = edge_set(reconstructed)

    inter = len(e_ref & e_new)
    union = len(e_ref | e_new)

    return {
        "reference_nodes": reference.number_of_nodes(),
        "reference_edges": reference.number_of_edges(),
        "reconstructed_nodes": reconstructed.number_of_nodes(),
        "reconstructed_edges": reconstructed.number_of_edges(),
        "common_edges": inter,
        "reference_only_edges": len(e_ref - e_new),
        "reconstructed_only_edges": len(e_new - e_ref),
        "edge_union": union,
        "edge_jaccard": inter / union if union else 1.0,
        "exact_edge_match": e_ref == e_new,
    }


# ============================================================================
# 6. ORIGINAL COMPONENT PURGE
# ============================================================================

def purge_smaller_components(
    G: nx.Graph,
    minimum_nodes: int = MIN_COMPONENT_SIZE
) -> nx.Graph:
    """
    Reproduce the historical script:
    retain all connected components having at least `minimum_nodes`.
    """
    if G.number_of_nodes() == 0:
        return G.copy()

    components = sorted(nx.connected_components(G), key=len, reverse=True)
    retained = [c for c in components if len(c) >= minimum_nodes]

    if not retained:
        return nx.Graph()

    keep = set().union(*retained)
    return nx.Graph(G.subgraph(keep).copy())


def largest_component_size(G: nx.Graph) -> int:
    if G.number_of_nodes() == 0:
        return 0
    return max((len(c) for c in nx.connected_components(G)), default=0)


# ============================================================================
# 7. ADAPTIVE UNWEIGHTED BETWEENNESS ATTACK
# ============================================================================

def sorted_nodes_by_betweenness(G: nx.Graph):
    """
    Exact historical logic:
      sorted(nx.betweenness_centrality(network).items(),
             key=lambda x: x[1], reverse=True)

    weight=None => unweighted shortest paths.
    """
    bc = nx.betweenness_centrality(G, normalized=True, weight=None)
    return sorted(bc.items(), key=lambda x: x[1], reverse=True)


def betweenness_attack(
    G_input: nx.Graph,
    number_of_measurements: int = N_ATTACK_MEASUREMENTS,
    minimum_component_size: int = MIN_COMPONENT_SIZE,
) -> pd.DataFrame:

    # Purge once before attack, as in the historical code.
    G = purge_smaller_components(
        nx.Graph(G_input.copy()),
        minimum_component_size
    )

    if G.number_of_nodes() == 0:
        raise ValueError("Nessuna componente sopravvive al purge iniziale.")

    gc0 = largest_component_size(G)
    n0 = G.number_of_nodes()

    rows = []

    # Historical code records 30 GC measurements and removes a node after the
    # first 29 measurements.
    for step in range(number_of_measurements):

        gc = largest_component_size(G)

        row = {
            "step": step,
            "nodes_removed": step,
            "fraction_gc_nodes_removed": step / gc0,
            "fraction_retained_nodes_removed": step / n0,
            "gc_size": gc,
            "relative_gc": gc / gc0,
            "node_removed_next": None,
            "betweenness_removed_next": np.nan,
        }

        rows.append(row)

        if step == number_of_measurements - 1 or G.number_of_nodes() == 0:
            break

        ranking = sorted_nodes_by_betweenness(G)
        if not ranking:
            break

        node, bc = ranking[0]

        rows[-1]["node_removed_next"] = node
        rows[-1]["betweenness_removed_next"] = bc

        G.remove_node(node)

    return pd.DataFrame(rows)


# ============================================================================
# 8. ROBUSTNESS SUMMARY: q50 AND AUC
# ============================================================================

def q50_from_curve(curve: pd.DataFrame) -> float:
    hit = curve.loc[curve["relative_gc"] <= 0.5]
    if hit.empty:
        return np.nan
    return float(hit.iloc[0]["fraction_gc_nodes_removed"])


def nodes_to_q50(curve: pd.DataFrame):
    hit = curve.loc[curve["relative_gc"] <= 0.5]
    if hit.empty:
        return np.nan
    return int(hit.iloc[0]["nodes_removed"])


def auc_common_fraction(
    curve: pd.DataFrame,
    qmax: float = AUC_QMAX,
    grid_points: int = AUC_GRID_POINTS,
) -> float:

    x = curve["fraction_gc_nodes_removed"].to_numpy(dtype=float)
    y = curve["relative_gc"].to_numpy(dtype=float)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    if x[-1] < qmax:
        raise ValueError(
            f"La curva arriva solo a q={x[-1]:.3f}; "
            f"AUC_QMAX={qmax:.3f}. Aumenta N_ATTACK_MEASUREMENTS o riduci AUC_QMAX."
        )

    grid = np.linspace(0.0, qmax, grid_points)
    y_grid = np.interp(grid, x, y)

    if hasattr(np, "trapezoid"):
        auc = np.trapezoid(y_grid, grid)
    else:
        auc = np.trapz(y_grid, grid)

    # Normalize by the x-domain length.
    # Larger = more robust; smaller = more vulnerable.
    return float(auc / qmax)


def summarize_curve(curve: pd.DataFrame, group: str) -> dict:
    return {
        "group": group,
        "gc0": int(curve.iloc[0]["gc_size"]),
        "nodes_to_50pct_gc": nodes_to_q50(curve),
        "q50": q50_from_curve(curve),
        "normalized_auc_q0_to_0.20": auc_common_fraction(curve),
    }


# ============================================================================
# 9. BUILD + ATTACK + AUC
# ============================================================================

def build_and_auc(df: pd.DataFrame, name: str):
    G = build_thresholded_network(df, name=name)
    curve = betweenness_attack(G)
    auc = auc_common_fraction(curve)
    return auc, curve, G


# ============================================================================
# 10. SUBJECT-LABEL PERMUTATION TEST
# ============================================================================

def permutation_auc_test(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    label_a: str,
    label_b: str,
    n_perm: int = N_PERMUTATIONS,
    seed: int = RANDOM_SEED,
):

    if list(region_matrix(df_a).columns) != list(region_matrix(df_b).columns):
        raise ValueError(f"Colonne regionali diverse tra {label_a} e {label_b}.")

    n_a = len(df_a)
    n_b = len(df_b)

    # Observed networks reconstructed with exactly the same method as permutations.
    auc_a, curve_a, _ = build_and_auc(df_a, label_a)
    auc_b, curve_b, _ = build_and_auc(df_b, label_b)

    delta_obs = auc_a - auc_b

    pooled = pd.concat([df_a, df_b], ignore_index=True)
    rng = np.random.default_rng(seed)

    deltas = np.full(n_perm, np.nan)
    auc_perm_a = np.full(n_perm, np.nan)
    auc_perm_b = np.full(n_perm, np.nan)

    all_idx = np.arange(len(pooled))
    failures = []

    for b in range(n_perm):
        idx = rng.permutation(all_idx)

        pseudo_a = pooled.iloc[idx[:n_a]].reset_index(drop=True)
        pseudo_b = pooled.iloc[idx[n_a:n_a + n_b]].reset_index(drop=True)

        try:
            auc1, _, _ = build_and_auc(
                pseudo_a, f"perm_{label_a}_{b+1}"
            )
            auc2, _, _ = build_and_auc(
                pseudo_b, f"perm_{label_b}_{b+1}"
            )

            auc_perm_a[b] = auc1
            auc_perm_b[b] = auc2
            deltas[b] = auc1 - auc2

        except Exception as exc:
            failures.append((b + 1, str(exc)))

        # Progress every 5%
        if (b + 1) % max(1, n_perm // 20) == 0:
            print(
                f"  {label_a} vs {label_b}: "
                f"{b+1}/{n_perm} permutations"
            )

    valid = np.isfinite(deltas)

    if valid.sum() == 0:
        raise RuntimeError("Tutte le permutation sono fallite.")

    # Primary two-sided permutation p-value
    p_two = (
        1 + np.sum(np.abs(deltas[valid]) >= abs(delta_obs))
    ) / (
        1 + valid.sum()
    )

    perm_df = pd.DataFrame({
        "permutation": np.arange(1, n_perm + 1),
        f"auc_{label_a}": auc_perm_a,
        f"auc_{label_b}": auc_perm_b,
        "delta_auc_A_minus_B": deltas,
    })

    result = {
        "comparison": f"{label_a}_vs_{label_b}",
        f"observed_auc_{label_a}": float(auc_a),
        f"observed_auc_{label_b}": float(auc_b),
        "observed_delta_auc_A_minus_B": float(delta_obs),
        "permutations_requested": int(n_perm),
        "permutations_valid": int(valid.sum()),
        "permutations_failed": int((~valid).sum()),
        "p_two_sided": float(p_two),
        "interpretation": (
            "Smaller AUC = greater network vulnerability; "
            "negative Delta AUC means group A is more vulnerable than group B."
        ),
        "first_failures": failures[:20],
    }

    return result, perm_df


# ============================================================================
# 11. MAIN
# ============================================================================

def main():

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n============================================================")
    print("  BETWENNESS ROBUSTNESS ANALYSIS - FINAL Z-SCORE MATRICES")
    print("============================================================\n")

    print(f"Working folder: {BASE_DIR}")
    print("Input matrices are ALREADY z-scored; no additional standardization.")
    print(f"Network threshold: r > {R_THRESHOLD}, two-sided p < {P_THRESHOLD}")
    print(f"Pre-attack minimum component size: {MIN_COMPONENT_SIZE}")
    print(f"Original attack measurements: {N_ATTACK_MEASUREMENTS}")
    print(f"AUC common interval: q = 0 to {AUC_QMAX}")
    print(f"RUN_PERMUTATIONS = {RUN_PERMUTATIONS}\n")

    # ----------------------------------------------------------------------
    # Load final z-score matrices
    # ----------------------------------------------------------------------
    raw = {}
    z_checks = []

    for group in ["SS", "MS", "SC", "MC"]:
        df = read_zscore_matrix(ZSCORE_FILES[group])
        raw[group] = df

        X = region_matrix(df)
        check = verify_already_zscored(df, group)
        z_checks.append(check)

        print(
            f"{group}: n={len(df)}, regions={X.shape[1]}, "
            f"missing={int(X.isna().sum().sum())}, "
            f"max|row mean|={check['abs_mean'].max():.3g}, "
            f"max|row SD-1|={check['abs_sd_minus_1'].max():.3g}"
        )

    pd.concat(z_checks, ignore_index=True).to_csv(
        OUTPUT_DIR / "zscore_input_check.csv",
        index=False
    )

    # Ensure identical region columns
    canonical = list(region_matrix(raw["SS"]).columns)
    for group in ["MS", "SC", "MC"]:
        if list(region_matrix(raw[group]).columns) != canonical:
            raise ValueError(
                f"L'ordine/nome delle regioni differisce tra SS e {group}."
            )

    # ----------------------------------------------------------------------
    # Reconstruct networks from the supplied final z-score matrices
    # ----------------------------------------------------------------------
    reconstructed = {}

    print("\n--- Reconstructing networks from final z-score matrices ---")

    for group in ["SS", "MS", "SC", "MC"]:
        G = build_thresholded_network(
            raw[group],
            name=f"{group}_zscore_reconstructed"
        )
        reconstructed[group] = G

        print(
            f"{group}: nodes={G.number_of_nodes()}, "
            f"edges={G.number_of_edges()}, "
            f"GC={largest_component_size(G)}"
        )

        nx.write_graphml(
            G,
            OUTPUT_DIR / f"{group}_zscore_reconstructed.graphml"
        )

    # ----------------------------------------------------------------------
    # Validate EXACTLY against historical GraphML
    # ----------------------------------------------------------------------
    validation_rows = []
    observed_summaries = []
    references = {}

    print("\n--- Exact GraphML validation ---")

    for group in ["SS", "MS", "SC", "MC"]:

        ref = load_original_graphml(
            REFERENCE_GRAPH_FILES[group],
            group
        )
        references[group] = ref

        comp = compare_graphs(
            ref,
            reconstructed[group]
        )
        comp["group"] = group
        validation_rows.append(comp)

        # Attack the historical reference graph
        curve = betweenness_attack(ref)
        summary = summarize_curve(curve, group)
        observed_summaries.append(summary)

        curve.insert(0, "group", group)
        curve.to_csv(
            OUTPUT_DIR / f"{group}_reference_betweenness_GC_curve.csv",
            index=False
        )

        curve.loc[
            curve["node_removed_next"].notna(),
            [
                "step",
                "node_removed_next",
                "betweenness_removed_next",
                "gc_size",
                "relative_gc",
            ]
        ].to_csv(
            OUTPUT_DIR / f"{group}_reference_betweenness_removal_order.csv",
            index=False
        )

        first7 = (
            curve.loc[
                curve["node_removed_next"].notna(),
                "node_removed_next"
            ]
            .head(7)
            .tolist()
        )

        print(
            f"{group}: ref edges={comp['reference_edges']}, "
            f"reconstructed edges={comp['reconstructed_edges']}, "
            f"common={comp['common_edges']}, "
            f"Jaccard={comp['edge_jaccard']:.6f}, "
            f"EXACT={comp['exact_edge_match']}"
        )
        print(
            f"    GC0={summary['gc0']}, "
            f"nodes-to-50%={summary['nodes_to_50pct_gc']}, "
            f"q50={summary['q50']:.5f}, "
            f"AUC={summary['normalized_auc_q0_to_0.20']:.5f}"
        )
        print(f"    first 7 removals: {first7}")

    validation_df = pd.DataFrame(validation_rows)
    validation_df.to_csv(
        OUTPUT_DIR / "reference_vs_final_zscore_validation.csv",
        index=False
    )

    pd.DataFrame(observed_summaries).to_csv(
        OUTPUT_DIR / "reference_robustness_summary.csv",
        index=False
    )

    exact_all = bool(validation_df["exact_edge_match"].all())

    print("\nExact reconstruction of all four historical networks:", exact_all)

    # ----------------------------------------------------------------------
    # Permutation tests
    # ----------------------------------------------------------------------
    if RUN_PERMUTATIONS:

        if REQUIRE_EXACT_REFERENCE_MATCH and not exact_all:
            raise RuntimeError(
                "STOP: almeno una rete ricostruita non coincide esattamente "
                "con il GraphML storico. Non eseguo le permutation."
            )

        print("\n============================================================")
        print("SUBJECT-LABEL PERMUTATION TESTS")
        print("============================================================")

        comparisons = [("SS", "MS"), ("SC", "MC")]
        summaries = []

        for i, (a, b) in enumerate(comparisons):

            print(
                f"\n{a} vs {b}: {N_PERMUTATIONS} permutations"
            )

            result, perm_df = permutation_auc_test(
                raw[a],
                raw[b],
                a,
                b,
                n_perm=N_PERMUTATIONS,
                seed=RANDOM_SEED + i,
            )

            summaries.append(
                {
                    k: v for k, v in result.items()
                    if k != "first_failures"
                }
            )

            perm_df.to_csv(
                OUTPUT_DIR / f"permutation_{a}_vs_{b}.csv",
                index=False
            )

            with open(
                OUTPUT_DIR / f"permutation_{a}_vs_{b}_result.json",
                "w",
                encoding="utf-8",
            ) as f:
                json.dump(
                    result,
                    f,
                    indent=2,
                    ensure_ascii=False
                )

            print(
                f"Observed AUC {a}: "
                f"{result[f'observed_auc_{a}']:.6f}"
            )
            print(
                f"Observed AUC {b}: "
                f"{result[f'observed_auc_{b}']:.6f}"
            )
            print(
                f"Delta AUC ({a}-{b}): "
                f"{result['observed_delta_auc_A_minus_B']:.6f}"
            )
            print(
                f"Two-sided permutation p: "
                f"{result['p_two_sided']:.6g}"
            )

        pd.DataFrame(summaries).to_csv(
            OUTPUT_DIR / "permutation_summary.csv",
            index=False
        )

    else:
        print("\nPermutation test NON eseguito.")
        print(
            "Se sopra tutte e quattro le reti riportano EXACT=True, "
            "la pipeline è validata."
        )
        print(
            "Per il test successivo modifica nello script:"
        )
        print("    RUN_PERMUTATIONS = True")
        print("e, per l'analisi finale:")
        print("    N_PERMUTATIONS = 5000")

    print(
        f"\nRisultati salvati in:\n{OUTPUT_DIR.resolve()}"
    )


def run_with_friendly_console():
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
        print(
            "\nNon chiudere la finestra: fai uno screenshot "
            "dell'errore e mandamelo."
        )
    finally:
        try:
            input("\nPremi INVIO per chiudere...")
        except Exception:
            pass


if __name__ == "__main__":
    run_with_friendly_console()
