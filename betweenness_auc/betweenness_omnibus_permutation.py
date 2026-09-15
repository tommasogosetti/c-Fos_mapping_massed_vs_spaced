#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Stratified omnibus permutation test for betweenness-targeted robustness.

Place this file in the SAME folder as:
- betweenness_robustness.py
- ss net.csv, ms net.csv, sc net.csv, mc net.csv
- the four original GraphML files

It reuses exactly the same network construction, GraphML validation,
betweenness attack and AUC functions from betweenness_robustness.py.

Primary omnibus statistic:
mean[(AUC_SS - AUC_MS), (AUC_SC - AUC_MC)]

Negative values indicate lower AUC (greater vulnerability) in spaced networks.
"""

from pathlib import Path
import json
import traceback
import numpy as np
import pandas as pd
import betweenness_robustness as br

BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "betweenness_robustness_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

N_PERMUTATIONS = 5000
RANDOM_SEED = 20260902


def load_groups():
    return {
        g: br.read_zscore_matrix(br.ZSCORE_FILES[g])
        for g in ["SS", "MS", "SC", "MC"]
    }


def validate_exact(groups):
    rows = []
    all_exact = True

    print("\n--- Exact GraphML validation ---")
    for g in ["SS", "MS", "SC", "MC"]:
        recon = br.build_thresholded_network(groups[g], name=f"{g}_omnibus")
        ref = br.load_original_graphml(br.REFERENCE_GRAPH_FILES[g], g)
        cmp = br.compare_graphs(ref, recon)
        cmp["group"] = g
        rows.append(cmp)

        print(
            f"{g}: ref={cmp['reference_edges']}, "
            f"reconstructed={cmp['reconstructed_edges']}, "
            f"Jaccard={cmp['edge_jaccard']:.6f}, "
            f"EXACT={cmp['exact_edge_match']}"
        )

        all_exact = all_exact and bool(cmp["exact_edge_match"])

    pd.DataFrame(rows).to_csv(
        OUTPUT_DIR / "omnibus_reference_validation.csv",
        index=False
    )

    if not all_exact:
        raise RuntimeError(
            "STOP: almeno una rete non coincide esattamente con il GraphML storico."
        )

    return all_exact


def auc_for(df, name):
    auc, _, _ = br.build_and_auc(df, name)
    return auc


def main():
    print("\n============================================================")
    print(" STRATIFIED OMNIBUS PERMUTATION TEST - BETWENNESS AUC")
    print("============================================================")
    print(f"\nWorking folder: {BASE_DIR}")
    print(f"Permutations: {N_PERMUTATIONS}")
    print("Primary statistic: mean[(AUC_SS-AUC_MS), (AUC_SC-AUC_MC)]")
    print("Negative = spaced networks are overall more vulnerable.")

    groups = load_groups()

    for g in ["SS", "MS", "SC", "MC"]:
        X = br.region_matrix(groups[g])
        print(
            f"{g}: n={len(groups[g])}, regions={X.shape[1]}, "
            f"missing={int(X.isna().sum().sum())}"
        )

    validate_exact(groups)

    auc_ss = auc_for(groups["SS"], "SS_obs")
    auc_ms = auc_for(groups["MS"], "MS_obs")
    auc_sc = auc_for(groups["SC"], "SC_obs")
    auc_mc = auc_for(groups["MC"], "MC_obs")

    delta_spatial_obs = auc_ss - auc_ms
    delta_cue_obs = auc_sc - auc_mc
    omnibus_obs = (delta_spatial_obs + delta_cue_obs) / 2.0

    print("\n--- Observed values ---")
    print(f"AUC SS = {auc_ss:.6f}")
    print(f"AUC MS = {auc_ms:.6f}")
    print(f"AUC SC = {auc_sc:.6f}")
    print(f"AUC MC = {auc_mc:.6f}")
    print(f"Delta spatial (SS-MS) = {delta_spatial_obs:.6f}")
    print(f"Delta cue (SC-MC) = {delta_cue_obs:.6f}")
    print(f"OMNIBUS mean delta AUC = {omnibus_obs:.6f}")

    spatial_pool = pd.concat([groups["SS"], groups["MS"]], ignore_index=True)
    cue_pool = pd.concat([groups["SC"], groups["MC"]], ignore_index=True)

    n_ss, n_ms = len(groups["SS"]), len(groups["MS"])
    n_sc, n_mc = len(groups["SC"]), len(groups["MC"])

    spatial_idx = np.arange(len(spatial_pool))
    cue_idx = np.arange(len(cue_pool))

    rng = np.random.default_rng(RANDOM_SEED)

    delta_spatial_perm = np.full(N_PERMUTATIONS, np.nan)
    delta_cue_perm = np.full(N_PERMUTATIONS, np.nan)
    omnibus_perm = np.full(N_PERMUTATIONS, np.nan)
    failures = []

    print("\n--- Running stratified omnibus permutations ---")

    for b in range(N_PERMUTATIONS):
        try:
            ps = rng.permutation(spatial_idx)
            pc = rng.permutation(cue_idx)

            p_ss = spatial_pool.iloc[ps[:n_ss]].reset_index(drop=True)
            p_ms = spatial_pool.iloc[ps[n_ss:n_ss+n_ms]].reset_index(drop=True)
            p_sc = cue_pool.iloc[pc[:n_sc]].reset_index(drop=True)
            p_mc = cue_pool.iloc[pc[n_sc:n_sc+n_mc]].reset_index(drop=True)

            a_ss = auc_for(p_ss, f"perm_SS_{b+1}")
            a_ms = auc_for(p_ms, f"perm_MS_{b+1}")
            a_sc = auc_for(p_sc, f"perm_SC_{b+1}")
            a_mc = auc_for(p_mc, f"perm_MC_{b+1}")

            ds = a_ss - a_ms
            dc = a_sc - a_mc
            om = (ds + dc) / 2.0

            delta_spatial_perm[b] = ds
            delta_cue_perm[b] = dc
            omnibus_perm[b] = om

        except Exception as exc:
            failures.append((b + 1, str(exc)))

        if (b + 1) % max(1, N_PERMUTATIONS // 20) == 0:
            print(f"  Omnibus: {b+1}/{N_PERMUTATIONS} permutations")

    valid = np.isfinite(omnibus_perm)

    if valid.sum() == 0:
        raise RuntimeError("Tutte le permutation sono fallite.")

    # Primary two-sided permutation p-value (+1 correction)
    p_two = (
        1 + np.sum(np.abs(omnibus_perm[valid]) >= abs(omnibus_obs))
    ) / (
        1 + valid.sum()
    )

    # Secondary directional p-value: spaced more vulnerable => negative statistic
    p_one = (
        1 + np.sum(omnibus_perm[valid] <= omnibus_obs)
    ) / (
        1 + valid.sum()
    )

    dist = pd.DataFrame({
        "permutation": np.arange(1, N_PERMUTATIONS + 1),
        "delta_spatial_SS_minus_MS": delta_spatial_perm,
        "delta_cue_SC_minus_MC": delta_cue_perm,
        "mean_delta_AUC": omnibus_perm,
    })
    dist.to_csv(
        OUTPUT_DIR / "omnibus_permutation_distribution.csv",
        index=False
    )

    result = {
        "test": "stratified_omnibus_betweenness_AUC",
        "statistic": "mean[(AUC_SS-AUC_MS),(AUC_SC-AUC_MC)]",
        "auc_SS": float(auc_ss),
        "auc_MS": float(auc_ms),
        "auc_SC": float(auc_sc),
        "auc_MC": float(auc_mc),
        "delta_spatial_SS_minus_MS": float(delta_spatial_obs),
        "delta_cue_SC_minus_MC": float(delta_cue_obs),
        "observed_mean_delta_AUC": float(omnibus_obs),
        "permutations_requested": int(N_PERMUTATIONS),
        "permutations_valid": int(valid.sum()),
        "permutations_failed": int((~valid).sum()),
        "p_two_sided_primary": float(p_two),
        "p_one_sided_spaced_more_vulnerable_secondary": float(p_one),
        "interpretation": (
            "Negative mean delta AUC indicates lower AUC and therefore "
            "greater vulnerability in spaced versus massed networks, "
            "averaged across spatial and cue tasks."
        ),
        "first_failures": failures[:20],
    }

    pd.DataFrame([{
        k: v for k, v in result.items() if k != "first_failures"
    }]).to_csv(
        OUTPUT_DIR / "omnibus_permutation_summary.csv",
        index=False
    )

    with open(
        OUTPUT_DIR / "omnibus_permutation_result.json",
        "w",
        encoding="utf-8"
    ) as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    print("\n============================================================")
    print(" OMNIBUS RESULT")
    print("============================================================")
    print(f"Observed mean delta AUC = {omnibus_obs:.6f}")
    print(f"Two-sided permutation p = {p_two:.6g}")
    print(f"One-sided p (spaced more vulnerable) = {p_one:.6g}")
    print(f"Valid permutations = {valid.sum()}/{N_PERMUTATIONS}")
    print(f"\nResults saved in:\n{OUTPUT_DIR.resolve()}")


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
