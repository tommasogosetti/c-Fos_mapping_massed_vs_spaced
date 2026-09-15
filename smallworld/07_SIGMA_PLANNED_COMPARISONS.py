#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
07_SIGMA_PLANNED_COMPARISONS.py
===============================

Purpose
-------
Formal group-level inference for the small-world coefficient sigma in the
2 x 2 experimental design:

    Spatial: SS vs MS
    Cue:     SC vs MC
    Spaced:  SS vs SC
    Massed:  MS vs MC

The script deliberately does NOT run an ordinary ANOVA on four group-level
sigma values and does NOT treat jackknife estimates as independent animals.
Instead, each inferential comparison is a SUBJECT-LABEL PERMUTATION TEST:
subjects from the two groups being compared are pooled, labels are reassigned
without replacement while preserving group sizes (8 vs 8), both networks are
rebuilt from the subject-level data, and sigma is recalculated.

Outputs also include leave-one-subject-out jackknife SE values for plotting
uncertainty around each observed group-level sigma.

IMPORTANT NETWORK RULE
----------------------
This script imports the already validated network-construction engine from:
    topology_smallworld_permutation_PYCHARM_FAST_RESUME.py

That engine reconstructs the manuscript networks using:
    positive Pearson correlations
    r > 0.65
    conventional TWO-SIDED p < 0.04
    topology on the unweighted giant component

The observed SS/MS/SC/MC networks must reproduce exactly:
    SS: 374 edges, GC 121
    MS: 237 edges, GC 114
    SC: 253 edges, GC 121
    MC: 430 edges, GC 113

Place this script in the SAME folder as:
    topology_smallworld_permutation_PYCHARM_FAST_RESUME.py
    ss net.csv
    ms net.csv
    sc net.csv
    mc net.csv

Primary outputs
---------------
Folder: sigma_group_comparisons_results/

    observed_sigma_high_precision.csv
    sigma_jackknife_summary.csv
    sigma_jackknife_all.csv
    sigma_pairwise_permutation_results.csv
    sigma_pairwise_permutation_results_holm.csv
    permutation_<comparison>.csv              [checkpoint/distribution]
    analysis_settings.txt

Run modes
---------
Set RUN_MODE below:
    TEST  = very short smoke test
    FINAL = manuscript analysis

FINAL defaults:
    1000 null graphs per observed network for displayed sigma
      50 null graphs per leave-one-out network for jackknife SE
      20 null graphs per reconstructed network in permutation tests
    5000 permutations per planned comparison

The permutation p-value is two-sided and uses the Monte Carlo +1 correction.
"""

from __future__ import annotations

from pathlib import Path
import importlib.util
import math
import os
import time

import numpy as np
import pandas as pd


# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent
SOURCE_SCRIPT = BASE_DIR / "topology_smallworld_permutation_PYCHARM_FAST_RESUME.py"
OUT_DIR = BASE_DIR / "sigma_group_comparisons_results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FILES = {
    "SS": BASE_DIR / "ss net.csv",
    "MS": BASE_DIR / "ms net.csv",
    "SC": BASE_DIR / "sc net.csv",
    "MC": BASE_DIR / "mc net.csv",
}

# Four biologically interpretable simple comparisons.  If you want only the
# protocol comparisons requested most directly by the paper, use:
# PLANNED_COMPARISONS = [("SS", "MS"), ("SC", "MC")]
PLANNED_COMPARISONS = [
    ("SS", "MS"),  # protocol effect within spatial task
    ("SC", "MC"),  # protocol effect within cue task
    ("SS", "SC"),  # task effect within spaced training
    ("MS", "MC"),  # task effect within massed training
]

RUN_MODE = "FINAL"  # "TEST" or "FINAL"
SEED = 20260913
SWAP_MULTIPLIER = 10
CHECKPOINT_EVERY = 25
RESUME = True

MODE_SETTINGS = {
    "TEST": {
        "n_perm": 20,
        "n_null_perm": 2,
        "n_null_observed": 10,
        "n_null_jackknife": 2,
    },
    "FINAL": {
        "n_perm": 5000,
        "n_null_perm": 20,
        "n_null_observed": 1000,
        "n_null_jackknife": 50,
    },
}

# Optional environment overrides, useful for sensitivity runs without editing.
# Examples (PowerShell):
#   $env:N_PERM="1000"
#   $env:N_NULL_PERM="50"
#   py 07_SIGMA_PLANNED_COMPARISONS.py


# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------

def load_source():
    if not SOURCE_SCRIPT.exists():
        raise FileNotFoundError(
            f"Missing {SOURCE_SCRIPT.name}. Put this script in the same folder "
            "as the validated permutation engine."
        )
    spec = importlib.util.spec_from_file_location("network_engine", SOURCE_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def holm_adjust(p_values: dict[str, float]) -> dict[str, float]:
    """Holm step-down family-wise correction."""
    valid = [(k, float(v)) for k, v in p_values.items() if np.isfinite(v)]
    valid.sort(key=lambda kv: kv[1])
    m = len(valid)
    out = {k: np.nan for k in p_values}
    running = 0.0
    for rank, (name, p) in enumerate(valid, start=1):
        adj = min(1.0, (m - rank + 1) * p)
        running = max(running, adj)
        out[name] = running
    return out


def jackknife_se(values: list[float]) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 2:
        return np.nan
    m = float(np.mean(x))
    return float(np.sqrt((n - 1) / n * np.sum((x - m) ** 2)))


def group_seed(group: str, extra: int = 0) -> int:
    offsets = {"SS": 11_000_000, "MS": 22_000_000,
               "SC": 33_000_000, "MC": 44_000_000}
    return int((SEED + offsets[group] + extra) % (2**32 - 1))


def safe_sigma(engine, df: pd.DataFrame, name: str, n_null: int, seed: int):
    rng = np.random.default_rng(seed)
    result = engine.calculate_metrics(
        df,
        name,
        n_random=n_null,
        swap_multiplier=SWAP_MULTIPLIER,
        rng=rng,
    )
    if result["random_networks_valid"] < n_null:
        raise RuntimeError(
            f"{name}: only {result['random_networks_valid']}/{n_null} valid null "
            "networks were generated."
        )
    return result


# -----------------------------------------------------------------------------
# OBSERVED + VALIDATION
# -----------------------------------------------------------------------------

def calculate_observed(engine, groups, n_null_observed):
    print("\n[1/3] High-precision observed sigma values")
    observed = {}
    rows = []
    for g in ["SS", "MS", "SC", "MC"]:
        t0 = time.time()
        res = safe_sigma(
            engine, groups[g], g,
            n_null=n_null_observed,
            seed=group_seed(g, 101),
        )
        observed[g] = res
        rows.append(res)
        print(
            f"  {g}: sigma={res['small_world_sigma']:.6f}, "
            f"edges={res['edges']}, GC={res['gc_nodes']}, "
            f"nulls={res['random_networks_valid']} "
            f"({time.time()-t0:.1f}s)"
        )

    # Strict manuscript-network validation.
    engine.OUTPUT_DIR = OUT_DIR
    engine.validate_observed(observed, stop_on_fail=True)

    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "observed_sigma_high_precision.csv", index=False)
    return observed, out


# -----------------------------------------------------------------------------
# JACKKNIFE UNCERTAINTY FOR FIGURE ERROR BARS
# -----------------------------------------------------------------------------

def calculate_jackknife(engine, groups, n_null_jackknife):
    print("\n[2/3] Leave-one-subject-out sigma jackknife")
    all_rows = []
    summary_rows = []

    for g in ["SS", "MS", "SC", "MC"]:
        vals = []
        df = groups[g].reset_index(drop=True)
        for i in range(len(df)):
            loo = df.drop(index=i).reset_index(drop=True)
            omitted = str(df.loc[i, "animal"])
            res = safe_sigma(
                engine,
                loo,
                f"{g}_omit_{omitted}",
                n_null=n_null_jackknife,
                seed=group_seed(g, 10_000 + i),
            )
            sigma = float(res["small_world_sigma"])
            vals.append(sigma)
            all_rows.append({
                "group": g,
                "omitted_subject": omitted,
                "sigma": sigma,
                "edges": res["edges"],
                "gc_nodes": res["gc_nodes"],
                "valid_nulls": res["random_networks_valid"],
            })

        se = jackknife_se(vals)
        summary_rows.append({
            "group": g,
            "n_leave_one_out": len(vals),
            "jackknife_mean_sigma": float(np.mean(vals)),
            "jackknife_SE_sigma": se,
            "leave_one_out_min": float(np.min(vals)),
            "leave_one_out_max": float(np.max(vals)),
        })
        print(f"  {g}: jackknife SE={se:.6f}")

    pd.DataFrame(all_rows).to_csv(
        OUT_DIR / "sigma_jackknife_all.csv", index=False
    )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "sigma_jackknife_summary.csv", index=False)
    return summary


# -----------------------------------------------------------------------------
# PAIRWISE SUBJECT-LABEL PERMUTATION TESTS
# -----------------------------------------------------------------------------

def observed_pair_delta_same_estimator(
    engine, group_a, group_b, df_a, df_b, n_null_perm, comparison_index
):
    # Use the same low-null estimator used inside permutations, so the observed
    # test statistic and null distribution have the same Monte Carlo precision.
    ra = safe_sigma(
        engine, df_a, f"{group_a}_obs_test",
        n_null=n_null_perm,
        seed=SEED + 700_000 + comparison_index * 10_000 + 1,
    )
    rb = safe_sigma(
        engine, df_b, f"{group_b}_obs_test",
        n_null=n_null_perm,
        seed=SEED + 700_000 + comparison_index * 10_000 + 2,
    )
    return (
        float(ra["small_world_sigma"]),
        float(rb["small_world_sigma"]),
        float(ra["small_world_sigma"] - rb["small_world_sigma"]),
    )


def run_pairwise_permutation(
    engine,
    group_a: str,
    group_b: str,
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    n_perm: int,
    n_null_perm: int,
    comparison_index: int,
):
    label = f"{group_a}_vs_{group_b}"
    checkpoint = OUT_DIR / f"permutation_{label}.csv"

    obs_a_test, obs_b_test, observed_delta = observed_pair_delta_same_estimator(
        engine, group_a, group_b, df_a, df_b, n_null_perm, comparison_index
    )

    pool = pd.concat([df_a, df_b], ignore_index=True)
    n_a = len(df_a)
    n_b = len(df_b)
    idx = np.arange(len(pool))

    records = []
    start = 0
    if RESUME and checkpoint.exists():
        try:
            old = pd.read_csv(checkpoint)
            old = old[old["permutation"] <= n_perm]
            old = old.sort_values("permutation").drop_duplicates("permutation")
            if len(old):
                records = old.to_dict("records")
                start = int(old["permutation"].max())
                print(f"    RESUME {label}: {start}/{n_perm} already complete")
        except Exception as exc:
            print(f"    WARNING: checkpoint unreadable ({exc}); starting fresh")
            records = []
            start = 0

    print(
        f"  {label}: observed test-estimator sigma "
        f"{group_a}={obs_a_test:.4f}, {group_b}={obs_b_test:.4f}, "
        f"delta={observed_delta:.4f}"
    )

    for b in range(start, n_perm):
        local_seed = int(
            (SEED + comparison_index * 100_000_000 + 1_000_003 * (b + 1))
            % (2**32 - 1)
        )
        rng = np.random.default_rng(local_seed)
        perm = rng.permutation(idx)
        pa = pool.iloc[perm[:n_a]].reset_index(drop=True)
        pb = pool.iloc[perm[n_a:n_a+n_b]].reset_index(drop=True)

        # Same RNG stream for the two reconstructed networks in this permutation.
        ra = engine.calculate_metrics(
            pa, f"{group_a}_perm_{b+1}",
            n_random=n_null_perm,
            swap_multiplier=SWAP_MULTIPLIER,
            rng=rng,
        )
        rb = engine.calculate_metrics(
            pb, f"{group_b}_perm_{b+1}",
            n_random=n_null_perm,
            swap_multiplier=SWAP_MULTIPLIER,
            rng=rng,
        )

        if (
            ra["random_networks_valid"] < n_null_perm
            or rb["random_networks_valid"] < n_null_perm
            or not np.isfinite(ra["small_world_sigma"])
            or not np.isfinite(rb["small_world_sigma"])
        ):
            continue

        delta = float(ra["small_world_sigma"] - rb["small_world_sigma"])
        records.append({
            "permutation": b + 1,
            f"sigma_{group_a}": float(ra["small_world_sigma"]),
            f"sigma_{group_b}": float(rb["small_world_sigma"]),
            "delta_A_minus_B": delta,
        })

        done = b + 1
        if done % CHECKPOINT_EVERY == 0 or done == n_perm:
            pd.DataFrame(records).sort_values("permutation").drop_duplicates(
                "permutation"
            ).to_csv(checkpoint, index=False)
            print(f"    {label}: {done}/{n_perm} [checkpoint]")

    dist = pd.DataFrame(records)
    vals = pd.to_numeric(
        dist.get("delta_A_minus_B", pd.Series(dtype=float)), errors="coerce"
    ).to_numpy(dtype=float)
    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        p_two = np.nan
    else:
        p_two = float(
            (1 + np.sum(np.abs(vals) >= abs(observed_delta)))
            / (1 + len(vals))
        )

    return {
        "comparison": label,
        "group_A": group_a,
        "group_B": group_b,
        "observed_sigma_A_test_estimator": obs_a_test,
        "observed_sigma_B_test_estimator": obs_b_test,
        "observed_delta_A_minus_B": observed_delta,
        "n_valid_permutations": int(len(vals)),
        "n_null_per_reconstructed_network": int(n_null_perm),
        "p_two_sided": p_two,
    }


def run_all_pairwise(engine, groups, n_perm, n_null_perm):
    print("\n[3/3] Planned pairwise subject-label permutation tests")
    rows = []
    for k, (a, b) in enumerate(PLANNED_COMPARISONS, start=1):
        t0 = time.time()
        row = run_pairwise_permutation(
            engine,
            a, b,
            groups[a], groups[b],
            n_perm=n_perm,
            n_null_perm=n_null_perm,
            comparison_index=k,
        )
        rows.append(row)
        print(
            f"    -> {row['comparison']}: p={row['p_two_sided']:.6g} "
            f"({time.time()-t0:.1f}s)"
        )

    raw = pd.DataFrame(rows)
    raw.to_csv(OUT_DIR / "sigma_pairwise_permutation_results.csv", index=False)

    pmap = dict(zip(raw["comparison"], raw["p_two_sided"]))
    adj = holm_adjust(pmap)
    corrected = raw.copy()
    corrected["p_Holm_across_planned_sigma_comparisons"] = corrected[
        "comparison"
    ].map(adj)
    corrected.to_csv(
        OUT_DIR / "sigma_pairwise_permutation_results_holm.csv", index=False
    )
    return corrected


def save_settings(settings):
    lines = [
        "SIGMA PLANNED COMPARISONS ANALYSIS",
        "==================================",
        f"RUN_MODE = {RUN_MODE}",
        f"N permutations per comparison = {settings['n_perm']}",
        f"Nulls per permuted network = {settings['n_null_perm']}",
        f"Nulls per observed network = {settings['n_null_observed']}",
        f"Nulls per jackknife network = {settings['n_null_jackknife']}",
        f"Seed = {SEED}",
        "Network rule = positive Pearson r > 0.65 AND two-sided p < 0.04",
        "Topology = unweighted giant component",
        "Sigma = (C/C_rand)/(L/L_rand)",
        "Pairwise test = subject-label permutation without replacement",
        "P values = two-sided Monte Carlo with +1 correction",
        "Multiple comparison correction = Holm across planned sigma comparisons",
        "Planned comparisons = " + ", ".join(f"{a} vs {b}" for a, b in PLANNED_COMPARISONS),
    ]
    (OUT_DIR / "analysis_settings.txt").write_text("\n".join(lines), encoding="utf-8")


def main():
    mode = RUN_MODE.strip().upper()
    if mode not in MODE_SETTINGS:
        raise ValueError(f"RUN_MODE must be TEST or FINAL, not {RUN_MODE!r}")

    settings = dict(MODE_SETTINGS[mode])
    settings["n_perm"] = int(os.environ.get("N_PERM", settings["n_perm"]))
    settings["n_null_perm"] = int(os.environ.get("N_NULL_PERM", settings["n_null_perm"]))
    settings["n_null_observed"] = int(os.environ.get("N_NULL_OBS", settings["n_null_observed"]))
    settings["n_null_jackknife"] = int(os.environ.get("N_NULL_JACK", settings["n_null_jackknife"]))

    print("=" * 78)
    print(" SMALL-WORLD SIGMA: PLANNED GROUP COMPARISONS")
    print("=" * 78)
    print(f"Mode: {mode}")
    print(f"Output: {OUT_DIR}")
    print(f"Permutations/comparison: {settings['n_perm']}")
    print(f"Nulls/permuted network: {settings['n_null_perm']}")
    print(f"Nulls/observed network: {settings['n_null_observed']}")
    print(f"Nulls/jackknife network: {settings['n_null_jackknife']}")

    engine = load_source()
    groups = {g: engine.read_group(p) for g, p in FILES.items()}

    observed, _ = calculate_observed(
        engine, groups, settings["n_null_observed"]
    )
    jack = calculate_jackknife(
        engine, groups, settings["n_null_jackknife"]
    )
    pairwise = run_all_pairwise(
        engine, groups,
        settings["n_perm"], settings["n_null_perm"]
    )

    # Figure-friendly table: high-precision observed sigma + jackknife SE.
    fig_rows = []
    for g in ["SS", "MS", "SC", "MC"]:
        j = jack.loc[jack["group"] == g].iloc[0]
        fig_rows.append({
            "group": g,
            "observed_sigma": observed[g]["small_world_sigma"],
            "jackknife_SE": j["jackknife_SE_sigma"],
        })
    pd.DataFrame(fig_rows).to_csv(
        OUT_DIR / "figure4E_sigma_with_jackknife_SE.csv", index=False
    )

    save_settings(settings)

    print("\nDONE")
    print("Most useful files:")
    print("  figure4E_sigma_with_jackknife_SE.csv")
    print("  sigma_pairwise_permutation_results_holm.csv")


if __name__ == "__main__":
    main()
