import logging
import coloredlogs
import argparse
import pandas as pd
import numpy as np
from utils import pls, bootstrap_test


def benjamini_hochberg_correction(pvals, alpha=0.05):
    m = len(pvals)
    sorted_indices = np.argsort(pvals)
    sorted_pvals = pvals[sorted_indices]
    thresholds = (np.arange(1, m + 1) / m) * alpha
    significant = np.zeros(m, dtype=bool)
    max_i = -1
    for i in range(m - 1, -1, -1):
        if sorted_pvals[i] <= thresholds[i]:
            max_i = i
            break
    if max_i >= 0:
        significant[sorted_indices[:max_i + 1]] = True
    return significant


def permutation_test_extended(xm, ym, n_permutations=1000):
    u_real, s_real, v_real = pls(xm, ym)
    n_lv = s_real.shape[0]
    s_real_vals = s_real[:n_lv]
    v_real_vals = v_real[:n_lv, :]

    s_perm = np.zeros((n_permutations, n_lv))
    v_perm = np.zeros((n_permutations, n_lv, v_real.shape[1]))

    for i in range(n_permutations):
        ym_permuted = ym[np.random.permutation(ym.shape[0])]
        u_perm, s_perm_i, v_perm_i = pls(xm, ym_permuted)
        s_perm[i, :len(s_perm_i)] = s_perm_i[:n_lv]
        v_perm[i, :, :] = v_perm_i[:n_lv, :]

    # Calcolo p-value per contrasti
    p_s = np.array([(np.sum(s_perm[:, i] > s_real_vals[i]) + 1) / (n_permutations + 1) for i in range(n_lv)])

    # Calcolo p-value raw per salienze
    p_v_raw = np.zeros_like(v_real_vals)
    for lv in range(n_lv):
        for region in range(v_real_vals.shape[1]):
            p_v_raw[lv, region] = (np.sum(v_perm[:, lv, region] > v_real_vals[lv, region]) + 1) / (n_permutations + 1)

    # Correzione FDR per ogni LV
    p_v_fdr = np.zeros_like(p_v_raw, dtype=bool)
    for lv in range(n_lv):
        p_v_fdr[lv, :] = benjamini_hochberg_correction(p_v_raw[lv, :])

    p_values_contrasts = {f'LV{i + 1}': p_s[i] for i in range(n_lv)}
    p_values_saliencies_fdr = {f'LV{i + 1}': p_v_fdr[i, :] for i in range(n_lv)}

    return p_values_contrasts, p_values_saliencies_fdr, s_real_vals, u_real, v_real, p_v_raw, n_lv


def main():
    logger = logging.getLogger(__name__)
    logging.basicConfig(format='[%(funcName)s] - %(asctime)s - %(message)s', level=logging.INFO)
    coloredlogs.install(level='DEBUG', logger=logger)

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="input csv file",
                        default='C:/Users/tomgo/Desktop/PLS_python/PLS/SS-SC-MS-MC/PLS_SS_SC_MS_MC.csv')
    parser.add_argument('-o', '--output', help="output base path",
                        default='C:/Users/tomgo/Desktop/PLS_python/PLS/output')
    parser.add_argument('-b', '--bootstrap', help="number of bootstrap iterations", type=int, default=1000)
    parser.add_argument('-p', '--procrustes', help="enable Procrustes rotation in bootstrap", action='store_true',
                        default=False)
    parser.add_argument('-c', '--columns', help="number of descriptor columns", type=int, default=2)
    parser.add_argument('-n', '--n_permutations', help="number of permutations for permutation test", type=int,
                        default=1000)

    args = parser.parse_args()

    logger.info('Reading data...')
    data = pd.read_csv(args.input)
    if 'Unnamed: 0' in data.columns:
        data.drop('Unnamed: 0', axis=1, inplace=True)

    dum = pd.get_dummies(data, columns=['group'])
    n = data['group'].nunique()

    x = dum.iloc[:, args.columns:(data.shape[1] - n)]
    y = dum.iloc[:, -n:]

    xm = x.to_numpy()
    ym = y.to_numpy()

    logger.info('Computing PLS...')
    u_real, s_real, v_real = pls(xm, ym)

    logger.info('PLS computed, now performing bootstrap...')
    vs = bootstrap_test(xm, ym, v_real, u_real, args.bootstrap, args.procrustes)

    regions = data.columns[args.columns:]
    logger.info(f"Original number of regions: {len(regions)}")
    logger.info(f"Shape of v (saliencies): {v_real.shape}")

    if len(regions) != v_real.shape[1]:
        logger.warning(f"Mismatch tra numero regioni ({len(regions)}) e dimensione saliencies ({v_real.shape[1]}).")
        if len(regions) > v_real.shape[1]:
            logger.info(f"Taglio regions da {len(regions)} a {v_real.shape[1]}")
            regions = regions[:v_real.shape[1]]
        else:
            logger.error("Numero regioni inferiore a dimensione saliencies, verificare dati!")
            raise ValueError("Mismatch tra regions e saliencies: regions troppo poche.")

    logger.info('Performing extended permutation test for contrasts and saliencies...')
    (p_contrasts, p_saliencies_fdr, s_real_vals, u_real, v_real, p_v_raw, n_lv) = permutation_test_extended(xm, ym,
                                                                                                            args.n_permutations)

    # Salvataggio contrasti con p-value
    contrast_df = pd.DataFrame({
        'LV': [f'LV{i + 1}' for i in range(n_lv)],
        'SingularValue': s_real_vals,
        'p_value': [p_contrasts[f'LV{i + 1}'] for i in range(n_lv)]
    })
    contrast_df.to_csv(args.output + '_contrasts_permutation.csv', index=False)

    # Salvataggio loadings contrasti
    lv_columns = [f'LV{i + 1}' for i in range(n_lv)]
    u_df = pd.DataFrame(u_real, columns=lv_columns)
    u_df.to_csv(args.output + '_contrasts_loadings.csv', index=False)

    # Salvataggio salienze per ogni LV con p-value raw e significatività FDR
    for lv in range(n_lv):
        saliency_df = pd.DataFrame({
            'Region': regions,
            'Saliency': v_real[lv, :],
            'p_value_raw': p_v_raw[lv, :],
            'Significant_FDR': p_saliencies_fdr[f'LV{lv + 1}']
        })
        saliency_df.to_csv(f"{args.output}_saliencies_LV{lv + 1}.csv", index=False)

    logger.info('Permutation test and FDR correction completed, results saved.')

    # Salvataggio bootstrap normalized saliencies
    vpd = pd.DataFrame((v_real / vs).T, index=regions, columns=lv_columns)
    vpd.to_csv(args.output + '_saliences.csv', index=True, header=True)

    # Salvataggio loadings contrasti
    upd = pd.DataFrame(u_real, columns=lv_columns)
    upd.to_csv(args.output + '_contrasts.csv', index=False, header=True)


if __name__ == "__main__":
    main()

import pandas as pd
import matplotlib.pyplot as plt

# Define the colors for the groups
group_colors = {
    'SS': '#03055B',  # Dark gray for HC
    'SC': 'darkorange',  # Dark blue for SS
    'MS': '#8B0A1A',
    'MC': 'darkgreen'
}

# Define the group mapping and color palette
group_mapping = {
    "Cortex": ["FP", "M1", "M2", "S1", "SS", "GA", "VA", "AuA", "ViA", "dACC", "Vacc",
               "PL", "IL", "OA", "AU", "RS", "AA", "TeAA", "Prh", "Ech","lpEC"],
    "Olfactory": ["aOn", "TT", "DPA"],
    "HPC": ["CA1", "CA2", "CA3", "DG", "FC", "IG", "Sub", "proSub"],
    "Amygdala": ["Clst", "dEPn", "vEPn","PC","PAA", "LAn", "CCA","BLA", "BMe", "PAn","aAA", "CA", "IAn", "MAn"],
    "BG": ["DLS", "DMS", "NACC SHELL", "NACC CORE", "FStr", "OT", "LSN", "SFn",
            "Gpe", "Gpi", "SubInn", "MGCn", "MSn", "DBn",
           "TnS", "BNST","STN"],
    "Thalamus": ["VALCT", "VMT", "VPLT", "VPMT", "PTTn", "MGC", "dLGC",
                 "LPT", "PCT", "PlnT", "EnT", "AVT", "AMn", "IADT",
                 "IMDT", "MDT", "SMT", "pRE ", "PT ", "PVT", "RE",
                 "XiT", "Rh", "CMT", "PCn", "CLT", "RT",
                 "IGLLGC", "vlLGC", "subGn", "MHb", "LHb"],
    "Hypothalamus": ["SON", "PVH", "IPVH", "ArH", "ADPO",
                     "AVPO", "AVPV", "DMH", "MePO",
                     "MPA", "PSN", "MPO", "SPZ", "SCN", "AHN", "SMN",
                     "VMHn", "PHn",
                     "LHA", "LPO",
                     "PSTN",
                      "TN ", "ZI ", "FF"],
    "Midbrain": ["SCO", "reSN",
                 "VTA",
                 "PAG",
                 "PreCN",
                 "APTN",
                 "cSN"]
}

color_palette = {
    "Cortex": "lime",  # Green for cortex
    "Olfactory": "darkcyan",  # Light blue for olfactory
    "HPC": "forestgreen",  # Dark green for HPC
    "Amygdala": "darkkhaki",  # Gold for amygdala
    "BG": "navy",  # Navy for BG
    "Thalamus": "red",  # Coral for thalamus
    "Hypothalamus": "darkorange",  # Orange for hypothalamus
    "Midbrain": "blueviolet"  # BlueViolet for midbrain
}

# Create a mapping for colors based on regions
region_colors = {}
for group, regions in group_mapping.items():
    for region in regions:
        region_colors[region] = color_palette[group]


# Plotting function for contrasts
def plot_contrast_histogram(data, title, row, output_file):
    fig, ax = plt.subplots(figsize=(5, 7))  # Set figure size

    # Assign colors based on the group
    colors = []
    for col in data.columns:
        if col == 'group_SS':
            colors.append(group_colors['SS'])
        elif col == 'group_SC':
            colors.append(group_colors['SC'])
        elif col == 'group_MS':
            colors.append(group_colors['MS'])
        elif col == 'group_MC':
            colors.append(group_colors['MC'])
        else:
            colors.append('gray')  # Default color for other columns

    bars = ax.bar(data.columns, data.iloc[row], color=colors, alpha=0.9, linewidth=0.7, edgecolor='black', width=0.8)

    ax.set_title(title, fontname='Century Gothic', fontsize=14, fontweight='bold')
    ax.set_xlabel('Regions', fontname='Century Gothic', fontsize=13, fontweight='bold')
    ax.set_ylabel('Contrast', fontname='Century Gothic', fontsize=13, fontweight='bold')
    ax.set_xticks(data.columns)
    ax.set_xticklabels(data.columns, rotation=45, fontname='Century Gothic', fontsize=10, fontweight='bold')
    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.legend(bars[:4], ['MC','MS','SC','SS'], loc='upper right', bbox_to_anchor=(1.5, 1.5))

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.show()

# Plotting function for saliences
def plot_saliences_histogram(data, title, row, output_file):
    fig, ax = plt.subplots(figsize=(20, 12))  # Set figure size
    bar_width = 0.8  # Set the width of the bars (less than 1 for spacing)
    x_positions = np.arange(len(data.columns))  # Get the x positions for the bars
    bars = ax.bar(x_positions, data.iloc[row], color=[region_colors[col] for col in data.columns],
                  alpha=0.9, linewidth=0.7, edgecolor='black', width=bar_width)

    ax.set_title(title, fontname='Century Gothic', fontsize=14)
    ax.set_xlabel('Regions', fontname='Century Gothic', fontsize=13, fontweight='bold')
    ax.set_ylabel('Salience', fontname='Century Gothic', fontsize=13, fontweight='bold')
    ax.set_xticks(x_positions)  # Set the x ticks to the positions of the bars
    ax.set_xticklabels(data.columns, rotation=90, fontname='Century Gothic', fontsize=6, fontweight='bold')
    ax.axhline(y=3, color='black', linestyle='--', linewidth=0.8, label='Threshold: ±3')
    ax.axhline(y=-3, color='black', linestyle='--', linewidth=0.8)
    ax.axhline(y=0, color='black', linewidth=0.8)

    # Create a custom legend
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_palette[group]) for group in color_palette]
    ax.legend(handles, color_palette.keys(), loc='upper right', bbox_to_anchor=(1.25, 1.25))

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.show()


# Load and plot data from each file
saliences_files = ['output_saliences.csv']
contrasts_files = ['output_contrasts.csv']

# Plot saliences histograms
for file in saliences_files:
    vpd = pd.read_csv(file)
    for i in range(vpd.shape[0]):
        output_file = file.replace('.csv', f'_saliences_{i + 1}.png')
        plot_saliences_histogram(vpd, f'Saliences Histogram {i + 1} from {file}', i, output_file)

# Plot contrasts histograms
for file in contrasts_files:
    upd = pd.read_csv(file)
    for i in range(upd.shape[0]):
        output_file = file.replace('.csv', f'_contrasts_{i + 1}.png')
        plot_contrast_histogram(upd, f'Contrasts Histogram {i + 1} from {file}', i, output_file)