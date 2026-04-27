#!/usr/bin/env python3

import logging
import coloredlogs
import argparse
import pandas as pd
import numpy as np
from utils import pls, bootstrap_test

def permutation_test(xm, ym, n_permutations=100000):
    """
    Esegue un test di permutazione per valutare la significatività della PLS.

    Args:
        xm: Matrice dei descrittori (attività cerebrale).
        ym: Matrice delle variabili dummy (appartenenza al gruppo).
        n_permutations: Numero di permutazioni da eseguire.

    Returns:
        p_value: P-value che indica la significatività della PLS.
    """
    # 1. Calcola la statistica di test sui dati reali
    u, s, v = pls(xm, ym)
    test_statistic_real = s[0]  # Ad esempio, usa la prima singular value come statistica di test

    # 2. Inizializza un array per memorizzare le statistiche di test delle permutazioni
    test_statistics_permutations = np.zeros(n_permutations)

    # 3. Esegui le permutazioni
    for i in range(n_permutations):
        # Mescola casualmente le etichette dei gruppi
        ym_permuted = ym[np.random.permutation(ym.shape[0])]

        # Esegui la PLS sui dati permutati
        u_permuted, s_permuted, v_permuted = pls(xm, ym_permuted)
        test_statistics_permutations[i] = s_permuted[0]  # Memorizza la statistica di test

    # 4. Calcola il p-value
    p_value = np.mean(test_statistics_permutations >= test_statistic_real)

    return p_value
def main():
    logger = logging.getLogger(__name__)
    logging.basicConfig(format='[%(funcName)s] - %(asctime)s - %(message)s', level=logging.INFO)
    coloredlogs.install(level='DEBUG', logger=logger)

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="input csv file", default='C:/Users/tomgo/Desktop/PLS_python/PLS/SC-HC/PLS_SC_HC.csv')
    parser.add_argument('-o', '--output', help="output base path", default='C:/Users/tomgo/Desktop/PLS_python/PLS/output')
    parser.add_argument('-b', '--bootstrap', help="number of bootstrap iterations", type=int, default=1000)
    parser.add_argument('-p', '--procrustes', help="enable Procrustes rotation in bootstrap", action='store_true', default=False)
    parser.add_argument('-c', '--columns', help="number of descriptor columns", type=int, default=2)

    # Aggiunta dell'argomento n_permutations
    parser.add_argument('-n', '--n_permutations', help="number of permutations for permutation test", type=int,
                        default=100000)

    args = parser.parse_args()

    logger.info('reading data...')
    data = pd.read_csv(args.input)  # csv file organized according to the example
    if 'Unnamed: 0' in data.columns:
        data.drop('Unnamed: 0', axis=1, inplace=True)  # remove index column if present


    # create dummy matrix containing the experimental group of each sample
    dum = pd.get_dummies(data, columns=['group'])
    n = data['group'].nunique()

    x = dum.iloc[:, args.columns:(data.shape[1] - 1)]
    y = dum.iloc[:, -n:]
    xm = x.to_numpy()
    ym = y.to_numpy()

    logger.info('computing PLS...')
    u, s, v = pls(xm, ym)

    logger.info('PLS computed, now performing bootstrap...')
    vs = bootstrap_test(xm, ym, v, u, args.bootstrap, args.procrustes)

    # Ensure the correct order of columns starting from "FP"
    regions = data.columns[3:]
    logger.info(f"Shape of v: {v.shape}")
    logger.info(f"Shape of vs: {vs.shape}")
    logger.info(f"Number of regions: {len(regions)}")

    # Calcola il p-value con il test di permutazione
    logger.info('Performing permutation test...')
    p_value = permutation_test(xm, ym, args.n_permutations)
    logger.info(f"P-value del test di permutazione: {p_value}")

    # Valuta la significatività
    alpha = 0.05  # Livello di significatività
    if p_value < alpha:
        logger.info("La differenza tra i gruppi è significativa.")
    else:
        logger.info("La differenza tra i gruppi non è significativa.")

    vpd = pd.DataFrame((v / vs), columns=regions)
    # Save the DataFrame without indices
    logger.info('saving output data...')
    vpd.to_csv(args.output + '_saliences.csv', index=False, header=True)

    upd = pd.DataFrame(u, columns=y.columns)
    upd.to_csv(args.output + '_contrasts.csv', index=False, header=True)

if __name__ == "__main__":
    main()

import pandas as pd
import matplotlib.pyplot as plt

# Define the colors for the groups
group_colors = {
    'HC': 'dimgrey',  # Dark gray for HC
    'SC': 'darkorange'  # Dark blue for SS
}

# Define the group mapping and color palette
group_mapping = {
    "Cortex": ["FP", "M1", "M2", "S1", "SS", "GA", "VA", "AuA", "ViA", "dACC", "Vacc",
               "PL", "IL", "OA", "AU", "RS", "AA", "TeAA", "Prh", "Ech","lpEC"],
    "Olfactory": ["aOn", "TT", "DPA", "PC",  "PAA"],
    "HPC": ["CA1", "CA2", "CA3", "DG", "FC", "IG", "Sub", "proSub"],
    "Amygdala": ["Clst", "dEPn", "vEPn", "LAn", "CCA","BLA", "BMe", "PAn"],
    "BG": ["DLS", "DMS", "NACC SHELL", "NACC CORE", "FStr", "OT", "LSN", "SFn",
           "aAA", "CA", "IAn", "MAn", "Gpe", "Gpi", "SubInn", "MGCn", "MSn", "DBn",
           "TnS", "BNST"],
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
                     "STN", "TN ", "ZI ", "FF"],
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
    bars = ax.bar(data.columns, data.iloc[row],
                  color=[group_colors['HC'] if col == 'group_HC' else group_colors['SC'] for col in data.columns],
                  alpha=0.9, linewidth=0.7, edgecolor='black', width=0.8)
    ax.set_title(title, fontname='Century Gothic', fontsize=14, fontweight='bold')
    ax.set_xlabel('Regions', fontname='Century Gothic', fontsize=13, fontweight='bold')
    ax.set_ylabel('Contrast', fontname='Century Gothic', fontsize=13, fontweight='bold')
    ax.set_xticks(data.columns)
    ax.set_xticklabels(data.columns, rotation=45, fontname='Century Gothic', fontsize=10, fontweight='bold')
    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.legend(bars[:2], ['HC', 'SC'], loc='upper right', bbox_to_anchor=(1.5, 1.5))
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
saliences_files = ['output_saliences_1.csv', 'output_saliences_2.csv']
contrasts_files = ['output_contrasts_1.csv', 'output_contrasts_2.csv']

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