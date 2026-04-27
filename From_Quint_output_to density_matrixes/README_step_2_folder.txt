“274_strutture” subfolder: collection of all the final files of the various mice brought to 274 structures, products of the workflow found in the folder “C-Fos whole brain mapping”

“aggiungi_sigle_togli_id.py” script: script that removes the column of the identification numbers of the various structures (kept intact until now from the original RefAtlasRegions files) and replaces it with another one containing the abbreviations that we have chosen to use for the various areas

“con_sigle” subfolder: output of the script above, containing all the mice files, each with the id column replaced by the abbreviations column

“abbreviazioni.txt” file: file used by the script above, which contains all the abbreviations that we use

“to_134.py” script: script that removes various areas and clusters others, bringing the files in the folder “con_sigle” to 134 structures

“134_strutture” folder: output of the script above, containing one file for each mouse with density etc. for 134 structures.

“aree_da_unire_eliminare_copiare_per_scendere_a_134.txt” file: file on which the script above is based to understand which areas to cluster, which to eliminate, and which to copy unchanged.
