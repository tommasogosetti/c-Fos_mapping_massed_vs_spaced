Each experimental subject has a dedicated folder. This folder contains: 

"png" subfolder: contains the images of the slices in .png format on which the analyses were conducted. These images are much smaller than the original ones.

"tiff" subfolder: contains the original images in .tiff format, with a scale of 0.44 microns/pixel, as exported from the microscope .czi propetary format scans. 

"scales.txt" file: correspondence between tiff-png and micron/pixel scale for each png.

"RefAtlasRegions" subfolder: results in terms of number of objects and area (in pixels) per brain region (one file is present for each slice). 

"MASK_RefAtlasRegions" subfolder: same, but data obtained with a mask useful for dividing nacc-shell from nacc-core and dls from dms 

"1_add_scale+area+density_to_csv.py" script: script which creates new files by adding scale, area, and density, obtained from the densities in the scales.txt file. 

"RefAtlas_with_density" subfolder: folder with the files produced by the script mentioned above 

"2_unite_fettine.py" script: script that unites the files of the various slices in "RefAtlas_with_density" into a single file, representative of the whole mouse, containing sums and averages of the values of the individual slices as well as the standard deviation of the density of an area among the various slices, as well as the ideal count of slices for that area.

"FXX_unito.csv" file: output of the script mentioned above

"count_ideali.txt" file: file from which the script mentioned above takes the list of ideal counts for the various areas 

"3_screma_sottolayer.py" script: starting from FXX_unito.csv, this script unites the areas that we decided to cluster together and eliminates various other areas deemed uninteresting 

"FXX_scremato.csv" file: file of the mouse with the clustered and selectedareas, output of the script mentioned above

"id_aree_da_considerare.txt" file: file from which the script mentioned above takes the list of areas to include 

"sottolayer_da_scremare.txt" file: file from which the script mentioned above takes the list of areas to cluster

"4_MASK_add_scale+area+density_to_csv" and "5_MASK_unite_fettine.py" scripts: these do the same as the homonymous non-MASK scripts but starting from the MASK_RefAtlasRegions folder instead of the RefAtlasRegions folder. 

"MASK_RefAtlas_with_density" folder and "MASK_FXX_unito.csv" file: output of the two MASK scripts mentioned above 

"6_add_mask_to_final_file.py" script: script that obtains a file by adding dls and dms and nacc-core and nacc-shell, obtained from MASK_FXX_unito.csv, to FXX_scremato.csv 

"FXX_final_for_variability" file: output of the script mentioned above, as well as the final output of the workflow described here, which we used to perform an initial analysis of variability. This file contains 274 areas (including nac-shell, nac-core, dls, and dms), with the area (in microns) of the areas, object count of the areas, density (in microns/pixel) of the areas, standard deviation, and slice count.