# RATTES mlr model component
<u>Overview:</u> this code is for multiple linear regression model development and predictions for the second RATTES model component. This code reads in input files derived from ResNet (Hurst et al., 2025), dam parameter datasets (flow-conditioned parameter grids (FCPG) from Barnhart et al. (2020) and LakeCat from Hill et al. (2018)), and the sediment-contributing drainage area (A_sedMLR) timeseries, parameterized using the initial trap efficiency and design capacity for all ResNet reservoirs. The code pre-processes the data to ensure all requirements for multiple linear regression are met. Then, several MLR models are built using different combinations of parameter data. Finally, the final MLR to make predictions is selected and parameterized to make sedimentation rate predictions for all unsurveyed RATTES reservoirs for each year between the first year after the prediction year (*yrp*), which is the year with the first available capacity data, and the end of the time loop (year 2050).

<u>Scripts:</u> R scripts and Python Jupyter Notebooks are utilized in this analysis. The R script **MLR_toshare.R** is the main script, which contains comments describing when the other scripts are run in a semi-automated manner. The two additional R scripts (**readin_lakecat_data.R** and **compute_wSAatdam.R**) have already been run, and their output files are "lakecat_ss_data_DATE.csv" and "ss_wbcomid_data_DATE.csv" for the former R script, and "wSAatdam_predSites_notTransformed_DATE.csv" for the latter R script.

<u>Input files:</u> The main input files are located in the folder "SYfiles_MonYEAR", including: "MLR_InputMMDDYY.csv", which is derived from ResNet, and "MLR_SAatDamInitialMMDDYY.csv", which contains the sediment-contributing drainage area for each reservoir in ResNet calculated using design trap efficiencies and capacities.

<u>Steps prior to running:</u>
- Ensure working directories are up to date in **MLR_toShare.R**. "version" (Line 23) will add the specified date string to the end of all output files; this date determines whether or not the additional 2 R scripts are run and should match the "version" in both additional R scripts. Line 31 should also be updated if input file names have changed.
- If the input file names have changed in the SYfiles_MonYEAR folder, then "version" must be updated in each R script. Once the version name changes, the additional 2 R scripts will re-run automatically, prompted by the main R Script.

<u>Quick notes on csv files:</u>
- "lakecat_ss_data_MMDDYY.csv"
    - this is one of the results dataframes from **readin_lakecat_data.R**, joining lakecat data to SYfiles_Jan2026/MLRInput_122925.csv file (sediment yield model output) for ResNet reservoirs. If a new version string is used, a new file will be created with the updated date string.
- "ss_wbcomid_data_MMDDYY.csv"
    - this is the other results dataframe from **readin_lakecat_data.R**, joining NHD WBCOMID with ResNet SID and lakecat data columns. Does not include all other ResNet data columns. If a new version string is used, a new file will be created with the updated date string.
- "wSAatdam_predSites_notTransformed_MMDDYY.csv"
    - this is the results dataframe from **compute_wSAatdam.R**, which computes the time-weighted effective sediment contributing drainage area using design capacities and trap efficiencies, for each timestep (calendar year), or A_sedMLR. This dataframe is used as input to the prediction loop. If a new version string is used, a new file will be created with the updated date string.

<u>Other notes:</u>
- Upon running **MLR_toShare.R**, additional folders with figures, results dataframes, model summaries, etc will be automatically created. See the comments in the code for more details on these derivations.


