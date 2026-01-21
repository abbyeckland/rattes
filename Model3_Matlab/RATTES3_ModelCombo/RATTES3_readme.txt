RATTES3_readme
RATTES 3 combine the sediment yield and multiple linear regression model components, back projects the multiple linear regression model results for instances where the only available survey data post-date dam closure, and completes model analysis. 
-Figures 2 and 3 in the paper are created using this code
-Table 1 in the paper is created using this code
-Supplemental delta data table is created using this code
-Supplemental figures are created using this code. 
-The RATTESqueries function prints data queries which exist in the paper
TO RUN:
-Using output from RATTES1 (MATLAB code) and RATTES2 (R code) run RATTES3 
-Input files: add the prefix "mlr_" from RATTES2 output and use to run RATTES3. The prefix is added to differentiate with final model results:
 (mlr_cap_010626.csv, mlr_cap_up_95_010626,mlr_cap_low_95_010626.csv, mlr_sv_010626.csv, mlr_sv_up_95_010626,mlr_sv_low_95_010626.csv) 
-creates .mat file outputs and csv output
