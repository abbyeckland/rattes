# update version of code to save to files:
version = "010626"   # udpate version date here

# import packages
library(dplyr)

### IMPORT DATA
# set wd to use in relative paths
setwd(".")  # make sure this matches working directory in MLR code

# import study site data for model building
ss_data <- read.csv("SYfiles_Jan2026/MLR_Input122925.csv", header=TRUE, stringsAsFactors = FALSE)  # update with current version

# Read in LakeCat data and merge to ss_data dataframe by column 'WBCOMID'
# subset ss_data by WBCOMID and SID only
column_names <- c("SID", "WBCOMID")
# subset the data frame by column names
ss_wbcomid <- ss_data[, column_names, drop = FALSE]

# drop columns that exist in every lakecat csv to ignore reading in duplicate columns
drop_columns <- c("CatAreaSqKm", "WsAreaSqKm", "inStreamCat")

# define the folder path containing the lakecat CSV files
lakecat_path <- "lakecat_data"

# list of all CSV files in the lakecat input data folder
files <- list.files(path = lakecat_path, pattern = "\\.csv$", full.names = TRUE)

# list of lakecat data columns to keep (only read-in relevant data to cut down on number of parameters)
columns <- c(
  'COMID', 'CatAreaSqKm', 'WsAreaSqKm', 'inStreamCat', # All tables
  'PctAg2006Slp20Ws', 'PctAg2006Slp10Ws', # AgMidHiSlopes
  'CBNFWs', 'FertWs', 'ManureWs', # Agricultural Nitrogen
  'BFIWs', # BFI (baseflow index, ratio of base flow to total flow as a %)
  'CanalDensWs', # CanalsDitches
  'CoalMineDensWs', # CoalMines
  'DamDensWs', 'DamNIDStorWs', 'DamNrmStorWs', # Dams
  'ElevWs', # Elevation
  'NPDESDensWs', 'SuperfundDensWs', 'TRIDensWs', # EPA_FRS (pollution)
  'PctFire2000Ws', 'PctFire2001Ws', 'PctFire2002Ws', 'PctFire2003Ws', 'PctFire2004Ws', 'PctFire2005Ws', 'PctFire2006Ws', 'PctFire2007Ws', 'PctFire2008Ws', 'PctFire2009Ws', 'PctFire2010Ws', # FirePerimeters
  'PctFrstLoss2001Ws', 'PctFrstLoss2002Ws', 'PctFrstLoss2003Ws', 'PctFrstLoss2004Ws', 'PctFrstLoss2005Ws', 'PctFrstLoss2006Ws', 'PctFrstLoss2007Ws', 'PctFrstLoss2008Ws', 'PctFrstLoss2009Ws', 'PctFrstLoss2010Ws', 'PctFrstLoss2011Ws', 'PctFrstLoss2012Ws', 'PctFrstLoss2013Ws', # Forest loss by year
  'HydrlCondWs', # GeoChemPhys3 (hydraulic conductivity of lithology)
  'CompStrgthWs', # GeoChemPhys4 (compressive strength of lithology)
  'PctImp2006Ws', 'PctImp2011Ws', # ImperviousSurfaces, 2006 and 2011
  'AgKffactWs', 'KffactWs', # Kffact (on agricultural lands and within watershed. From USLE)
  'PctCarbResidWs', 'PctNonCarbResidWs', 'PctAlkIntruVolWs', 'PctSilicicWs', 'PctExtruVolWs', 'PctColluvSedWs', 'PctGlacTilClayWs', 'PctGlacTilLoamWs', 'PctGlacTilCrsWs', 'PctGlacLakeCrsWs', 'PctGlacLakeFineWs', 'PctHydricWs', 'PctEolCrsWs', 'PctEolFineWs', 'PctSalLakeWs', 'PctAlluvCoastWs', 'PctCoastCrsWs', 'PctWaterWs', # Lithology
  'MineDensWs', # Mines
  'NABD_DensWs', 'NABD_NIDStorWs', 'NABD_NrmStorWs', # NABD (dam density, volume for georeferenced dams)
  'PctOw2006Ws', 'PctIce2006Ws', 'PctUrbOp2006Ws', 'PctUrbLo2006Ws', 'PctUrbMd2006Ws', 'PctUrbHi2006Ws', 'PctBl2006Ws', 'PctDecid2006Ws', 'PctConif2006Ws', 'PctMxFst2006Ws', 'PctShrb2006Ws', 'PctGrs2006Ws', 'PctHay2006Ws', 'PctCrop2006Ws', 'PctWdWet2006Ws', 'PctHbWet2006Ws', # NLCD2006 (land cover)
  'PctOw2011Ws', 'PctIce2011Ws', 'PctUrbOp2011Ws', 'PctUrbLo2011Ws', 'PctUrbMd2011Ws', 'PctUrbHi2011Ws', 'PctBl2011Ws', 'PctDecid2011Ws', 'PctConif2011Ws', 'PctMxFst2011Ws', 'PctShrb2011Ws', 'PctGrs2011Ws', 'PctHay2011Ws', 'PctCrop2011Ws', 'PctWdWet2011Ws', 'PctHbWet2011Ws', # NLCD2011 (land cover)
  'PctNonAgIntrodManagVegWs', # NonAgIntrodManagVeg (introduced, non agricultural invasive species)
  'Pestic97Ws', # Pesticides97
  'Precip8110Ws', 'Tmax8110Ws', 'Tmean8110Ws', 'Tmin8110Ws', # PRISM_1981_2010
  'RdDensWs', # RoadDensity
  'RdCrsWs', # RoadStreamCrossings
  'RunoffWs', # Runoff
  'ClayWs', 'SandWs', # STATSGO_Set1
  'OmWs', 'PermWs', 'RckdepWs', 'WtDepWs', # STATSGO_Set2
  'HUDen2010Ws', 'PopDen2010Ws', # USCensus2010
  'WetIndexWs', # WetIndx
  "ElevCat"
)

# loop through each CSV file in the folder, merging data with ss_wbcomid dataframe
for (file in files) {
  # Read the CSV file into a DataFrame
  df <- read.csv(file)
  df$COMID <- as.integer(df$COMID)
  
  # Step 1: Keep only relevant columns that exist in both `columns` and the dataframe
  relevant_columns <- intersect(columns, names(df))
  df <- df[, relevant_columns, drop = FALSE]
  
  # Step 2: Drop repetitive columns if they already exist in both `ss_wbcomid` and `df`
  for (col in drop_columns) {
    if (col %in% names(ss_wbcomid) && col %in% names(df)) {
      df <- df[, !(names(df) %in% col), drop = FALSE]
    }
  }
  
  # Step 3: Merge with ss_wbcomid DataFrame on 'COMID' and 'WBCOMID'
  ss_wbcomid <- merge(ss_wbcomid, df, by.x = "WBCOMID", by.y = "COMID", all.x = TRUE)
}

# merge remaining ss_data with lakecat data
lakecat_ss_data <- merge(ss_data, ss_wbcomid, by="SID", all=TRUE)


# save dataframe to file
file_path <- paste0("lakecat_ss_data_", version, ".csv")
write.csv(lakecat_ss_data, file=file_path, row.names=FALSE)  # this dataframe contains all the lakecat data merged with all the study site data

file_path <- paste0("ss_wbcomid_data_", version, ".csv")
write.csv(ss_wbcomid, file=file_path, row.names=FALSE)  # this dataframe contains just the study site SIDs, WBCOMIDs, and lakecat data

