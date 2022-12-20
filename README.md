# Characterization of Oregon's Shelf Groundfishes and Associated Fisheries

The work contained here is for my M.S. thesis research and its associated manuscripts. This project was funded by Oregon Sea Grant and focuses on the groundfishes targeted by the Oregon nearshore trawl fleet.

### Data
The **NMFS survey data** raw and modified files are available in [here](data/NMFS_data/)
- [species_data_raw](data/NMFS_data/species_data_raw.csv/) contains the original data on species caught in a sample
- [station_data_raw](data/NMFS_data/station_data_raw.csv/) contains the original data for each haul/sample
- This data was retrieved from the [NMFS FRAM Data Warehouse](https://www.webapps.nwfsc.noaa.gov/data/map) for 1977-2018

The **Oregon logbook, fish ticket, and vessel data** used are _confidential_ and are therefore unavailable for public use
- This data was obtained with permission from the Oregon Department of Fish and Wildlife

Climate indices used are also available
- [NPGO](data/Environmental_data/NPGO.csv/)
- [PDO](data/Environmental_data/PDO/)
- [Grain size](data/Environmental_data/grain_df.Rdata) and [lithology](data/Environmental_data/lith_df.Rdata) though larger files can be sourced from elsewhere (some included in sediment folder)
- All other environmental data is found in the NMFS station data


### The effects of climate, oceanography, and habitat on the distribution and abundance of northern California Current continental shelf groundfishes
#### [Manuscript](https://doi.org/10.1111/fog.12553)
#### Non-metric multidimensional scaling (NMS)
NMS was used to investigate changes in community composition over time
- First make samples by species and environmental variables by species matrices [here](code/NMS_matrices.R/)
- Then run NMS for each survey [here](code/NMS_analysis.R/)

#### Generalized additive models (GAM)
##### Stationary formulation
- Separate the data into species specific data frames, ensuring all zero tows retained
- Initial analyses used basic [Gaussian GAMs](code/GAM_analysis.R/)
- Final analyses used [zero-inflated Poisson GAMs](code/ZIP_analysis.R/)
- Additional functions needed to run these scripts are available in [here](code/functions/)

##### Nonstationary formulation/threshold GAM (TGAM)
- Only basic formulations were used to run threshold GAMs, meaning no environmental terms
- Analyses are available [here](code/TGAM_analysis.R/) and the bootstrapping is available [here](code/TGAM_bootstrap.R/)
- There are relevant [TGAM functions](code/functions/TGAM_function.R/) used in the above scripts

### Comparing fishery-independent and fishery-dependent data for analysis of the distributions of Oregon shelf groundfishes
#### [Manuscript](https://doi.org/10.1016/j.fishres.2022.106553)
This section uses the _confidential_ ODFW data. If access to the data is available, it can be converted to a useable format [here](code/Access_import.R/) and then [QC](code/Logbook-Ticket_QC.R/)'d.
#### Visualization of Logbook Data
- General [effort](code/effort_visualization.R/) of the fleet and for [specific species](code/Logbook_visualization.R/) can be seen
- Data is gridded to ensure confidentiality
- Changes in [vessel use](code/Vessel_regression.R/) can also be visualized

#### Local index of collocation (LIC)
- This was calculated to look at overlap between the NMFS survey and Oregon logbook data
- See [Carroll et al. (2019)](https://doi.org/10.1111/geb.12984) for more information on other metrics
- The method used in this manuscript is calculated [here](code/LIC.R/)
