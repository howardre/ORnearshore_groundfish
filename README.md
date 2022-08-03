# Characterization of Oregon's Shelf Groundfishes and Associated Fisheries

This repo is for my M.S. thesis research and its associated manuscripts. There are two sub-projects/chapters.

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
- All other environmental data is found in the NMFS station data


### The effects of climate, oceanography, and habitat on the distribution and abundance of northern California Current continental shelf groundfishes
#### Non-metric multidimensional scaling (NMS)
NMS was used to investigate changes in community composition over time
- First make samples by species and environmental variables by species matrices [here](code/NMS_matrices.R/)
- Then run NMS for each survey [here](code/NMS_analysis.R/)

#### Generalized additive models (GAM)
##### Stationary formulation
Work in progress: there is no available code yet.

##### Nonstationary formulation/threshold GAM (TGAM)
Work in progress: only some code is available
- [TGAM function](code/TGAM_function.R/) WIP
- [TGAM bootstrap](code/TGAM_bootstrap.R/)

### Utility of combining fishery-independent and fishery-dependent data for spatiotemporal analyses in Oregon's groundfish fishery
This section is a work-in-progress and uses the _confidential_ ODFW data. If access to the data is available, it can be converted to a useable format [here](code/Access_import.R/) and then [QC](code/Logbook-Ticket_QC.R/)'d.
