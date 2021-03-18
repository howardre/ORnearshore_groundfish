subset_species <- function(species, catch, tows){
  OR_subset <- catch[catch$scientific_name == species, ]
  match_id <- match(tows$trawl_id, OR_subset$trawl_id)
  tows$lncpue <- OR_subset$lncpue_n[match_id]
  tows$lncpue[is.na(tows$lncpue)] <- 0
  selected_species <- select(tows, julian, year, lncpue, latitude, longitude, PDO, NPGO,
                             bottom_temp, depth_m, program)
  selected_species <- na.omit(selected_species)
  selected_species$pres <- 1 * (selected_species$lncpue > 0)
  return(selected_species)
}

subset_species_temp <- function(species, catch, tows){
  OR_subset <- catch[catch$scientific_name == species, ]
  match_id <- match(tows$trawl_id, OR_subset$trawl_id)
  tows$lncpue <- OR_subset$lncpue_n[match_id]
  tows$lncpue[is.na(tows$lncpue)] <- 0
  tows$count <- as.integer(tows$lncpue)
  selected_species <- select(tows, julian, year, lncpue, latitude, longitude, PDO, NPGO,
                             bottom_temp, depth_m, count, program)
  selected_species <- na.omit(selected_species)
  selected_species$pres <- 1 * (selected_species$lncpue > 0)
  selected_species <- selected_species[selected_species$bottom_temp < 11,] # Filter out outlier temperature tows
  return(selected_species)
}
