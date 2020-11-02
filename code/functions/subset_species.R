subset_species <- function(species, catch, tows){
  OR_subset <- catch[catch$scientific_name == species, ]
  match_id <- match(tows$trawl_id, OR_subset$trawl_id)
  tows$lncpue <- OR_subset$lncpue_n[match_id]
  tows$lncpue[is.na(tows$lncpue)] <- 0
  selected_species <- select(tows, julian, year, lncpue, latitude, longitude)
  selected_species <- na.omit(selected_species)
  selected_species$pres <- 1 * (selected_species$lncpue > 0)
  return(selected_species)
}
