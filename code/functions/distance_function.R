distance_function <- function(start_lat, start_lon, end_lat, end_lon)
{
  med_lat <- (start_lat + end_lat) / 2
  rad_lat <- (pi * med_lat) / 180
  shrink <- cos(rad_lat)
  delta_lat <- end_lat - start_lat
  delta_lon <- start_lon - end_lon
  mpermile <- 111195
  distance <- mpermile * sqrt((delta_lon * shrink) ^ 2 + (delta_lat) ^ 2)
  distance
}
