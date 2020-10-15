##############################################################################################
setwd('/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code/')

##############################################################################################
library(RODBC)

##############################################################################################
# Requires 32-bit R, cannot use 64-bit
# Connect to the 2007 Access database
con <- odbcConnectAccess2007("../data/ODFW_data/CiannelliDataRequest.accdb")
sqlTables(con)$BottomTrawlTickets # extract the fish tickets from the Access database
sqlTables(con)$BottomTrawlLogs # extract the logbooks from the Access database
sqlTables(con)$VesselData # extract the vessel data from the Access database

##############################################################################################
# Read the tables
fish_tickets <- sqlFetch(con, "BottomTrawlTickets")
str(fish_tickets) # display structure of the table
save(fish_tickets, file = "fish_tickets")
write.csv(fish_tickets,
  "../data/ODFW_data/fish_tickets.csv",
  row.names = FALSE)

logbooks <- sqlFetch(con, "BottomTrawlLogs")
str(logbooks) # display structure of the table
save(logbooks, file = "logbooks")
write.csv(logbooks,
  "../data/ODFW_data/logbooks.csv",
  row.names = FALSE)

vessel_data <- sqlFetch(con, "VesselData")
str(vessel_data) # display structure of the table
save(vessel_data, file = "vessel_data")
write.csv(vessel_data,
  "../data/ODFW_data/vessel_data.csv",
  row.names = FALSE)

odbcCloseAll() # close the database
