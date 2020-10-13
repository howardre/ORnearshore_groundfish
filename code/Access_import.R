setwd('/Users/howar/Documents/Oregon State/ORnearshore_groundfish')

# requires 32-bit R, cannot use 64-bit
library(RODBC)
con <- odbcConnectAccess2007("C:\\Users\\howar\\Documents\\Oregon State\\ORnearshore_groundfish\\data\\ODFW_data\\CiannelliDataRequest.accdb")
sqlTables(con)$BottomTrawlTickets

fish_tickets <- sqlFetch(con, "BottomTrawlTickets")
str(fish_tickets)
save(fish_tickets, file = "fish_tickets")
write.csv(fish_tickets,"C:\\Users\\howar\\Documents\\Oregon State\\ORnearshore_groundfish\\data\\ODFW_data\\fish_tickets.csv",row.names=FALSE)

logbooks <- sqlFetch(con, "BottomTrawlLogs")
str(logbooks)
save(logbooks, file = "logbooks")
write.csv(logbooks,"C:\\Users\\howar\\Documents\\Oregon State\\ORnearshore_groundfish\\data\\ODFW_data\\logbooks.csv",row.names=FALSE)

sqlTables(con)$VesselData
vessel_data <- sqlFetch(con, "VesselData")
str(vessel_data)
save(vessel_data, file = "vessel_data")
write.csv(vessel_data,"C:\\Users\\howar\\Documents\\Oregon State\\ORnearshore_groundfish\\data\\ODFW_data\\vessel_data.csv",row.names=FALSE)

odbcCloseAll()
