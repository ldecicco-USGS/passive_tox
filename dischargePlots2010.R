library(dataRetrieval)

siteList <- read.csv("M:/QW Monitoring Team/GLRI toxics/Data Analysis/Passive samplers/SiteList_sampleDates.csv",
                     stringsAsFactors=FALSE)

padVariable <- function(x,padTo){
  numDigits <- nchar(x)
  if (padTo != numDigits){
    if ((padTo-numDigits)>0) {
      leadingZeros <- paste(rep("0",(padTo-numDigits)),collapse="",sep="")
      x <- paste(leadingZeros,x,sep="")
    }
  }
  return(x)
}

siteList$USGS.Station.ID <- as.character(sapply(siteList$USGS.Station.ID, function(x) padVariable(x=x,padTo=8)))
siteList$USGS.Station.ID["40851385" == siteList$USGS.Station.ID] <- "040851385"

sites <- siteList$USGS.Station.ID["000000NA" != siteList$USGS.Station.ID]

siteList$dischargeSites <- siteList$USGS.Station.ID
siteList$dischargeSites["04119400" == siteList$USGS.Station.ID] <- "04119000"
siteList$dischargeSites["04126010" == siteList$USGS.Station.ID] <- "04125550"
siteList$dischargeSites["04128500" == siteList$USGS.Station.ID] <- "04127997"
siteList$dischargeSites["04132052" == siteList$USGS.Station.ID] <- "04127997"
siteList$dischargeSites["04135020" == siteList$USGS.Station.ID] <- "04133501"
siteList$dischargeSites["04043000" == siteList$USGS.Station.ID] <- "04041500"
siteList$dischargeSites["04157000" == siteList$USGS.Station.ID] <- "04151500"
siteList$dischargeSites["04095090" == siteList$USGS.Station.ID] <- "04094000"
siteList$dischargeSites["04032000" == siteList$USGS.Station.ID] <- "04036000"

siteList$distanceToSite <- rep(0,nrow(siteList))
siteList$distanceToSite["04119400" == siteList$USGS.Station.ID] <- 18
siteList$distanceToSite["04126010" == siteList$USGS.Station.ID] <- 17.5
siteList$distanceToSite["04128500" == siteList$USGS.Station.ID] <- 9.5
siteList$distanceToSite["04132052" == siteList$USGS.Station.ID] <- 25.5
siteList$distanceToSite["04135020" == siteList$USGS.Station.ID] <- 10.3
siteList$distanceToSite["04043000" == siteList$USGS.Station.ID] <- 14.8
siteList$distanceToSite["04157000" == siteList$USGS.Station.ID] <- 12.3
siteList$distanceToSite["04095090" == siteList$USGS.Station.ID] <- 4.6
siteList$distanceToSite["04032000" == siteList$USGS.Station.ID] <- 11.5

# for (i in 1:length(sites)){
#   if(nchar(sites[i]) == 8){
#     cat(sites[i],"\n")
#     possibleError <- tryCatch(
#       whatUVAvail <- getDataAvailability(siteNumber=sites[i]),
#       error=function(e) e
#     )
#     
#     if(inherits(possibleError, "error")){
#       whatUVAvail <- rep(NA,11)
#     } else {
#       
#       whatUVAvail <- whatUVAvail["uv" == whatUVAvail$service & "00060" == whatUVAvail$parameter_cd,]
#       
#       if(nrow(whatUVAvail) == 0){
#         whatUVAvail <- rep(NA,11)
#       }
#       
#     }
#     if(1 == i){
#       dischargeAvailability <- data.frame(c(sites=sites[i],whatUVAvail),stringsAsFactors=FALSE)
#     } else {
#       dischargeAvailability <- rbind(dischargeAvailability, c(sites=sites[i],whatUVAvail))
#     }
#   }
# }

# Just checking:
# noData <- dischargeAvailability$sites[is.na(dischargeAvailability$parameter_cd)]
# for (i in 1:length(noData)){
#   if(nchar(noData[i]) == 8){
#     cat(noData[i],"\n")
#     possibleError <- tryCatch(
#       whatUVAvail <- getDataAvailability(siteNumber=noData[i]),
#       error=function(e) e
#     )
#     
#     if(inherits(possibleError, "error")){
#       whatUVAvail <- rep(NA,11)
#     } else {
#       if(1 == i){
#         noDischarge <- data.frame(cbind(sites=rep(noData[i],nrow(whatUVAvail)),whatUVAvail),stringsAsFactors=FALSE)
#       } else {
#         noDischarge <- rbind(noDischarge, data.frame(cbind(sites=rep(noData[i],nrow(whatUVAvail)),whatUVAvail,stringsAsFactors=FALSE)))
#       }
#     }
#   }
# }

sampleStartDate <- as.Date(strptime(siteList$Date.Deployed["000000NA" != siteList$USGS.Station.ID],format="%m/%d/%Y"))
sampleEndDate <- as.Date(strptime(siteList$Date.Retrieved["000000NA" != siteList$USGS.Station.ID],format="%m/%d/%Y"))

retrievalStartWindow <- 5
retrievalEndWindow <- 5


dischargeData <- list()

for (i in 1:length(sites)){
  cat(sites[i],"\n")
  StartDate <- sampleStartDate[i] - retrievalStartWindow
  EndDate <- sampleEndDate[i] + retrievalEndWindow
  if(is.na(EndDate)) EndDate <- "2010-11-01"
  if(siteList$dischargeSites[siteList$USGS.Station.ID == sites[i]] == siteList$USGS.Station.ID[siteList$USGS.Station.ID == sites[i]]) {
    site <- sites[i]
  } else {
    site <- siteList$dischargeSites[siteList$USGS.Station.ID == sites[i]]
  }
  
  possibleError <- tryCatch(
    discharge <- retrieveUnitNWISData(site,"00060",
                                      as.character(StartDate),
                                      as.character(EndDate)
                                      ),
    error=function(e) e
  )
  
  if(inherits(possibleError, "error")){
    discharge <- NA
    cat("Error in retrieving data at", site,"\n")
  }
  dischargeData[[sites[i]]] <- discharge
  
}



dischargeDaily <- retrieveNWISData("04029990","00060",
                                       as.character(sampleStartDate["04029990" == sites] - retrievalStartWindow),
                                       as.character(sampleEndDate["04029990" == sites] + retrievalEndWindow)
)
names(dischargeDaily)["datetime"  == names(dischargeDaily)] <- "dateTime"
dischargeData[["04029990"]] <- dischargeDaily

dischargeDaily <- retrieveNWISData("04108660","00060",
                                   as.character(sampleStartDate["04108660" == sites] - retrievalStartWindow),
                                   as.character(sampleEndDate["04108660" == sites] + retrievalEndWindow)
)
names(dischargeDaily)["datetime"  == names(dischargeDaily)] <- "dateTime"
dischargeData[["04108660"]] <- dischargeDaily


rm(possibleError,whatUVAvail, dischargeDaily)

save(dischargeData, file="dischargeData.RData")
save(siteList, file="siteList.RData")
#########################################################################
load("//igsarmewfsapa/projects/QW Monitoring Team/GLRI toxics/Data Analysis/Passive samplers/dischargeData.RData")
load("//igsarmewfsapa/projects/QW Monitoring Team/GLRI toxics/Data Analysis/Passive samplers/siteList.RData")

sampleStartDate <- as.Date(strptime(siteList$Date.Deployed["000000NA" != siteList$USGS.Station.ID],format="%m/%d/%Y"))
sampleEndDate <- as.Date(strptime(siteList$Date.Retrieved["000000NA" != siteList$USGS.Station.ID],format="%m/%d/%Y"))

retrievalStartWindow <- 5
retrievalEndWindow <- 5

sites <- siteList$USGS.Station.ID["000000NA" != siteList$USGS.Station.ID]

sortedNames <- sort(siteList$Site.ID)

pdf("passiveSamplerDischarges.pdf")
for (i in sortedNames){
  
  siteID <- siteList$USGS.Station.ID[i == siteList$Site.ID]
  
  discharge <- dischargeData[[siteID]]
  
  mainTitle <- paste(i,siteID,sep=": ")
  
  if ("04029990" == siteID || "04108660" == siteID) {
    names(discharge)["datetime"  == names(discharge)] <- "dateTime"
    discharge$dateTime <- as.POSIXct(discharge$dateTime)
  } else if (siteList$dischargeSites[siteList$USGS.Station.ID == siteID] != siteList$USGS.Station.ID[siteList$USGS.Station.ID == siteID]){
    mainTitle <- paste(siteList$Site.ID[siteID == siteList$USGS.Station.ID],": ",siteID,"\nDischarge at:", siteList$dischargeSites[siteList$USGS.Station.ID == siteID],
                       ", ", siteList$distanceToSite[siteList$USGS.Station.ID == siteID], " miles from the site",sep="")
  }
  
  if(is.na(discharge) || nrow(discharge) == 0 || "000000NA" == siteID){
    plot(1, type="n", axes=F, xlab="", ylab="",main=paste(siteList$Site.ID[siteID == siteList$USGS.Station.ID],siteID,sep=": "))
    box()
  } else {
    
    start <- as.POSIXct(sampleStartDate[siteID == sites])
    end <- as.POSIXct(sampleEndDate[siteID == sites])
    
    if(is.na(end)) end <- max(discharge$dateTime, na.rm=TRUE)
    
    plot(discharge$dateTime, discharge[,4],type="l",ylim=c(0,max(discharge[,4],na.rm=TRUE)),
         xlab="Date",ylab="Discharge[cfs]",main=mainTitle
    )
    
    rowsToShade <- which(discharge$dateTime >= start & discharge$dateTime <= end)
    polygon(c(discharge$dateTime[rowsToShade[1]],discharge$dateTime[rowsToShade],discharge$dateTime[rowsToShade[length(rowsToShade)]]),
            c(-1000,discharge[rowsToShade,4],-1000),col="lightgrey")
    abline(v=start,col="blue",lwd=2)
    abline(v=end,col="blue",lwd=2)
    text(start+5*24*60*60, grconvertY(.08, from = "npc", to = "user"),labels=strptime(start,format="%Y-%m-%d"))
    text(end-5*24*60*60, grconvertY(.08, from = "npc", to = "user"),labels=strptime(end,format="%Y-%m-%d"))
    box()
    
  }
  
}
dev.off()