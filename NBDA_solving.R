######## NBDA analyses ###########
# modelling the diffusion of solving behaviour # 

# 1) Install and load libraries -------------------------------------------

# load NBDa via devtools
library(devtools)
# install_github("whoppitt/NBDA")
library(NBDA)
library(lubridate)

# 2) Data preparation -----------------------------------------------------


# 2.1) Load network data and puzzle data --------------------------------------------------

setwd("Data/")

load("Guett.gmm.RDA")
load("Mill.gmm.RDA")

load("Mill.puzzle.data.RDA")
load("Guett.puzzle.data.RDA")

puzzle.data.spring21.combined <- rbind(Guett.puzzle.data, Mill.puzzle.data)

# 2.2. Load breeder data --------------------------------------------------

# fledge.data 
fledge.data <- read.delim("fledge.dates.txt", sep=" ")
head(fledge.data)
# fledge.data:
# Box: Name of nest box
# PIT: alphanumeric code of individual (unique)
# Date.Time: yyyy-mm-dd HH:MM:SS
# Faceplated: 'yes' if we used a faceplate to determine fledge date and time, 'no' if we just went back to check on the day whether they had fledged


breeders <- read.delim("Breeders 2021.txt")
head(breeders)
# breeders:
# box: name of nest box
# Species: GRETI for great tit, BLUTI for blue tit
# Ring: unique metal leg ring of individual
# Tag: PIT tag of individual
# Who: 'Female' if breeding female, 'Male' if breeding male', 'Chick' for offspring
# Clutch.size: empty (only NAs)
# Chick.weight: weigth of chick when 15 days old (in g)
# Fledged: empty column
# Fledge.order: empty column
# Notes: notes
# Location: 'Castle' for one field site, 'Guett' for the other field site
# Assigned.side: 'rigth' of 'left' for adults only (which side they were restricted to on the puzzle box)
# Treatment: either 'same' if both parents assigned to same side, or 'opposite' if parents were assigned to opposite sides
# Assigned.day: day (in number of days since the 1st of April) on which the parent was restricted to one side
# Failed: 'yes' if the brood died before fledging

PITs <- read.delim("Species_PITs.txt")
# PITS:
# relevant columns are:
# Species: GRETI for great tit
# Pit: alphanumeric PIt tag


# Make a list of great tit chicks' PIT tags
GT.chicks <- subset(breeders$Tag, breeders$Who=="Chick" & breeders$Species=="GRETI" & breeders$Tag!="")
GT.chicks <- GT.chicks[!is.na(GT.chicks)]

# load all species list
species <- read.delim("Species_PITs.txt")
# make a list of all great tit tags
GT.list <- unique(subset(species$Pit, species$Species=="GRETI"))

# 2.3. Make a list of IDs to include --------------------------------------

# we include all birds who have been seen at the network feeders at least 10 times
# regardless of whether they solved the puzzle box or not

# Mill

Mill.IDs <- colnames(Mill.gmm$gbi[,colSums(Mill.gmm$gbi)>=10])
length(Mill.IDs)
# 234 tagged birds were recorded at the network feeders at least 10 times
# reduce to great tits only
Mill.IDs <- subset(Mill.IDs, Mill.IDs %in% GT.list)
length(Mill.IDs)
# 170 great tits with at least 10 sightings

# Guett
Guett.IDs <- colnames(Guett.gmm$gbi[,colSums(Guett.gmm$gbi)>=10])
length(Guett.IDs)
# 79 tagged birds were recorded at the network feeders at least 10 times

# reduce to great tits only
Guett.IDs <- subset(Guett.IDs, Guett.IDs %in% GT.list)
length(Guett.IDs)
# 61 great tits with at least 10 sightings

# 2.4. Create vertical social network -------------------------------------

# this function creates a matrix among all great tits with entries of 1 between parents and their offspring, and 0s otherwise
create.vert.net <- function(IDs){
  mat <- matrix(0, nrow=length(IDs), ncol=length(IDs))
  rownames(mat) <- IDs
  colnames(mat) <- IDs
  
  # chick IDs
  IDs.sub <- IDs[IDs %in% breeders$Tag & IDs %in% GT.chicks]
  
  for(i in IDs.sub){
    # extract the info for the bird i
    i.box <- subset(breeders$Box, breeders$Tag==i)
    # extract the parents of that chick
    parents <- subset(breeders$Tag, breeders$Box==i.box & breeders$Who %in% c("Male", "Female"))
    parents <- subset(parents, parents %in% IDs)
    # set entries between parents and offspring to 1 in the direction of the chick learning from parents
    # but not the other way round
    mat[i, parents] <- 1
    
  }  
  
 return(mat)   
}

Mill.vert.matrix <- create.vert.net(IDs = Mill.IDs)
Guett.vert.matrix <- create.vert.net(IDs = Guett.IDs)


# 2.5. Create oblique social network --------------------------------------

# creates a matrix with entries of 1 between non-parent adults and offspring
create.obl.net <- function(IDs){
  mat <- matrix(0, nrow=length(IDs), ncol=length(IDs))
  rownames(mat) <- IDs
  colnames(mat) <- IDs
  
  # chick IDs
  IDs.sub <- IDs[IDs %in% breeders$Tag & IDs %in% GT.chicks]
  
  for(i in IDs.sub){
    # extract the info for the bird i
    i.box <- subset(breeders$Box, breeders$Tag==i)
    parents <- subset(breeders$Tag, breeders$Who %in% c("Male", "Female") & breeders$Box== i.box)
    
    # remove parents and chicks from the ID lists to get adults
    adults <- subset(IDs, !(IDs %in% c(parents, GT.chicks, IDs.sub)))
    
    # set entries between non-parent adults and offspring to 1 in the direction of the chick learning from adults
    # but not the other way round
    mat[i, adults] <- 1
    
  }  
  
  return(mat)   
}

Mill.obl.matrix <- create.obl.net(IDs = Mill.IDs)
Guett.obl.matrix <- create.obl.net(IDs = Guett.IDs)


# 2.6. Create siblings social network -----------------------------------

# creates a matrix with entries of 1 among siblings
create.sib.net <- function(IDs){
  mat <- matrix(0, nrow=length(IDs), ncol=length(IDs))
  rownames(mat) <- IDs
  colnames(mat) <- IDs
  
  # chick IDs
  IDs.sub <- IDs[IDs %in% breeders$Tag & IDs %in% GT.chicks]
  
  for(i in IDs.sub){
    # extract the info for the bird i
    i.box <- subset(breeders$Box, breeders$Tag==i)
    # extract the parents of that chick
    siblings <- subset(breeders$Tag, breeders$Box==i.box & breeders$Who %in% c("Chick"))
    siblings <- subset(siblings, siblings %in% IDs)
    # set entries between siblings to 1 
    mat[i, siblings] <- 1
    mat[siblings, i] <- 1
    
  }  
  diag(mat) <- 0
  return(mat)   
}

Mill.sib.matrix <- create.sib.net(IDs = Mill.IDs)
Guett.sib.matrix <- create.sib.net(IDs = Guett.IDs)

length(Mill.sib.matrix[Mill.sib.matrix==1])
length(Guett.sib.matrix[Guett.sib.matrix==1])

# 2.7. Create horizontal social network -----------------------------------

# creates a matrix with entries of 1 among peers (non-siblings only)
create.hor.net <- function(IDs, sib.net){
  mat <- matrix(0, nrow=length(IDs), ncol=length(IDs))
  rownames(mat) <- IDs
  colnames(mat) <- IDs
  
  # chick IDs
  IDs.sub <- IDs[IDs %in% breeders$Tag & IDs %in% GT.chicks]
  
  # set connections among all chicks to 1
  # they can learn in both directions
  mat[IDs.sub, IDs.sub] <- 1
  
  # finally, set the sibling connections to 0
  mat[which(sib.net>0)] <-0 
    
  return(mat)   
}

Mill.hor.matrix <- create.hor.net(IDs = Mill.IDs, sib.net=Mill.sib.matrix)
Guett.hor.matrix <- create.hor.net(IDs = Guett.IDs, sib.net=Guett.sib.matrix)




# 2.8. Load demonstrators -------------------------------------------------

load("demonstrators.RDA")
demos

# these are all knowledgeable birds from Mill and Guett up until the 14th of May 2021


# 2.9. Extract date of acquisition for fledgies ---------------------------
# considering the third solve as the date of acquisition (for those who have solved at least 10 times in total)


# This function extracts the data, but we provide the created data set below 

# extract.acquisition.date <- function(puzzle.data){
#   sites <- sort(unique(puzzle.data$Location))
#   # extract those with at least three solves
#   
#   puzzle.data.sub <- subset(puzzle.data, puzzle.data$Event %in% c("right", "left"))
#   solvers <- as.data.frame(table(puzzle.data.sub$PIT))
#   solvers <- as.vector(subset(solvers$Var1, solvers$Freq>=10))
#   # remove demonstrators
#   solvers <- subset(solvers, !(solvers %in% demos ) & solvers %in% GT.list)
#   
#   # no extract the date of the first solves at each of the puzzle sites (if applicable)
#   # we consider the bird a learner at the site that it first produced the solution
#   # in all other sites, it will be filtered from the diffusion (it can still transmit behaviour but it is not considered in the model)
#   data.all <- as.data.frame(matrix(ncol=length(sites)+1))
#   
#   colnames(data.all) <- c("PIT", sites)
#   
#   for(i in solvers){
#     sub.i <- subset(puzzle.data.sub, puzzle.data.sub$PIT==i ) # consider the first solve as acquisition date
#     # ensure it's sorted according to date time
#     sub.i <- sub.i[order(sub.i$Date.Time),]
#     date.acq.i <- sub.i[1, "Date.Time"]
#     site.acq.i <- sub.i[1, "Location"]
#     
#     data.all[which(solvers==i), "PIT"] <- i
#     data.all[which(solvers==i), "first.solved"] <- date.acq.i
#     data.all[which(solvers==i), "site.solved"] <- site.acq.i
#     
#     # extract for each other site at which point it started solving
#     other.sites <- setdiff(sites, site.acq.i)
#     
#     for(j in other.sites){ # only consider the solves after it has learned
#       sub.i.j <- subset(puzzle.data.sub, puzzle.data.sub$PIT==i & puzzle.data.sub$Location==j & puzzle.data.sub$Date.Time>date.acq.i)
#       if(length(sub.i.j$PIT)==0){
#         date.first.solve.j <- NA
#       } else {
#         date.first.solve.j <- sub.i.j[1, "Date.Time"] # pick the first solve at each site
#       }
#       
#       data.all[which(solvers==i), j] <- date.first.solve.j
#     }
#     
#   }
#  
#   
#   
#   
#   return(data.all)
# }
# 
# 
# Mill.acquisition.date <- extract.acquisition.date(puzzle.data = Mill.puzzle.data)
# 
# Guett.acquisition.date <- extract.acquisition.date(puzzle.data = Guett.puzzle.data)
# 
# # remove site Guett 1 as no one has ever learned
# Guett.acquisition.date <- Guett.acquisition.date[,colnames(Guett.acquisition.date)!="Guett_1"]
# 
# save(Mill.acquisition.date, file="Mill.acquisition.date.RDA")
# save(Guett.acquisition.date, file="Guett.acquisition.date.RDA")

load("Mill.acquisition.date.RDA")
load("Guett.acquisition.date.RDA")

# columns explained:
head(Mill.acquisition.date)
# PIT: PIT code of individual
# Mill_1: if there is an entry, the individual has performed its first solve at location Mill_1 - and the date and time of the first solve is shown as yymmddHHMMSS. If it has learned elsewhere, it has NA as entry
# Mill_2-Mill_4: equivalent to Mill_1
# first.solved: info from Mill_1-M ill_4 in one column (when the bird has produced its first solve)
# site_solved: at which location it has produced its first solve


# 3) Prepare social networks ----------------------------------------------
# We generate networks for each 48 hour period separately after the 14th of May (this goes across 14 weeks, fledglings start entering the population in week 4)
library(asnipe)

# we provide details on how we have created the social networks, but input files can be loaded directly below

# these mark the start and end times of the 48 hour periods (and the corresponding week number)
# start.end.Mill <- rbind(c(210524090000, 210526090000, 1),
#                         c(210531093000, 210602093000, 2),
#                         c(210607083000, 210609083000, 3),
#                         c(210614090000, 210616090000, 4),
#                         c(210621090000, 210623090000, 5),
#                         c(210628084500, 210630084500, 6),
#                         c(210705093000, 210707093000, 7))
# 
# 
# start.end.Guett <- rbind(c(210526103000, 210528103000, 1),
#                          c(210602110000, 210604110000, 2),
#                          c(210609101500, 210611101500, 3),
#                          c(210616120000, 210618120000, 4),
#                          c(210623101500, 210625101500, 5),
#                          c(210630101500, 210702101500, 6),
#                          c(210707110000, 210709110000, 7))
# 
# # we write a function that automatically creates the social networks 
# # we create 8 networks (always with two weeks of data together, with exception of week 1 and week 7)
# # net1: week1
# # net2: week1+2
# # net3: week2+3
# # net4: week3+4
# # net5: week4+5
# # net6: week5+6
# # net7: week6+7
# # net8: week7
# 
# create.networks <- function(gmm, IDs, start.end) {
#   list <- NULL
#   
#   # add a column defining the middle time of the group event
#   midtime <-
#     as.POSIXct(as.character(gmm$metadata$Start),
#                format = "%y%m%d%H%M%S",
#                origin = "1970-01-01") +
#     abs((
#       as.POSIXct(
#         as.character(gmm$metadata$Start),
#         format = "%y%m%d%H%M%S",
#         origin = "1970-01-01"
#       ) - as.POSIXct(
#         as.character(gmm$metadata$End),
#         format = "%y%m%d%H%M%S",
#         origin = "1970-01-01"
#       )
#     ) / 2)
#   
#   gmm$metadata$Time <- midtime
#   
#   # subset gbi to IDs only
#   gmm$gbi <- gmm$gbi[, colnames(gmm$gbi) %in% IDs]
#   
#   weeks.for.net <- rbind(c(1),
#   c(1,2),
#   c(2,3),
#   c(3,4),
#   c(4,5),
#   c(5,6),
#   c(6,7),
#   c(7))
#   
#   for (i in 1:8) {
#     which.weeks <- unique(weeks.for.net[i,])
#     start.end.sub <- subset(start.end, start.end[,3] %in% which.weeks)
#     start.time <- start.end.sub[1,1]
#     end.time <- start.end.sub[length(which.weeks),2]
#    
#      n.i <- get_network(
#       association_data = gmm$gbi,
#       data_format = "GBI",
#       association_index = "SRI",
#       identities = colnames(gmm$gbi),
#       which_identities = IDs,
#       times = as.numeric(format(gmm$metadata$Time, format =
#                                   "%y%m%d%H%M%S")),
#       start_time = start.time,
#       end_time = end.time
#     )
#     list[[i]] <- n.i
#     
#   }
# return(list)  
# }
# 
# # run for Mill site
# Mill.networks <-
#   create.networks(gmm = Mill.gmm,
#                   IDs = Mill.IDs,
#                   start.end = start.end.Mill)
# # run for Guett site
# Guett.networks <-
#   create.networks(gmm = Guett.gmm,
#                   IDs = Guett.IDs,
#                   start.end = start.end.Guett)
# 
# # we again provide both objects
# save(Mill.networks, file = "Mill.networks.RDA")
# save(Guett.networks, file = "Guett.networks.RDA")

load("Mill.networks.RDA")
load("Guett.networks.RDA")

# both objects contain 8 lists with a network each


# 4) NBDA: learning to solve ------------------------------------------------------------------

# In preparation, we first need to extract the scrounging data
# 4.1. Extract scrounging data --------------------------------------------

# we set the doors to close 2 seconds after birds solve, meaning that we assume that birds scrounged if arriving within 2 seconds after a solve
scrounging.window <- 2

# next we extract all scrounging events 
scrounge.fun <- function(data.puzzle, location){
  data.puzzle.comb <- NULL
  for(k in location){
    # consider a scrounging opportunity if the bird arrived in the same second as the solve or one second after
    # given the solve was done by a different individual
    solve.times.k <- subset(data.puzzle$Date.Time, data.puzzle$Location==k & data.puzzle$Event %in% c("left", "right"))
    # extract the times 2 seconds, 1 second after the solve and on the solve time
    around.solve.times.k <- c(as.POSIXct(solve.times.k, format="%y%m%d%H%M%S")+scrounging.window, as.POSIXct(solve.times.k, format="%y%m%d%H%M%S")+(scrounging.window-1))
    around.solve.times.k <- format(around.solve.times.k, format="%y%m%d%H%M%S")
    # add the solves times (if birds arrived in the same second that was solved in)
    around.solve.times.k <- c(around.solve.times.k, solve.times.k)
    around.solve.times.k <- sort(around.solve.times.k)
    
    # extract all events within two seconds of a solve
    scrounging.df <- subset(data.puzzle, data.puzzle$Date.Time %in% around.solve.times.k & data.puzzle$Location==k)
    
    # these still include other solves - reduce it to arrivals and displacements within 2s of solves
    scrounging.df <- subset(scrounging.df, scrounging.df$Event %in% c("displacement", "arrival"))
    scrounging.df$Event <- rep("scrounging", length(scrounging.df$Event))
    
    # add them to the data as scrounging events
    data.puzzle.comb <- rbind(data.puzzle.comb, scrounging.df)    
    
  }
  data.puzzle.comb <- data.puzzle.comb[order(as.numeric(rownames(data.puzzle.comb))),]
  return(data.puzzle.comb)
}

# run the function on both Mill and Guettingen data (spring 2021)
Mill.scrounge <- scrounge.fun(data.puzzle = Mill.puzzle.data, location = unique(Mill.puzzle.data$Location))
Guett.scrounge <- scrounge.fun(data.puzzle = Guett.puzzle.data[Guett.puzzle.data$Location!="Guett_1",], location = unique(Guett.puzzle.data$Location[Guett.puzzle.data$Location!="Guett_1"]))

scrounge.data.combined <- rbind(Mill.scrounge, Guett.scrounge)

# 4.3. Prepare NBDA data objects for solving ------------------------------------------

# we first need to extract how many days after fledgling birds performed their first solve on average
acquisition.dates.fledgies <- rbind(Mill.acquisition.date[,c(1,6)], Guett.acquisition.date[,c(1,4)])
acquisition.dates.fledgies <- subset(acquisition.dates.fledgies, acquisition.dates.fledgies$PIT %in% GT.chicks)
for(i in acquisition.dates.fledgies$PIT){
  fledge.date.i <- subset(fledge.data$Date.Time, fledge.data$PIT==i)
  time.till.solve.i <- abs(as.numeric(difftime(as.POSIXct(fledge.date.i), as.POSIXct(as.character(acquisition.dates.fledgies$first.solved[acquisition.dates.fledgies$PIT==i]),  format="%y%m%d%H%M%S"))))
  acquisition.dates.fledgies[which(acquisition.dates.fledgies$PIT==i),"time.since.fl"] <- time.till.solve.i
}

av.time.since.fl.till.first.solve <- mean(acquisition.dates.fledgies$time.since.fl)
av.time.since.fl.till.first.solve
# [1] 15.05643
# on average, fledglings started solving 15 days after fledging

#assign in which order the siblings learned

for(i in unique(fledge.data$Box)){
  siblings <- subset(fledge.data$PIT, fledge.data$Box==i)
  
  siblings <- siblings[siblings %in% acquisition.dates.fledgies$PIT]
  
  if(length(siblings)>0){
    acquisition.dates.fledgies[which(acquisition.dates.fledgies$PIT %in% siblings),"Box"] <- i
  }
}
  
acquisition.dates.fledgies <- acquisition.dates.fledgies[order(acquisition.dates.fledgies$first.solved),]

for(i in unique(acquisition.dates.fledgies$Box)){
  acq.sub.i <- subset(acquisition.dates.fledgies, acquisition.dates.fledgies$Box==i)
  order.i <- seq(1:length(acq.sub.i$PIT))
  acquisition.dates.fledgies[which(acquisition.dates.fledgies$PIT %in% acq.sub.i$PIT), "order"] <- order.i
}




# the following function creates NBDA data objects

prepare.NBDA.object.solving <-
  function(vert.net,
           obl.net,
           hor.net,
           sib.net,
           acquisition.dates,
           demos,
           site,
           box,
           network.data,
           network.data.raw) {
    # we first remove the ones that were not part of the network
  acquisition.dates <- subset(acquisition.dates, acquisition.dates$PIT %in% rownames(vert.net))
  
  
  # extract the actual learners (that performed their third solve at this box)
  learners <- subset(acquisition.dates$PIT, acquisition.dates$site.solved==box)
  
  # for the NBDA, all are learners and filtered individuals combined as 'learners' but then the ones that learned elsewhere are filtered out
  sub.box <- subset(acquisition.dates, acquisition.dates$PIT %in% learners)
  
  # we need to sort them according to their date of acquisition
  # for that we first plug all the dates of acquisition into the respective box column
  sub.box[,box] <- ifelse(is.na(sub.box[,box]), sub.box$first.solved, sub.box[,box])
  
  # now order
  sub.box <- sub.box[order(sub.box[,box]), ]
  
  # add those who learned before the 14th of May to the demos 
  demos.comb <- c(demos, sub.box$PIT[which(as.numeric(sub.box$first.solved)<210514000000)])
  
  # and remove them from the box data and learners vector
  sub.box <- subset(sub.box, !(sub.box$PIT %in% demos.comb))
  learners <- learners[!(learners%in% demos.comb)]
  
  # create the order of acquisition vector
  # a vector with the a number referring to the birds position in the social network (in the order that they learned the behaviour)
  OAc <- NULL
  for(i in sub.box$PIT){
    pos.i <- which(rownames(vert.net)==i)
    OAc[which(sub.box$PIT==i)] <- pos.i
  }
  
  # as time of acquisition, we use the time (in days) since the 14th of May (first fledgie fledged)
  TAc <- as.numeric(difftime(as.POSIXct(sub.box[,box], format="%y%m%d%H%M%S"), as.POSIXct("210514000000", format="%y%m%d%H%M%S"), unit="days"))
  
  # now we need to create filters - we filter adults
  
  # add adults
  filter.IDs <- learners[!(learners %in% GT.chicks)]
  
  
  # next, we prepare the networks
  # we use a dynamic frame, which means that networks are updated at every acquisition event
  # the vertical network has entries of 1 from the fledge date until the date of independence, after independence, the entries change to the association strength based on SRI
  # the the oblique and horizontal network contain the association strength from the day of independence onwards
  

  #If we have dynamic (time-varying) networks we need to create a four dimensional array of size
  #no. individuals x no.individuals x no.networks x number of time periods
  #and provide an assMatrixIndex vector that specifies in which time period each acquisition event occured
  
  # the association array contains 4x num.acquisition events matrices
  num.acquisition.events <- length(OAc)
  
 matrix.array.full <- array(data = NA, dim=c(nrow(vert.net), ncol(vert.net), 4, num.acquisition.events))
 

 # first the vertical network
 for(i in 1:num.acquisition.events){
   event.date.time <- TAc[i]
   
   # for each acquisition event, we first assign the network data to be used
   # we take the network data that flanks the event
   # for the first few acquisition events, there is not network data before, so we just take the closest network data after
   # same for the final few acquisition events, for which we just take the one right before
   
   # Mill
   # network 1: 10.375-12.375
   # network 2: 17.39583-19.39583
   # network 3: 24.35417-26.35417
   # network 4: 31.375-33.375
   # network 5: 38.375-40.375
   # network 6: 45.36458-47.36458
   # network 7: 52.39583-54.39583
   
   # Guett
   # network 1: 12.4375-14.4375
   # network 2: 19.45833-21.45833
   # network 3: 26.42708-28.42708
   # network 4: 33.5-35.5
   # network 5: 40.42708-42.42708
   # network 6: 47.42708-49.42708
   # network 7: 54.45833-56.45833
   
   
   if(site=="Mill"){
     if (event.date.time < 10.375) {
       net.to.use <- 1
     } else if (12.375 < event.date.time &
                event.date.time < 17.39583) {
       net.to.use <- 2
     } else if (19.39583 < event.date.time &
                event.date.time <  24.35417) {
       net.to.use <- 3
     } else if (26.35417 < event.date.time &
                event.date.time <  31.375) {
       net.to.use <- 4
     } else if (33.375 < event.date.time &
                event.date.time <  38.375) {
       net.to.use <- 5
     } else if (40.375 < event.date.time &
                event.date.time <  45.36458) {
       net.to.use <- 6
     } else if (47.36458 < event.date.time &
                event.date.time <  52.39583) {
       net.to.use <- 7
     } else if (54.39583 < event.date.time) {
       net.to.use <- 8
     }
   } else if(site== "Guett"){
     if (event.date.time < 12.4375) {
       net.to.use <- 1
     } else if (14.4375 < event.date.time &
                event.date.time < 19.45833) {
       net.to.use <- 2
     } else if (21.45833 < event.date.time &
                event.date.time <  26.42708) {
       net.to.use <- 3
     } else if (28.42708 < event.date.time &
                event.date.time <  33.5) {
       net.to.use <- 4
     } else if (35.5 < event.date.time &
                event.date.time <  40.42708) {
       net.to.use <- 5
     } else if (42.42708 < event.date.time &
                event.date.time < 47.42708) {
       net.to.use <- 6
     } else if (49.42708 < event.date.time &
                event.date.time <  54.45833) {
       net.to.use <- 7
     } else if (56.45833 < event.date.time) {
       net.to.use <- 8
     }
   }
   
   # this social netowrk is based on the SRI values for the given acquisition event
   SRI.network <- network.data[[net.to.use]]
   
   # now we prepare the vertical social network
   vert.net.temp <- vert.net*SRI.network
  
  # now we plug this into the matrix array at the right slot
  matrix.array.full[,,1,i] <- vert.net.temp
  

   obl.net.temp <- obl.net
   
   # we first create the non-sib network
   non.sib.net <- sib.net-1
   non.sib.net <- abs(non.sib.net)
   # now we simply need to set the rows and columns of non-fledglings to 0
   non.sib.net[!(rownames(non.sib.net)%in% GT.chicks),] <- 0
   non.sib.net[,!(colnames(non.sib.net)%in% GT.chicks)] <- 0
   

   hor.net.temp.nonsib <- SRI.network*non.sib.net
   hor.net.temp.sib <- SRI.network*sib.net
   obl.net.temp <- SRI.network*obl.net
     
   # now we plug this into the matrix array at the right slot
   matrix.array.full[,,2,i] <- obl.net.temp # set of matrices 2 are for oblique
   matrix.array.full[,,3,i] <- hor.net.temp.sib # set of matrices 3 are for horizontal from siblings
   matrix.array.full[,,4,i] <- hor.net.temp.nonsib # set of matrices 4 are for horizontal from non-siblings
   
 }
 # now we need to define an index vector to specify which time period each matrix belonged to
 # each acquisition event has a separate period, so we just create a vector with consecutive numbers
 
 assMatrixIndex<-c(1:num.acquisition.events)
 
  
 # we now generate the presence matrix - based on whether each individual was present in the past 7 days (at network feeders or puzzles)
 # if yes, it gets an entry of 1, if not, it gets an entry of 0
 PresMat<-matrix(1, nrow=length(rownames(vert.net.temp)), ncol=num.acquisition.events)
 
 # add a column to the network metadata with the box names
 network.data.raw.temp <- network.data.raw
 network.data.raw.temp$metadata$Box <- NA
 
 if(site=="Mill"){
   network.data.raw.temp$metadata$Box <- substr(network.data.raw.temp$metadata$Location, 1, 6)
 } else if(site=="Guett"){
   network.data.raw.temp$metadata$Box <- substr(network.data.raw.temp$metadata$Location, 1, 7)
 }
 
 network.data.raw.temp$metadata <- subset(network.data.raw.temp$metadata, network.data.raw.temp$metadata$Box==box)
 

   for(j in sort(sub.box[,box])){ # for each acquisition event
     # subset the raw group by individual matrices to only include the week before the acquisition event
     which.rows <- subset(network.data.raw.temp$metadata, network.data.raw.temp$metadata$Start <= as.numeric(j) & network.data.raw.temp$metadata$Start > (as.Date(as.POSIXct(j, format="%y%m%d%H%M%S"))-7))
     present.in.net.at.j <- names(which(colSums(network.data.raw.temp$gbi[as.numeric(rownames(which.rows)),])>0))
    # vector with the PIT codes of the birds that had been seen at least once at the network feeder in the 7 days before the acquistion event 
    
   # next, we extract whether the bird was registered on the puzzle box antenna in the 7 days before the acquisition event j
   if(site=="Mill"){
     puzzle.data.raw <- Mill.puzzle.data
   } else if(site=="Guett"){
     puzzle.data.raw <- Guett.puzzle.data
   }
   # subset to the 7 days before the acquisition event
   puzzle.sub.j <- subset(puzzle.data.raw, puzzle.data.raw$Date < as.POSIXct(j, format="%y%m%d%H%M%S") & puzzle.data.raw$Date > as.Date(j, format="%y%m%d%H%M%S")-7)
   present.at.puzzle.at.j <- unique(puzzle.sub.j$PIT)
   
   # now extract which rownames were neither present at the network feeders nor at the puzzle box
  # set all of those not present to 0 in the presence matrix 
   PresMat[which(!(rownames(vert.net.temp) %in% c(present.in.net.at.j, present.at.puzzle.at.j))),which(sort(sub.box[,box])==j)] <- 0

   }
 
 
 
 # now we extract transmission weights 
 # for each knowledgeable individual, we take the solving rate after they acquired the behaviour 

 count <- 1
 weights <- NULL
 for(o in rownames(vert.net)){
   if(o %in% demos.comb){
     acq.date.o <- 2105140000 
   } else {
     acq.date.o <- subset(acquisition.dates$first.solved, acquisition.dates$PIT==o)
   }
   
  solves.after.acq.o <-   subset(
       puzzle.data.spring21.combined$Event,
       puzzle.data.spring21.combined$PIT == o &
         puzzle.data.spring21.combined$Event %in% c("left", "right") &
         as.numeric(puzzle.data.spring21.combined$Date.Time) > acq.date.o
     )
   num.solves.o <- length(solves.after.acq.o)
   weights[count] <- num.solves.o
   count <- count+1
 }
 
 # now we divide the number of solves by the total number of experimental days (57)
 weights <- weights/57
 
  # extract demonstrators
  demos.box <- intersect(rownames(vert.net), demos.comb)

  # demos needs to be a 0/1 vector of the length of the included IDs
  demos.box.vector <- rep(0, length(rownames(vert.net)))
  demos.box.vector[which(rownames(vert.net)%in% demos.box)] <- 1
  
  # now we have to create the filter - it is a character vector consisting of the name of the oNBDA data object and the position of the individual in the network
  filter <- paste(box, which(rownames(vert.net)%in% filter.IDs), sep="_")
  
  # Finally, we prepare individual-level variables
  # we prepare two matrices with rows corresponding to the individuals and columns corresponding to the acquisition events
  # matrix for scrounging
  scrounge.mat <-  matrix(NA, ncol=length(TAc), nrow=length(rownames(SRI.network)))
  parents.solve.mat <-  matrix(NA, ncol=length(TAc), nrow=length(rownames(SRI.network)))
  sibs.knowledgeable.mat<-  matrix(NA, ncol=length(TAc), nrow=length(rownames(SRI.network)))
  
  # first, we extract the number of scrounging event an individual performed from the start of the experiment until the event occurred (we include this as a time-varying variable)
  for(d in 1:num.acquisition.events){
    # extract day of acquisition
    d.TAc <- TAc[d]
    d.date <- as.POSIXct("210514000000", format="%y%m%d%H%M%S")+seconds_to_period(d.TAc*24*60*60)
    
    # subset the scrounge data to before the acquisition event
   if(site=="Mill"){
     scrounge.data <- Mill.scrounge
   } else {
     scrounge.data <- Guett.scrounge
   }
    # subset the scrounge data
    d.scrounge.sub <- subset(scrounge.data, scrounge.data$Date<d.date)
    scrounge.vec <- NULL
    # we extract the number of scrounging events up until the date of acquisition
        for(q in rownames(SRI.network)){
      q.scrounge.events <- length(subset(d.scrounge.sub$Event, d.scrounge.sub$PIT==q))
      scrounge.vec[which(rownames(SRI.network)==q)] <- q.scrounge.events
        }
    # add it to the scrounge matrix
    scrounge.mat[,d] <- scrounge.vec
    
    # now we extract the number of parent solves between fledging and the acquisition date
    if(site=="Mill"){puzzle.d <- Mill.puzzle.data
    } else {
      puzzle.d <- Guett.puzzle.data
    }
    
    # subset the puzzle data to the acquisition date
    d.puzzle.sub <- subset(puzzle.d, puzzle.d$Date<d.date & puzzle.d$Event %in% c("right", "left"))
    parent.solve.vec <- NULL
    sibs.knowledgeable.vec <- NULL
    # now loop through all individuals
    for(q in rownames(SRI.network)){
      if(q %in% GT.chicks){
        # extract its parents
        box.q <- subset(breeders$Box, breeders$Tag==q)
        fledge.q <- subset(fledge.data$Date.Time, fledge.data$PIT==q)
        # extract the date of acquisition for the fledgie
        if(q %in% acquisition.dates$PIT){
          acq.date.q <- subset(acquisition.dates$first.solved, acquisition.dates$PIT==q)
        } else {
          acq.date.q <- as.Date(fledge.q) + 15
        }
        parents <- subset(breeders$Tag, breeders$Box==box.q & breeders$Who %in% c("Male", "Female"))
        q.puzzle.events <- length(subset(d.puzzle.sub$Event, d.puzzle.sub$PIT %in% parents & d.puzzle.sub$Date>fledge.q & d.puzzle.sub$Date<acq.date.q))
        
        # if it's a chick, then extract how many sibligns were knowledgeable at time d
        
        knowledgeable.sibs.d <- length(subset(acquisition.dates.fledgies$PIT, as.POSIXct(acquisition.dates.fledgies$first.solved, format="%y%m%d%H%M%S")<=d.date & acquisition.dates.fledgies$Box ==box.q))
        
        length(subset(acquisition.dates.fledgies$PIT, as.POSIXct(acquisition.dates.fledgies$first.solved, format="%y%m%d%H%M%S")<=d.date & acquisition.dates.fledgies$Box ==box.q))
        
        
        
      } else {
        q.puzzle.events <- 0
        knowledgeable.sibs.d <- 0
      }
      
      
      parent.solve.vec[which(rownames(SRI.network)==q)] <- q.puzzle.events
      sibs.knowledgeable.vec[which(rownames(SRI.network)==q)] <- knowledgeable.sibs.d
    }
    parents.solve.mat[,d] <- parent.solve.vec
    sibs.knowledgeable.mat[,d] <- sibs.knowledgeable.vec
  }
  
  #scrounge.mat.log <- log(scrounge.mat+1)
  
  scrounge.mat <- scrounge.mat-mean(scrounge.mat)
  scrounge.mat <- scrounge.mat/100
  
#  scrounge.mat.sd <- (scrounge.mat-mean(scrounge.mat))/100 # per 100 solves, centered around 0
#  scrounge.mat.sd <- (scrounge.mat-38.4281) # these correspond to overall mean and stdev
  
  #parents.solve.mat.log <- log(parents.solve.mat+1)
  
  parents.solve.mat <- parents.solve.mat-mean(parents.solve.mat)
  parents.solve.mat <- parents.solve.mat/100
  
#  parents.solve.mat.sd <- (parents.solve.mat-mean(parents.solve.mat))/1000
#  parents.solve.mat.sd <- (parents.solve.mat-292.4787)
  # convert to per 100 solves/scrounges
  
  #scrounge.mat.sd <- (scrounge.mat-0)/(1198-0)-0.5 # (x/min)/(max-min)-0.5
  # 1198 corresponds to the max scrounges across all sites
  #parents.solve.mat.sd <- (parents.solve.mat-0)/(4760-0)-0.5
  # 4760 is the max number of parental solves across all sites
  
  # take the log of both matrices
 # scrounge.mat.sd <- log(scrounge.mat)
#  scrounge.mat.sd[scrounge.mat.sd==-Inf] <- 0
  
#  parents.solve.mat.sd <- log(parents.solve.mat)
#  parents.solve.mat.sd[parents.solve.mat.sd==-Inf] <- 0
  

  # standardize both the scorunging events and the number of parental solves 
  #scrounge.mat.sd <- (scrounge.mat-min(scrounge.mat))/(max(scrounge.mat)-min(scrounge.mat))-0.5
  #parents.solve.mat.sd <- (parents.solve.mat-min(parents.solve.mat))/(max(parents.solve.mat)-min(parents.solve.mat))-0.5
  
  assign(paste("parents.solves", box, sep="_"),  parents.solve.mat, envir = .GlobalEnv)
  assign(paste("scrounges", box, sep="_"),  scrounge.mat, envir = .GlobalEnv)
 # assign(paste("kn.sibs", box, sep="_"),  sibs.knowledgeable.mat, envir = .GlobalEnv)
  ILVs <- c(paste("parents.solves",box, sep="_"), paste("scrounges",box, sep="_"))
  

    object <- nbdaData(label=box,
                       assMatrix = matrix.array.full,
                       asoc_ilv = c(paste("parents.solves", box, sep="_"), paste("scrounges", box, sep="_")),
                       #   multi_ilv = "ILVabsent",
                       int_ilv = c(paste("parents.solves", box, sep="_"), paste("scrounges", box, sep="_")),
                       orderAcq = OAc,
                       timeAcq = TAc,
                       demons = demos.box.vector,
                       assMatrixIndex = assMatrixIndex,
                       presenceMatrix = PresMat,
                       weights = weights,
                       asocialTreatment ="timevarying"
    )
  
  
 
  
  object.filtered <- filteredNBDAdata(nbdadata=object, filter="id", exclude=filter)
  
  # little summary
  summary.df <- rbind.data.frame(c("learners",length(learners[!(learners %in% filter.IDs)])),
                   c("filtered", length(filter)),
                   c("total.acq.steps", length(OAc)))
  colnames(summary.df) <- c("who", "number")
  object.to.return <- NULL
  object.to.return$nbda.data <- object.filtered
  object.to.return$summary <- summary.df
  object.to.return$parent.solves <- parents.solve.mat
  object.to.return$scrounges <- scrounge.mat
  object.to.return$kn.sibs <- sibs.knowledgeable.mat
  
  
  return(object.to.return)
  
}


# run to create all NBDA data objects
# (Note that Mill 1 is removed because there was only one learner)
# NOTE: running these are computationally intense - we therefore provide the NBDA data objects below to load directly
Mill2.nbda.object <- prepare.NBDA.object.solving(vert.net = Mill.vert.matrix,
                                         obl.net = Mill.obl.matrix,
                                         hor.net = Mill.hor.matrix,
                                         sib.net = Mill.sib.matrix,
                                         network.data = Mill.networks,
                                         network.data.raw = Mill.gmm,
                                         acquisition.dates = Mill.acquisition.date,
                                         demos = demos,
                                         site="Mill",
                                         box="Mill_2"
                                         )
#save(Mill2.nbda.object, file="Output/Mill2.nbda.object.RDA")
load("Output/Mill2.nbda.object.RDA")

Mill2.nbda.object$summary
# who number
# 1        learners     20
# 2        filtered      9
# 3 total.acq.steps     29

Mill3.nbda.object <- prepare.NBDA.object.solving(vert.net = Mill.vert.matrix,
                                         obl.net = Mill.obl.matrix,
                                         hor.net = Mill.hor.matrix,
                                         sib.net = Mill.sib.matrix,
                                         network.data = Mill.networks,
                                         network.data.raw = Mill.gmm,
                                         acquisition.dates = Mill.acquisition.date,
                                         demos = demos,
                                         site="Mill",
                                         box="Mill_3")

#save(Mill3.nbda.object, file="Output/Mill3.nbda.object.RDA")
load("Output/Mill3.nbda.object.RDA")

Mill3.nbda.object$summary
#               who number
# 1        learners     18
# 2        filtered      8
# 3 total.acq.steps     26


Mill4.nbda.object <- prepare.NBDA.object.solving(vert.net = Mill.vert.matrix,
                                         obl.net = Mill.obl.matrix,
                                         hor.net = Mill.hor.matrix,
                                         sib.net = Mill.sib.matrix,
                                         network.data = Mill.networks,
                                         network.data.raw = Mill.gmm,
                                         acquisition.dates = Mill.acquisition.date,
                                         demos = demos,
                                         site="Mill",
                                         box="Mill_4")
#save(Mill4.nbda.object, file="Output/Mill4.nbda.object.RDA")
load("Output/Mill4.nbda.object.RDA")

Mill4.nbda.object$summary
# who number
# 1        learners      2
# 2        filtered      3
# 3 total.acq.steps      5


Guett2.nbda.object <- prepare.NBDA.object.solving(vert.net = Guett.vert.matrix,
                                         obl.net = Guett.obl.matrix,
                                         hor.net = Guett.hor.matrix,
                                         sib.net = Guett.sib.matrix,
                                         network.data = Guett.networks,
                                         network.data.raw = Guett.gmm,
                                         acquisition.dates = Guett.acquisition.date,
                                         demos = demos,
                                         site="Guett",
                                         box="Guett_2")
#save(Guett2.nbda.object, file="Output/Guett2.nbda.object.RDA")
load("Output/Guett2.nbda.object.RDA")

Guett2.nbda.object$summary
# who number
# 1        learners      5
# 2        filtered      2
# 3 total.acq.steps      7

Guett3.nbda.object <- prepare.NBDA.object.solving(vert.net = Guett.vert.matrix,
                                         obl.net = Guett.obl.matrix,
                                         hor.net = Guett.hor.matrix,
                                         sib.net = Guett.sib.matrix,
                                         network.data = Guett.networks,
                                         network.data.raw = Guett.gmm,
                                         acquisition.dates = Guett.acquisition.date,
                                         demos = demos,
                                         site="Guett",
                                         box="Guett_3")
#save(Guett3.nbda.object, file="Guett3.nbda.object.RDA")
load("Guett3.nbda.object.RDA")

Guett3.nbda.object$summary
# who number
# 1        learners      8
# 2        filtered      3
# 3 total.acq.steps     11


# 4.4. Prepare constraints vector matrix ----------------------------------



# now we need to specify which models should be run, in the so-called constraints vector matrix
# we have a function that autoamtically generates it from the nbda data object
create.constraints.Vect.Matrix <- function(NBDA_data_object, num_networks, num_ILVs){
  suppressWarnings(
    if("ILVabsent" %in% NBDA_data_object@asoc_ilv){
      num.ILV.asoc <- 0
    } else {num.ILV.asoc <- length(NBDA_data_object@asoc_ilv)})
  
  suppressWarnings(
    if("ILVabsent" %in% NBDA_data_object@int_ilv){
      num.ILV.int<- 0
    } else {num.ILV.int<- length(NBDA_data_object@int_ilv)})
  
  suppressWarnings(
    if("ILVabsent" %in% NBDA_data_object@multi_ilv){
      num.ILV.multi <- 0
    } else {num.ILV.multi <- length(NBDA_data_object@multi_ilv)})
  
  vector <- seq(1:(num_networks+num.ILV.asoc+num.ILV.int+num.ILV.multi))
  
  count <- 0 # create an object 'count', which starts on 0
  
  constraintsVect <- matrix(nrow = 10000000, ncol=(num_networks+num.ILV.asoc+num.ILV.int+num.ILV.multi)) # create a matrix to save the combination of parameters in
  constraintsVect[1,] <- vector # the first row gets filled with a sequence from 1:8 (all parameters will be estimated, none are set to 0)
  
  for (i in 1:(length(vector)-1)){ # a loop for each number of parameters to be estimated
    array <- combn(vector, i, FUN = NULL, simplify = TRUE) # for each number of paramters to be estiamted (e.g. 2) create all possible combinations of numbers between 1:12 (e.g. 2&8, 1&5 etc)
    
    for (j in 1:length(array[1,])){ # for each of those combinations
      vector2 <- seq(1:((num_networks+(num.ILV.asoc+num.ILV.int+num.ILV.multi))-i)) # create a second vector with 11-i free spaces
      position <- array[,j] # for each created combination
      count <- count+1 # add +1 to the count
      
      for (k in position){ # at each possible position
        vector2 <- append(vector2, 0, after=k-1) # add a 0 (e.g. 1 0 2 3 ...; 1 2 0 3 4 5 ...; 1 2 3 0 4 5 ....)
      }
      constraintsVect[count+1,] <- vector2 # and save the resulting order in a matrix
    }
  }
  
  
  constraintsVect <- as.matrix(as.data.frame(na.omit(constraintsVect))) # remove all NAs from the matrix
  
  # extract which columns are networks
  col.networks <- c(1:num_networks)
  
  col.names <- NULL
  
  if(num.ILV.asoc!=0){
    col.names <- rep("asoc", num.ILV.asoc)
  }
  
  if(num.ILV.int!=0){
    col.names <- c(col.names, rep("int", num.ILV.int))
  }
  
  if(num.ILV.multi!=0){
    col.names <- c(col.names, rep("multi", num.ILV.multi))
  }
  
  colnames(constraintsVect) <- c(rep("network", num_networks), col.names)
  
  constraintsVect <- as.matrix(as.data.frame(constraintsVect))
  
  # extract the models containing any social network
  
  social.models <- rep(NA, length(constraintsVect[,1]))

  for (k in 1:length(constraintsVect[,1])){
    sum <- sum(constraintsVect[k,1:num_networks])
    if(sum!=0){
      social.models[k] <- k
    }
  }
  social.models <- as.vector(na.omit(social.models))

  social.models.matrix <- constraintsVect[social.models,]

  # if multiplicative models are fit, we need to adjust the matrix
  # if the multiplicative slots are filled, it automatically fits the parameter for asoc and social (just constrained to be the same)
  # meaning that we can remove it from the asoc and int slot

  if(num.ILV.multi!=0){
    social.models.retain <- rep(NA, length(social.model.matrix[,1]))
    multi.models <- rep(NA, length(social.models.matrix[,1]))
    for (k in 1:length(social.models.matrix[,1])){
      sum <- sum(social.models.matrix[k,which(colnames(social.models.matrix)=="multi")])
      sum2 <- sum(social.models.matrix[k, c(which(colnames(social.models.matrix)=="asoc"),which(colnames(social.models.matrix)=="int"))])
      if(sum!=0 & sum2==0){ # if multi models are fit and int and asoc are set to 0
        multi.models[k] <- k # then retain the model
      } else if (sum==0){
        social.models.retain[k] <- k
      }
    }

    multi.models <- as.vector(na.omit(multi.models))
    social.models.retain <- as.vector(na.omit(social.models.retain))

    models.to.retain <- c(multi.models, social.models.retain)

    # these models are retained
    retain.matrix.soc <- social.models.matrix[models.to.retain,]

    social.models.matrix <- retain.matrix.soc
  }

  # extract the models containing no social network

  asocial.models <- rep(NA, length(constraintsVect[,1]))

  for (k in 1:length(constraintsVect[,1])){
    sum <- sum(constraintsVect[k,1:num_networks])
    if(sum==0){
      asocial.models[k] <- k
    }
  }
  asocial.models <- as.vector(na.omit(asocial.models))

  asocial.models.matrix <- constraintsVect[asocial.models,]

  # cols.asoc <- which(colnames(constraintsVect)=="asoc")
  # 
  # asocial.retain <- rep(NA, length(asocial.models))
  # if(length(asocial.models)>0){
  #   for (k in 1:length(asocial.models)){
  #     sum <- sum(asocial.models.matrix[k,which(colnames(constraintsVect)!="asoc")])
  #     if(sum==0){
  #       asocial.retain[k] <- k
  #     }
  #   }  
  # }
  

  asocial.retain <- asocial.models.matrix
  asocial.retain <- as.vector(na.omit(asocial.retain))

#  asocial.models.to.retain <- asocial.models.matrix[asocial.retain, ]
#  asocial.models.to.retain.matrix <- as.matrix(asocial.models.to.retain, nrow=1, ncol=5)
  constraintsVectMatrix <- rbind(social.models.matrix,asocial.models.matrix)
  rownames(constraintsVectMatrix) <- NULL

  # add the Null model (without social learning, and no ILVs)
  constraintsVectMatrix <- rbind(constraintsVectMatrix, rep(0, length(constraintsVectMatrix[1,])))

  row.names(constraintsVectMatrix) <- NULL
  return(constraintsVectMatrix)
  
  
}


# we can run it on any NBDA data object
# it does not matter which NBDA data object we choose - the matrix is the same for all
# we simply specify the number of networks used and the number of ILVs 
constraintsVectMatrix <- create.constraints.Vect.Matrix(NBDA_data_object = Mill2.nbda.object$nbda.data, 
                                                        num_networks = 4, 
                                                        num_ILVs = 2)


colnames(constraintsVectMatrix) <- c("vert_net", "obl_net", "hor_net_sib", "hor_net_nonsib", "asoc_parent_solves", "asoc_scrounges", "soc_parent_solves", "soc_scrounges")

# or load directly
#save(constraintsVectMatrix, file="Output/constraintsVectMatrix.RDA")
load("Output/constraintsVectMatrix.RDA")

# we have a look at the output
head(constraintsVectMatrix)


# each row represents a model, each column a parameters
# the first four columns refer to the four networks
# columns 5-6 to the ILVs influencing asocial learning
# columns 7-8 to the ILVs influencing social learning
# if a parameter is set to 0, it is not estimated in that model
# if it is a number >0, then the parameter is estimated 
# we could in theory constrain model parameters to be the same 
# by setting equal numbers (e.g. vertical network =1, oblique network =1)
# but here, we want to estimate all parameters independently (hence, consecutive numbers for each additional parameter)


# we ran TADA with two baselines (constant and weibull)
# we here provide the code to do it:
# we first double up on the constraitsvector and create a vector with the baseline to use



# constraintsVectMatrix <- rbind(constraintsVectMatrix,constraintsVectMatrix)
# 
# dim(constraintsVectMatrix)
# 
# baseline <- c(rep("constant", 256), rep("weibull", 256))

# we then run the tada AIC function
# AIC.table <-
#   tadaAICtable(
#     nbdadata = list(
#       Mill2.nbda.object$nbda.data,
#       Mill3.nbda.object$nbda.data,
#       #     Mill4.nbda.object$nbda.data, # with only two learners, we remove Mill4
#       Guett2.nbda.object$nbda.data,
#       Guett3.nbda.object$nbda.data
#     ),
#     constraintsVectMatrix = constraintsVectMatrix,
#     baselineVect = baseline,
#     writeProgressFile = F,
#     cores = 6
#   )
# and look at the baseline support - we get overwhelming support for weibull baseline, so we rerun the main analysis with just weibull as the baseline
# baselineSupport(AIC.table)

# support numberOfModels
# constant 0.1457853            256
# weibull  0.8542147            256

# 4.5. Run TADA -----------------------------------------------------------
load("Output/constraintsVectMatrix.RDA")
baseline <- rep("weibull", 256)

AIC.table <-
  tadaAICtable(
    nbdadata = list(
      Mill2.nbda.object$nbda.data,
      Mill3.nbda.object$nbda.data,
 #     Mill4.nbda.object$nbda.data, # with only two learners, we remove Mill4
      Guett2.nbda.object$nbda.data,
      Guett3.nbda.object$nbda.data
),
    constraintsVectMatrix = constraintsVectMatrix, 
baselineVect = baseline,
    writeProgressFile = F,
 cores = 6
  )

#save(Output/AIC.table, file="tadaAICtable.RDA")

# or load the TADA object there directly
load("Output/tadaAICtable.RDA")

# 4.6. Compute summary results --------------------------------------------

# print the AIC tables (models are ordered according to performance by AICc)
print(AIC.table@printTable)
write.csv(AIC.table@printTable, file="Output/AIC.table.csv")

# get the network support
networkSupport <- networksSupport(AIC.table)
round(networkSupport,3)

# network 1: vertical learning
# network 2: oblique learning
# network 3: horizontal learning from siblings
# network 4: horizontal learning from non-siblings

# we can see that combination 1:2:3:0 gets most support (73 %), which suggests a combination of vertical, oblique and horizontal learning from siblings

# support numberOfModels
# 0:0:0:0   0.000             16
# 0:0:0:1   0.000             16
# 0:0:1:0   0.000             16
# 0:0:1:2   0.000             16
# 0:1:0:0   0.000             16
# 0:1:0:2   0.000             16
# 0:1:2:0   0.025             16
# 0:1:2:3   0.012             16
# 1:0:0:0   0.000             16
# 1:0:0:2   0.000             16
# 1:0:2:0   0.000             16
# 1:0:2:3   0.000             16
# 1:2:0:0   0.000             16
# 1:2:0:3   0.000             16
# 1:2:3:0   0.646             16
# 1:2:3:4   0.318             16

# 4.7 Create plot with network support ----------------------------------

old.combo <- rownames(networkSupport)

# for plotting, we want to rename the network combinations
# extract which networks are not included
temp <- gregexpr("0", old.combo) 
new.net.combo.names <- c()
for(i in 1:length(temp)){
  # remove the network that is 0 in that model
  nets <- c(1,3,5,7)
  nets <- nets[!nets %in% unlist(temp[i])]

  
  nets.new <- NULL
  
  
  
  
  if(1 %in% nets){
    nets.new <- c(nets.new, 1)
  } else {
    nets.new <- c(nets.new, 0)
  }
  if(3 %in% nets){
    nets.new <- c(nets.new, 2)
  }else {
    nets.new <- c(nets.new, 0)
  }
  if(5 %in% nets){
    nets.new <- c(nets.new, 3)
  }else {
    nets.new <- c(nets.new, 0)
  }
  if(7 %in% nets){
    nets.new <- c(nets.new, 4)
  }else {
    nets.new <- c(nets.new, 0)
  }


  
  new.net.combo.names[[i]] <- paste(unlist(nets.new), collapse = ":")
  
}

networkSupport <- cbind(networkSupport,as.character(as.data.frame(new.net.combo.names)))
colnames(networkSupport) <- c("support", "num.models", "net.comb.new")

# reorder by weight:
networkSupport <- networkSupport[ order(-networkSupport$support),]
networkSupport$net.comb.new <- factor(networkSupport$net.comb.new, levels = networkSupport$net.comb.new)

library(ggplot2)

# plot the network support:
 tiff("C:/Users/sonja/Desktop/Konstanz/Breeding Season 2021/Fledgie social learning/Figures/NBDA_support.tiff", width=800, height=500, units="px")



ggplot(aes(x=net.comb.new, y=support), data = networkSupport)+
  geom_bar(stat="identity", fill="grey95", col="black")+
  theme_classic()+
  ylim(c(0,0.8))+
  theme(axis.title.x = element_text(size = rel(1.3), angle = 00), 
        axis.title.y = element_text(size = rel(1.3)),
        axis.text.x = element_text(size = rel(1.3), angle = 45, vjust = 1, hjust=1))+
  xlab("network combination")+
  ylab("support (summed AICc weights)")+
#  ggtitle("A) Solving")+
  theme(title =element_text(size=14, face = "bold"), 
        axis.title=element_text(size=12,face="bold"))+
  theme(plot.title = element_text(hjust=-0.13))+

  annotate(geom="text", x=11.7, y=0.78, label="networks",
           color="black", fontface =2, cex=7)+
  annotate(geom="text", x=11.58, y=0.73, label="1: vertical",
           color="black", cex=6)+
  annotate(geom="text", x=11.6, y=0.68, label="2: oblique",
           color="black", cex=6)+
  annotate(geom="text", x=12.55, y=0.63, label="3: horizontal (siblings)",
           color="black", cex=6)+
  annotate(geom="text", x=12.9, y=0.58, label="4: horizontal (non-siblings)",
           color="black", cex=6)+
annotate(geom="text", x=11.9, y=0.53, label="0: no network",
         color="black", cex=6)

dev.off()


# 4.8 variable support and MLEs -------------------------------------------------------------------

# extract support for each variable
variable_support <- variableSupport(AIC.table, includeAsocial = T)
round(variable_support,3)

#           s1 s2 s3   s4 ASOC:parents.solves_Mill_2 ASOC:scrounges_Mill_2 SOCIAL:parents.solves_Mill_2 SOCIAL:scrounges_Mill_2
# support 0.964  1  1 0.33                      0.184                 0.219                         0.24                   0.411


## we can extract model averaged estimates for each variable
# the mean is often skewed, so the median may be more informative
MLE_med  <- modelAverageEstimates(AIC.table,averageType = "median")
round(MLE_med,2)


# s1                           s2                           s3                           s4 ASOCIALparents.solves_Mill_2 
# 51.74                         9.18                      5399.58                         0.00                         0.00 
# ASOCIALscrounges_Mill_2  SOCIALparents.solves_Mill_2       SOCIALscrounges_Mill_2 
# 0.00                         0.00                         0.00

# 4.9. Extract effect sizes  -------------------------------------------------------------------

# acording to the AIC table, model 214 is best performing including vertical, oblique and horizontal learning from siblings, and no ILVs
constraintsVectMatrix[214,]

best.model.no <- 214

# we constrain the nbda data to only include the parameters from the best model
Mill2.constrained.NBDA <-
  constrainedNBDAdata(Mill2.nbda.object$nbda.data, constraintsVect = constraintsVectMatrix[best.model.no, ])

Mill3.constrained.NBDA <-
  constrainedNBDAdata(Mill3.nbda.object$nbda.data, constraintsVect = constraintsVectMatrix[best.model.no, ])

# Mill4.constrained.NBDA <-
# constrainedNBDAdata(Mill4.nbda.object$nbda.data, constraintsVect = constraintsVectMatrix[best.model.no, ])

Guett2.constrained.NBDA <-
  constrainedNBDAdata(Guett2.nbda.object$nbda.data, constraintsVect = constraintsVectMatrix[best.model.no, ])

Guett3.constrained.NBDA <-
  constrainedNBDAdata(Guett3.nbda.object$nbda.data, constraintsVect = constraintsVectMatrix[best.model.no, ])


# run TADA on the best model
best.model <- tadaFit(nbdadata = list(
  Mill2.constrained.NBDA,
  Mill3.constrained.NBDA, 
#  Mill4.constrained.NBDA,
  Guett2.constrained.NBDA,
  Guett3.constrained.NBDA
), type="social",
baseline = "weibull")

# extract the percentage of birds that learned through social learning - by event
prop.solve.social.byevent <- oadaPropSolveByST.byevent( nbdadata = list(
  Mill2.constrained.NBDA,
  Mill3.constrained.NBDA,
#  Mill4.constrained.NBDA ,
  Guett2.constrained.NBDA,
  Guett3.constrained.NBDA 
), model = best.model)
prop.solve.social.byevent
# for each acquisition event, we now have an estimate for the likelihood of social learning  (P network)

# for each acquisition event, we add the information whether it was the first sibling of a cohort, or other siblings
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
# Mill2
for(i in 1:20){
  order.events <- c(8,9,10,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,
                   9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                   2, 3, 4, 6, 7,
                   3, 5, 6, 7, 8, 9, 10, 11)
  ID.order <- Mill2.nbda.object$nbda.data@orderAcq[order.events[i]]
  learner.PIT <- Mill.IDs[ID.order]
  # look up whether it was the first of its sibling cohort to learn or not
  order.i <- subset(acquisition.dates.fledgies$order, acquisition.dates.fledgies$PIT==learner.PIT)
  prop.solve.social.byevent[i, "order"] <- order.i
  prop.solve.social.byevent[i, "PIT"] <- learner.PIT
}

# Mill3
for(i in 21:38){
  order.events <- c(8,9,10,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,
                    9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                    2, 3, 4, 6, 7,
                    3, 5, 6, 7, 8, 9, 10, 11)
  ID.order <- Mill3.nbda.object$nbda.data@orderAcq[order.events[i]]
  learner.PIT <- Mill.IDs[ID.order]
  # look up whether it was the first of its sibling cohort to learn or not
  order.i <- subset(acquisition.dates.fledgies$order, acquisition.dates.fledgies$PIT==learner.PIT)
  if(length(order.i)==0){
    order.i <- NA
  }
  prop.solve.social.byevent[i, "order"] <- order.i
  prop.solve.social.byevent[i, "PIT"] <- learner.PIT
}

# Guett2
for(i in 39:43){
  order.events <- c(8,9,10,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,
                    9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                    2, 3, 4, 6, 7,
                    3, 5, 6, 7, 8, 9, 10, 11)
  ID.order <- Guett2.nbda.object$nbda.data@orderAcq[order.events[i]]
  learner.PIT <- Guett.IDs[ID.order]
  # look up whether it was the first of its sibling cohort to learn or not
  order.i <- subset(acquisition.dates.fledgies$order, acquisition.dates.fledgies$PIT==learner.PIT)
  if(length(order.i)==0){
    order.i <- NA
  }
  prop.solve.social.byevent[i, "order"] <- order.i
  prop.solve.social.byevent[i, "PIT"] <- learner.PIT
}
  
# Guett3
for(i in 44:51){
  order.events <- c(8,9,10,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,
                    9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                    2, 3, 4, 6, 7,
                    3, 5, 6, 7, 8, 9, 10, 11)
  ID.order <- Guett3.nbda.object$nbda.data@orderAcq[order.events[i]]
  learner.PIT <- Guett.IDs[ID.order]
  # look up whether it was the first of its sibling cohort to learn or not
  order.i <- subset(acquisition.dates.fledgies$order, acquisition.dates.fledgies$PIT==learner.PIT)
  if(length(order.i)==0){
    order.i <- NA
  }
  prop.solve.social.byevent[i, "order"] <- order.i
  prop.solve.social.byevent[i, "PIT"] <- learner.PIT
}

# 
prop.solve.social.byevent$order.cat <- as.numeric(prop.solve.social.byevent$order)
prop.solve.social.byevent$order.cat[prop.solve.social.byevent$order.cat!=1] <- 'other'
prop.solve.social.byevent$order.cat[prop.solve.social.byevent$order.cat==1] <- 'first'

# for plotting
plot.prop <- rbind.data.frame(cbind(prop.solve.social.byevent$`P(Network 1)`,"vertical",prop.solve.social.byevent$order.cat),
                 cbind(prop.solve.social.byevent$`P(Network 2)`,"oblique",prop.solve.social.byevent$order.cat),
                 cbind(prop.solve.social.byevent$`P(Network 3)`,"horizontal_siblings",prop.solve.social.byevent$order.cat)
                  )

colnames(plot.prop) <- c("support", "network", "sibling_category")
plot.prop$support <- as.numeric(plot.prop$support)
# 

tiff("C:/Users/sonja/Desktop/Konstanz/Breeding Season 2021/Fledgie social learning/Figures/sibling_cohorts.tiff", width = 1800, height = 1520, units="px", res = 300)

ggplot(aes(x=network, y=support, col=sibling_category, fill=sibling_category), data=plot.prop)+
  geom_boxplot()+
  theme_bw()+
#  scale_color_grey()+
  scale_fill_manual(values = c("white", "grey60"))+
  scale_color_manual(values = c("grey20", "grey20"))+
  ylab("Probability of social learning through network")+
  xlab("Network")+
 # ggtitle("B)")+
  theme(text = element_text(size=12))+
  theme(plot.title = element_text(hjust = 0, size=14))+
  theme(plot.title.position = "plot")+
  scale_x_discrete(limits=c("vertical", "oblique", "horizontal_siblings"),
                   labels=c("vertical", "oblique", "horizontal (siblings)"))+
  labs(col="sibling category")+
  guides(col="none", fill=guide_legend(title="sibling category"))

dev.off()

# get 95% CIs for plot.prop
library(dplyr)

# Assuming your data frame is called df
result <- plot.prop %>%
  group_by(network, sibling_category) %>%
  summarise(
    mean = mean(support),
    se = sd(support) / sqrt(n()),
    ci_lower = mean - qt(0.975, df = n() - 1) * se,
    ci_upper = mean + qt(0.975, df = n() - 1) * se
  )

print(as.data.frame(result))

# network sibling_category       mean         se   ci_lower  ci_upper
# 1 horizontal_siblings            first 0.00000000 0.00000000 0.00000000 0.0000000
# 2 horizontal_siblings            other 0.76704941 0.05823825 0.64755439 0.8865444
# 3             oblique            first 0.66626426 0.07001320 0.52106577 0.8114628
# 4             oblique            other 0.14129353 0.04881827 0.04112672 0.2414603
# 5            vertical            first 0.33016799 0.07012209 0.18474368 0.4755923
# 6            vertical            other 0.08977781 0.02661045 0.03517768 0.1443779

# then averaged across all acquisition events
prop.solve.social <- oadaPropSolveByST( nbdadata = list(
  Mill2.constrained.NBDA,
  Mill3.constrained.NBDA,
#  Mill4.constrained.NBDA ,
  Guett2.constrained.NBDA,
  Guett3.constrained.NBDA 
), model = best.model)
prop.solve.social

# P(Network 1) P(Network 2) P(Network 3)  P(S offset) 
# 0.19819      0.37805      0.42113      0.00000

1-sum(prop.solve.social)

# distinguish into first siblings and other siblings
# we have 23 'first' siblings
median(subset(prop.solve.social.byevent$`P(Network 1)`, prop.solve.social.byevent$order.cat=="first"))
# 22.8%
median(subset(prop.solve.social.byevent$`P(Network 2)`, prop.solve.social.byevent$order.cat=="first"))
# 76.6% oblique from no-parent adults
median(subset(prop.solve.social.byevent$`P(Network 3)`, prop.solve.social.byevent$order.cat=="first"))
# 0 percent horizontal from siblings


# subsequent siblings
median(subset(prop.solve.social.byevent$`P(Network 1)`, prop.solve.social.byevent$order.cat=="other"))
# 1.5%
median(subset(prop.solve.social.byevent$`P(Network 2)`, prop.solve.social.byevent$order.cat=="other"))
# 3.7% oblique from no-parent adults
median(subset(prop.solve.social.byevent$`P(Network 3)`, prop.solve.social.byevent$order.cat=="other"))
# 93.9 percent horizontal from siblings

# how many asocial events?
1-sum(prop.solve.social)
#[1] 0.00263

cbind.data.frame(best.model@varNames, best.model@outputPar, best.model@se)
# best.model@varNames best.model@outputPar best.model@se
# 1         Scale (1/rate):         5.686044e-04           NaN
# 2                   Shape         1.145528e-03  1.402370e-03
# 3 1 Social transmission 1         5.174167e+01  6.965600e+01
# 4 2 Social transmission 2         9.183973e+00  1.153759e+01
# 5 3 Social transmission 3         5.399585e+03  6.719099e+03

# we look at the baseline rate:

curve(dweibull(x, best.model@outputPar[2], best.model@outputPar[1]), from=1, to=30)
# dweibull takes arguments shape and then scale (second and first parameter of the output)
# the baseline indicates that the likelihood of learning was highest in the first days after fledging and then decreases and flattens out over time

# now we look at s parameters - extracting confidence intervals
# for plotting the profile likelihoods, we use the argument 'which'- this refers the parameter number in the best model and is detailed in the output just above (e.g. 1 = social transmission 1, which would result in which=1 for the code below)


# the output parameters from above provide some good starting points to set the intervals for the confidence intervals
# due to convergence issues (zigzag lines), we only extract the lower bounds
best.model@optimisation

# for vertical learning (which=1)
plotProfLik(which=1, model=best.model, range=c(0, 400), resolution = 20)
# we first zoom in on the lower limit a bit more
plotProfLik(which=1, model=best.model, range=c(0,5), resolution = 10)
# then the upper limit 
plotProfLik(which=1, model=best.model, range=c(20,  200), resolution = 10)
# we find that there are convergence issues for the upper limit 
# we therefore can only obtain the lower CI for the vertical network
CI.s1 <- profLikCI(which=1, model=best.model, lowerRange = c(1,3))
# Lower CI Upper CI 
# 1.991223       NA 


# for oblique learning (which=2)
plotProfLik(which=2, model=best.model, range=c(0,3500), resolution = 20)
plotProfLik(which=2, model=best.model, range=c(0,5), resolution = 10)
plotProfLik(which=2, model=best.model, range=c(3210,3220), resolution = 10)
CI.s2 <- profLikCI(which=2, model=best.model, lowerRange = c(0,1), upperRange = c(3212,3215))
# Lower CI     Upper CI 
# 0.3870949 3212.0323769

# for horizontal learning (from siblings) (which=3)
plotProfLik(which=3, model=best.model, range=c(0, 20000), resolution = 10)
# we zoom in on the lower limit
plotProfLik(which=3, model=best.model, range=c(270, 276), resolution = 10)
# even though it looks like the model is unable to draw a line, it seems like it can estimate the upper limit too. we zoom in
plotProfLik(which=3, model=best.model, range=c(11750, 11820), resolution = 10)
CI.s3 <- profLikCI(which=3, model=best.model, lowerRange = c(270,274), upperRange = c(11760, 11780))
# Lower CI   Upper CI 
# 272.3554 11770.7164 

# Now we want to extract the boundS for the different s parameters in percentages, which are easier to interpret
# note that this only works for the best performing model that includes all four networks
extract.percentages <- function(bound, s){
  
  # first we determine which bound (upper or lower)
  # lower uses the first element in the CI object, upper the second
  if(bound=="lower"){
    which.ci <- 1
  } else if(bound=="upper"){
    which.ci <- 2
  }
  
  # we need to specify the offset vectors
  # for extracting s, we add the CI as an offset, all else is set to 0
  if(s==1){
    offset <- c(CI.s1[which.ci], rep(0,7))

  } else if(s==2){
    offset <- c(0, CI.s2[which.ci], rep(0,6))

  } else if(s==3) {
    offset <- c(0, 0, CI.s3[which.ci], rep(0,5))

  } 
  
  # we need to define the constraints vector
  # we set the network in question to 0, but as it has an offset specified, it sets the value to that offset
  # the other networks and the ILV are estimated in the model
  if(s==1){
    constraintsV <- c(0,1,2,0,0,0,0,0)
  } else if(s==2){
    constraintsV <- c(1,0,2,0,0,0,0,0)
  } else if(s==3){
    constraintsV <- c(1,2,0,0,0,0,0,0)
  } 
  
  # we define the constrained NBDA data objects with the defined constraint and the offset
  Mill2.best.model.data.bound <- constrainedNBDAdata(
  nbdadata =
   Mill2.nbda.object$nbda.data,
  constraintsVect = constraintsV,
  offset = offset
)

Mill3.best.model.data.bound <- constrainedNBDAdata(
  nbdadata =
    Mill3.nbda.object$nbda.data,
  constraintsVect = constraintsV,
  offset = offset
)

# Mill4.best.model.data.bound <- constrainedNBDAdata(
#   nbdadata =
#     Mill4.nbda.object$nbda.data,
#   constraintsVect = constraintsVectMatrix[best.model.no, ],
#   offset = offset
# )

Guett2.best.model.data.bound <- constrainedNBDAdata(
  nbdadata =
    Guett2.nbda.object$nbda.data,
  constraintsVect = constraintsV,
  offset = offset
)

Guett3.best.model.data.bound <- constrainedNBDAdata(
  nbdadata =
    Guett3.nbda.object$nbda.data,
  constraintsVect = constraintsV,
  offset = offset
)

# now we fit an social model -as we have set the network in question to 0 in the constraintsVect, it takes the value of the offset. the other 2 networks are estimated

bestModelbound <-
  tadaFit(
    list(
      Mill2.best.model.data.bound,
      Mill3.best.model.data.bound,
#      Mill4.best.model.data.bound ,
      Guett2.best.model.data.bound,
      Guett3.best.model.data.bound 
    ) ,
    type = "social",
baseline = "weibull"
  )


prop.solve.social.bound <-
  oadaPropSolveByST(
    model = bestModelbound ,
    nbdadata = list(
      Mill2.best.model.data.bound,
      Mill3.best.model.data.bound,
#      Mill4.best.model.data.bound ,
      Guett2.best.model.data.bound,
      Guett3.best.model.data.bound 
    )
  )

# the value we are after is now in the offset (position 4)

return(paste("For network ", s, ", the ", bound, " bound is at ", 100*round(prop.solve.social.bound[3],3), "%.", sep="" ))

}

# s1 lower
extract.percentages(bound="lower", s=1)
#[1] "For network 1, the lower bound is at 14.0%."
#extract.percentages(bound="upper", s=1)
# NA


# s2 lower
extract.percentages(bound="lower", s=2)
# [1] "For network 2, the lower bound is at 32.2%."
extract.percentages(bound="upper", s=2)
# [1] "For network 2, the upper bound is at 38.3%."

# s3 lower
extract.percentages(bound="lower", s=3)
# [1] "For network 3, the lower bound is at 41.6%."
extract.percentages(bound="upper", s=3)
#  "For network 3, the upper bound is at 43.9%."



# 5) Plot social networks -------------------------------------------------


  # for illustrative purposes, we are plotting the networks from MI

  # subset gbi to IDs only
Mill.gmm.plot <- Mill.gmm  
Mill.gmm.plot$gbi <- Mill.gmm.plot$gbi[, colnames(Mill.gmm.plot$gbi) %in% Mill.IDs]

     net.plot <- get_network(
      association_data = Mill.gmm.plot$gbi,
      data_format = "GBI",
      association_index = "SRI",
      identities = colnames(Mill.gmm.plot$gbi),
      which_identities = Mill.IDs)
     
# extract ages
IDs <- colnames(net.plot)
chicks <- subset(breeders$Tag, breeders$Who=="Chick")
age.vec.plot <- NULL

for(i in IDs){
  if(i %in% chicks){
    age.vec.plot[which(IDs==i)] <- "J"
  } else {
    age.vec.plot[which(IDs==i)] <- "A"
  }
}

# extract knowledge status at the end of the experiment
knowledgeable <- unique(c(demos, acquisition.dates.fledgies$PIT, Mill.acquisition.date$PIT))

knowledge.plot <- NULL

for(i in IDs){
  if(i %in% knowledgeable){
    knowledge.plot[which(IDs==i)] <- "yes"
  } else {
    knowledge.plot[which(IDs==i)] <- "no"
  }
}

# for better visibility, we set the edge weights under 0.01 to 0
net.plot2 <- net.plot
net.plot2[net.plot2<0.01] <- 0
     
shape <- age.vec.plot
color <- knowledge.plot
library(igraph)
library(qgraph)

plot.network <- function(which.network){
  g.net <- igraph::graph_from_adjacency_matrix(net.plot2, mode = "undirected",
                                               weighted = TRUE, diag = TRUE)
  

  E(g.net)$width <- E(g.net)$weight
  
  color[color=="no"] <- "#00177e" # not knowledgeable
  color[color=="yes"] <- "#b23714" # knowledgeable
  
  shape[shape=="J"] <- "circle" # grey
  shape[shape=="A"] <- "square" # white
  
  V(g.net)$colour <- color
  V(g.net)$shape <- shape
  
  e <- as_edgelist(g.net,names=FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g.net),
                                         area=8*(vcount(g.net)^2),repulse.rad=(vcount(g.net)^2.9))
  
  # now we extract which edges are 0 in the which.matrix (e.g. all non-vertical connections), and set those as transparent.
  
  # Step 1: Get edge list (igraph uses numeric IDs)
  edge_list <- as_edgelist(g.net)
  
  # Step 2: Extract corresponding alpha values from your matrix
  # You need to map the edge pairs to the matrix values
  edge_alphas <- apply(edge_list, 1, function(e) {
    which.network[e[1], e[2]]
  })
  
  # Step 3: Set edge colors based on alpha
  # 0 → transparent, >0 → black opaque
  edge_colors <- ifelse(edge_alphas == 0,
                        rgb(0, 0, 0, alpha = 0),   # transparent
                        rgb(0.2, 0.2, 0.2, alpha = 1))  # opaque dark grey
  
  
  
  
  p <-  igraph::plot.igraph( g.net,
                             edge.color = edge_colors,
                             vertex.size = 6,
                             edge.curved = 0.25,
                             edge.color = "grey20" ,
                             vertex.color = V(g.net)$colour,
                             vertex.shape = V(g.net)$shape,
                             vertex.label = NA,
                             vertex.frame.colour = "grey20",
                             edge.width = E(g.net)$width*20,
                             frame = FALSE,
                             layout=l,
                             asp = 1,
                             #     mark.groups = list.D,
                             #    mark.border =NA,
                             #  mark.col=col.adj
                             edge.width = E(g.net)$weight)
  
return(p)  
}

# 
setwd("C:/Users/sonja/Desktop/Konstanz/Breeding Season 2021/Fledgie social learning/git/fledgie_learning/Output/Figures")

tiff("Figure_Networks_vertical.tiff", units="in", width=5, height=5, res=300, compression = 'lzw')
vert.plot <- plot.network(which.network = Mill.vert.matrix)
dev.off()

tiff("Figure_Networks_oblique.tiff", units="in", width=5, height=5, res=300, compression = 'lzw')
obl.plot <- plot.network(which.network = Mill.obl.matrix)
dev.off()

tiff("Figure_Networks_sibling.tiff", units="in", width=5, height=5, res=300, compression = 'lzw')
sib.plot <- plot.network(which.network = Mill.sib.matrix)
dev.off()

tiff("Figure_Networks_nonsibling.tiff", units="in", width=5, height=5, res=300, compression = 'lzw')
nonsib.plot <- plot.network(which.network = Mill.hor.matrix)
dev.off()



