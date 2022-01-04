
# 1) Install and load libraries -------------------------------------------

# load NBDa via devtools
library(devtools)
# install_github("whoppitt/NBDA")
library(NBDA)


# 2) Data preparation -----------------------------------------------------


# 2.1) Load network data and puzzle data --------------------------------------------------

setwd("C:/Users/swild/Desktop/Konstanz/Breeding Season 2021/Fledgie social learning/Data")

load("Guett.gmm.RDA")
load("Mill.gmm.RDA")

load("Mill.puzzle.data.RDA")
load("Guett.puzzle.data.RDA")

# load puzzle data from previous winter
load("puzzle.data.winter20.RDA")



# 2.2. Load breeder data --------------------------------------------------

fledge.data <- read.delim("fledge.dates.txt", sep=" ")
breeders <- read.delim("Breeders 2021.txt")
PITs <- read.delim("Species_PITs.txt")

# Make a list of great tit chicks
GT.chicks <- subset(breeders$Tag, breeders$Who=="Chick" & breeders$Species=="GRETI" & breeders$Tag!="")
GT.chicks <- GT.chicks[!is.na(GT.chicks)]

GT.list <- subset(PITs$Pit, PITs$Species=="GRETI")

# 2.3. Extract independence day for fledgies ------------------------------
######
# extract the day where chick independent from parents
# i.e. the first day where the number of visits with parents is smaller or equal than the number of visits with other adults (when parents were not present)
# or if the chick has been seen alone at least once

`%notin%` <- Negate(`%in%`)


independent.chicks <- function(data, site){ 
  
  birds.visited <- colnames(data$gbi)
  chick.visited <- subset(birds.visited, birds.visited%in%GT.chicks)
  df.chick.independent <- as.data.frame(matrix(NA, nrow=length(chick.visited), ncol=7))
  colnames(df.chick.independent) <- c("chick_tag","fledge_date", "day_independent", "days_till_independence", "box", "site", "weight")
  df.chick.independent$chick_tag <- chick.visited
  days <- unique(substr(data$metadata$Location, nchar(data$metadata$Location[1])-5, nchar(data$metadata$Location[1])))
  data$metadata$days <- substr(data$metadata$Location, nchar(data$metadata$Location[1])-5, nchar(data$metadata$Location[1]))
  for (i in chick.visited){
    box <- subset(breeders$Box, breeders$Tag==i)
    site.i <- subset(breeders$Location, breeders$Tag==i)
    weight <- subset(breeders$Chick.weight, breeders$Tag==i)
    parents <- c(subset(breeders$Tag, breeders$Box==box & breeders$Who=="Female"),subset(breeders$Tag, breeders$Box==box & breeders$Who=="Male"))
    parents.seen <- subset(parents, parents%in%birds.visited) # parents who were actually seen at the feeder
    
    for (j in sort(unique(days), decreasing=T)){ # run it from the last  day to the first to get the earliest day where chicks were independent
      j.date <- as.Date(as.character(j), format=c("%y%m%d"))
      fledge.date <- as.Date(fledge.data[which(fledge.data$PIT==i),"Date.Time"], format="%Y-%m-%d")
      days.since.fledge <- j.date-fledge.date
      if(days.since.fledge<30){
        sub.data <-  data$gbi[as.numeric(rownames(data$metadata[data$metadata$days %in% j,])),] # subset the date j
        sub.data <-sub.data[which(sub.data[,i]>0),] # subset to groups that the chick i had been part of
      
        
        if(!is.null(nrow(sub.data))){
          alone <- length(subset(sub.data, rowSums(sub.data)==1)[,1])
          
          if(length(parents.seen)!=0){
            w.parents.sub <- as.matrix(sub.data[,parents.seen])
            w.parents.sub <- subset(w.parents.sub, rowSums(w.parents.sub)>0)
            w.parents <- length(w.parents.sub[,1])
         #   rows.w.parents <- rownames(w.parents.sub)
          } else {w.parents <- 0}
          
          # extract adult non-parents
          non.parents <- setdiff(birds.visited[birds.visited %notin% GT.chicks], c(chick.visited, parents))
          w.non.parents.sub <- sub.data[,non.parents] # extract columns with non parents
          w.non.parents <- length(subset(w.non.parents.sub, rowSums(w.non.parents.sub)>0)[,1])
          
          # if the number of visits with other adults is equal to or exceeds the number of visits with parents, or the chick has visited the feeder alone at least twice
          
          if(w.non.parents >= w.parents & sum(alone+w.non.parents+w.parents)!=0 | alone >=2) {
            df.chick.independent[which(chick.visited==i),2] <- as.character(format(fledge.date, format="%y%m%d"))
            df.chick.independent[which(chick.visited==i),3] <- j 
            df.chick.independent[which(chick.visited==i),4] <- as.numeric(days.since.fledge)
            df.chick.independent[which(chick.visited==i),5] <- box
            df.chick.independent[which(chick.visited==i),6] <- site.i
            df.chick.independent[which(chick.visited==i),7] <- weight
            
          }
        }
        
      }
      
    } 
    
  }
  df.chick.independent <- df.chick.independent[order(df.chick.independent$box), ]
  df.chick.independent <- df.chick.independent[df.chick.independent$site %in% site,]
  return(df.chick.independent) 
}

independence.Guett <- independent.chicks(data=Guett.gmm, site="Guett")
independence.Mill <- independent.chicks(data=Mill.gmm, site="Castle")


# 2.4. Make a list of IDs to include --------------------------------------

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

# 2.5. Create vertical social network -------------------------------------

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
    # SW: double check direction is correct
    mat[i, parents] <- 1
    
  }  
  
 return(mat)   
}

Mill.vert.matrix <- create.vert.net(Mill.IDs)
Guett.vert.matrix <- create.vert.net(Guett.IDs)


# 2.6. Create oblique social network --------------------------------------

create.obl.net <- function(IDs){
  mat <- matrix(0, nrow=length(IDs), ncol=length(IDs))
  rownames(mat) <- IDs
  colnames(mat) <- IDs
  
  # chick IDs
  IDs.sub <- IDs[IDs %in% breeders$Tag & IDs %in% GT.chicks]
  
  for(i in IDs.sub){
    # extract the info for the bird i
    i.box <- subset(breeders$Box, breeders$Tag==i)
    # extract the adult non-parents of that chick
    adults <- subset(breeders$Tag, breeders$Who %in% c("Male", "Females") & breeders$Box != i.box)
    adults <- subset(adults, adults %in% IDs)
    # set entries between non-parent adults and offspring to 1 in the direction of the chick learning from adults
    # but not the other way round
    # SW: double check direction is correct
    mat[i, adults] <- 1
    
  }  
  
  return(mat)   
}

Mill.obl.matrix <- create.obl.net(IDs = Mill.IDs)
Guett.obl.matrix <- create.obl.net(IDs = Guett.IDs)



# 2.7. Create horizontal social network -----------------------------------

create.hor.net <- function(IDs){
  mat <- matrix(0, nrow=length(IDs), ncol=length(IDs))
  rownames(mat) <- IDs
  colnames(mat) <- IDs
  
  # chick IDs
  IDs.sub <- IDs[IDs %in% breeders$Tag & IDs %in% GT.chicks]
  
  # set connections among all chicks to 1
  # they can learn in both directions
  mat[IDs.sub, IDs.sub] <- 1
  
  return(mat)   
}

Mill.hor.matrix <- create.hor.net(IDs = Mill.IDs)
Guett.hor.matrix <- create.hor.net(IDs = Guett.IDs)


# 2.8. Extract the demonstrators ----------------------------------------

# first we extract all the knowledgeable birds until the day that the first fledgling left the nest
min(fledge.data$Date.Time[fledge.data$Date.Time!=""])
# "2021-05-14 02:00:00" # the first fledgie fledged on the 14.05.2021
# everybody who knows how to solve up to this point (with at least 3 solves) is considered a demonstrator

# first we need to consider those who have learned during the winter experiment (20/21), but have not produced any solutions in spring up until this day
# we first subset the data to solves
puzzle.data.winter20.combined <- rbind(puzzle.data.winter20$Guett, puzzle.data.winter20$Mill)
puzzle.data.winter20.combined <- subset(puzzle.data.winter20.combined, puzzle.data.winter20.combined$event %in% c("right", "left"))
num.solves.winter20 <- as.data.frame(table(puzzle.data.winter20.combined$PIT))
solvers.winter20 <- as.vector(subset(num.solves.winter20$Var1, num.solves.winter20$Freq>=3))

# now we extract the knowledgeable birds from the breeding season 2021 up until the 14th of May 2021
puzzle.data.spring21.combined <- rbind(Guett.puzzle.data, Mill.puzzle.data)
puzzle.data.spring21.combined <- subset(puzzle.data.spring21.combined, puzzle.data.spring21.combined$Event %in% c("right", "left") & puzzle.data.spring21.combined$Date.Time < "210514000000")
num.solves.spring21 <- as.data.frame(table(puzzle.data.spring21.combined$PIT))
solvers.spring21 <- as.vector(subset(num.solves.spring21$Var1, num.solves.spring21$Freq>=3))

# combining the two
demos <- unique(c(solvers.winter20, solvers.spring21))
# these are all knowledgeable birds from Mill and Guett up until the 14th of May 2021


# 2.9. Extract date of acquisition for fledgies ---------------------------

extract.acquisition.date <- function(puzzle.data){
  sites <- sort(unique(puzzle.data$Location))
  # extract those with at least three solves
  
  puzzle.data.sub <- subset(puzzle.data, puzzle.data$Event %in% c("right", "left"))
  solvers <- as.data.frame(table(puzzle.data.sub$PIT))
  solvers <- as.vector(subset(solvers$Var1, solvers$Freq>=3))
  # remove demonstrators
  solvers <- subset(solvers, !(solvers %in% demos ))
  
  # no extract the date of the first solves at each of the puzzle sites (if applicable)
  # we consider the bird a learner at the site that it first produced the solution
  # in all other sites, it will be filtered from the diffusion (it can still transmit behaviour but it is not considered in the model)
  data.all <- NULL
  for(i in solvers){
  acq.dates <- NULL
      for(j in sites){
      sub.i.j <- subset(puzzle.data.sub, puzzle.data.sub$PIT==i & puzzle.data.sub$Location==j)
      if(length(sub.i.j$PIT)==0){
        date.first.solve.j <- NA
      } else {
        date.first.solve.j <- min(sub.i.j$Date.Time)  
      }
      acq.dates[which(sites==j)] <- date.first.solve.j
      }
  i.data <- c(i, acq.dates)
  data.all <- rbind(data.all, i.data)
  }
  data.all <- as.data.frame(data.all)
  rownames(data.all) <- NULL
  colnames(data.all) <- c("PIT", sites)
  return(data.all)
}


Mill.acquisition.date <- extract.acquisition.date(puzzle.data = Mill.puzzle.data)

Guett.acquisition.date <- extract.acquisition.date(puzzle.data = Guett.puzzle.data)


# 3) NBDA ------------------------------------------------------------------


# 3.1. Prepare NBDA data objects ------------------------------------------

prepare.NBDA.object <- function(networks, acquisition.dates, demos){
  
}

