
# 1) Install and load libraries -------------------------------------------

# NBDA



# 2) Data preparation -----------------------------------------------------


# 2.1) Load network data --------------------------------------------------

setwd("C:/Users/swild/Desktop/Konstanz/Breeding Season 2021/Fledgie social learning/Data")

load("Guett.gmm.RDA")
load("Mill.gmm.RDA")


# 2.2. Load breeder data --------------------------------------------------

fledge.data <- read.delim("fledge.dates.txt", sep=" ")
breeders <- read.delim("Breeders 2021.txt")

# Make a list of great tit chicks
GT.chicks <- subset(breeders$Tag, breeders$Who=="Chick" & breeders$Species=="GRETI" & breeders$Tag!="")
GT.chicks <- GT.chicks[!is.na(GT.chicks)]


# 2.3. Extract independence day for fledgies ------------------------------

# read in GBI data

GBI_Mill_1_13 <- read.delim("week1-13.gbi.Mill.txt", sep=",", check.names = FALSE)
head(GBI_Mill_1_13)


# read in chick list
chicks.df <- read.delim("greti.chicks.2020.txt", sep="\t")
chicks <- as.character(na.omit(subset(chicks.df$Tag, subset=chicks.df$Who=="Chick")))

######
# extract the day where chick independent from parents
# i.e. the first day where the number of visits with parents is smaller or equal than the number of visits with other adults (when parents were not present)
# or if the chick has been seen alone at least once

`%notin%` <- Negate(`%in%`)

data <- Guett.gmm

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


