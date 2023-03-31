##### Fungi Perfecti ####
##### March 31, 2023 ####
##### Version: Public GitHub version    ####

### Objective:  Analyze Cryo-library Culture Inventory###

#Load packages
library(tidyverse)
library(lubridate)

#Import Indexed Data
rawCryo <- as.data.frame(read_csv("Cryo_In -- Public.csv")) 
#Import Cryo_Out Data
rawCryoOut <- as.data.frame(read_csv("Cryo_Out -- Public.csv")) 

##Column Rename Function!
### Takes df and desired separator as arguments.
### Iterates through all columns, replaces spaces with chosen separator.
### i.e "Column 1" to "Column.1" or "Column_1", etc.
### default sep = "."
col_rename <- function (df, sep=".") {
  s <- unique(colnames(df));
  #print(s)
  newNames <- c()
  
  for (i in seq_along(s)) {
    newName <- gsub(" ", sep, s[i]);
    newNames <- c(newNames, newName)
    #print(newName)
  }
  
  colnames(df) <- newNames
  #print(newNames)
  return(df)
  
}

#Reshape Cryo In Data
cryo <- col_rename(rawCryo)%>% #rename columns
  #lubridate
  mutate(Inoc.Date = mdy(Inoc.Date),
         Date.Frozen = mdy(Date.Frozen))

#Reshape Cryo Out data
cryoOut <- col_rename(rawCryoOut)%>% #rename columns
  #lubridate
  mutate(Inoc.Date = mdy(Inoc.Date),
         Date.Frozen = mdy(Date.Frozen),
         Date.Thawed = mdy(Date.Thawed),
         #Calculate Frozen Duration, and Viability for tubes that have been removed from LN2
         Days.Frozen.when.Thawed = Date.Thawed - Date.Frozen,
         Viability = Viable/n.Plates * 100,
         Grade = if_else(Viability > 2/3 *100, "Pass", "Fail"))

#Calculate Number of Tubes Thawed and % Viability
##for Every Unique Species, Strain, and Pvalue (irrespective of Thaw Date)
nThawed <- cryoOut%>%
  group_by(Accession.Number, Genus, Species, Strain, Inoc.Date, Pvalue)%>% 
  summarise(n.Tubes.Thawed = sum(n.Tubes.Thawed),
            Viability = sum(Viable, na.rm = T)/sum(n.Plates) * 100)

#Add Thawed Tubes and Viability Data to Cryo DF
##Calculate Number of Tubes Currently in LN2
##Calculate Incubation Time and Frozen Duration
cryo <- cryo%>%
  left_join(nThawed)%>% 
  mutate(n.Tubes = n.Tubes.Frozen-n.Tubes.Thawed,
         Incubation = Date.Frozen - Inoc.Date,
         Days.Frozen = today() - Date.Frozen)%>%
  filter(is.na(n.Tubes) |
           n.Tubes != 0)  

#Summarize Cryo Out Data
cryoOutSumm <- cryoOut%>%
  group_by(Accession.Number, Genus, Species, Strain, Inoc.Date, Pvalue)%>%
  summarise(Viability.Tested.Max.Days.Frozen = max(Days.Frozen.when.Thawed, 
                                                   na.rm = T),
            n.Tubes.Thawed = sum(n.Tubes.Thawed, na.rm = T),
            n.Plates = sum(n.Plates, na.rm = T),
            Viability = sum(Viable, na.rm = T)/sum(n.Plates) * 100)

#Genus Summ Table
genusTable <- cryo%>% 
  filter(!is.na(Viability))%>%
  ungroup()%>%
  group_by(Genus)%>%
  summarise(Species = n_distinct(Species),
            Strains = n_distinct(Strain),
            Batches = n_distinct(paste(Species, Strain, Inoc.Date, sep = ' ')),
            Samples = sum(n.Tubes, na.rm = T))
#write.csv(genusTable, "Cryo Genus Summary_v2.0.csv", row.names = F)

speciesTable <- cryo%>% 
  filter(!is.na(Viability))%>%
  ungroup()%>%
  group_by(Genus, Species)%>%
  summarise(Strains = n_distinct(Strain),
            Batches = n_distinct(paste(Species, Strain, Inoc.Date, sep = ' ')),
            Samples = sum(n.Tubes, na.rm = T))
#write.csv(speciesTable, "Cryo Species Summary_v2.0.csv", row.names = F)
