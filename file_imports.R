## Information ##
# This is an example of a simple MSDial import for HILICPos + HILICNeg runs from the QE.
# Written by K.Heal and R.Lionheart

library(tidyverse)
library(tidyr)
options(scipen=999)

# Functions ---------------------------------------------------------------

header.true <- function(df) {
  colnames(df) <- make.names(as.character(unlist(df[1,])))
  df <- df[-1, ]
  
  return(df)
}

filter.unknowns <- function(df) {
  df <- df %>%  
    filter(Metabolite.name != 'Unknown') %>%
    select(-c(Formula, Ontology, INCHIKEY, SMILES, Isotope.tracking.parent.ID, Isotope.tracking.weight.number, MS1.isotopic.spectrum, MS.MS.spectrum)) %>%
    mutate(Ordered.ID = 1:n()) %>%
    select(Ordered.ID, everything())
}

remove.csv <- function(full.names) {
  no.path <- substr(full.names, 1, nchar(full.names)-4)
  no.ID <-   gsub("\\_.*","", no.path)
  
  return(no.path)
}


# File uploads ------------------------------------------------------------

filenames <- list.files(path = 'data_raw', pattern = '*.csv')

names <- remove.csv(filenames)

for(i in names) {
  filepath <- file.path('data_raw', paste(i,".csv", sep = ""))
  assign(i, read.csv(filepath, skip = 3))
}



# Trim, set header, filter unknowns ---------------------------------------

SN.pos <- header.true(SN_HILICPos_EddyTransect)
SN.pos <- filter.unknowns(SN.pos)

RT.pos <- header.true(RT_HILICPos_EddyTransect)
RT.pos <- filter.unknowns(RT.pos)

Area.pos <- header.true(Area_HILICPos_EddyTransect)
Area.pos <- filter.unknowns(Area.pos)

MZ.pos <- header.true(Mz_HILICPos_EddyTransect)
MZ.pos <- filter.unknowns(MZ.pos)

SN.neg <- header.true(SN_HILICNeg_EddyTransect)
SN.neg <- filter.unknowns(SN.neg)

RT.neg <- header.true(RT_HILICNeg_EddyTransect)
RT.neg <- filter.unknowns(RT.neg)

Area.neg <- header.true(Area_HILICNeg_EddyTransect)
Area.neg <- filter.unknowns(Area.neg)

MZ.neg <- header.true(Mz_HILICNeg_EddyTransect)
MZ.neg <- filter.unknowns(MZ.neg)

## Add parameter tags and index columns
SN.pos <- SN.pos %>%
  select(Ordered.ID, -Post.curation.result, everything())

RT.pos <- RT.pos %>%
  select(Ordered.ID, -Post.curation.result, everything())

Area.pos <- Area.pos %>%
  select(Ordered.ID, -Post.curation.result, everything())

MZ.pos <- MZ.pos %>%
  select(Ordered.ID, -Post.curation.result, everything())

# Negative
SN.neg <- SN.neg %>%
  select(Ordered.ID, -Post.curation.result, everything())

RT.neg <- RT.neg %>%
  select(Ordered.ID, -Post.curation.result, everything())

Area.neg <- Area.neg %>%
  select(Ordered.ID, -Post.curation.result, everything())

MZ.neg <- MZ.neg %>%
  select(Ordered.ID, -Post.curation.result, everything())


# Change variable classes -------------------------------------------------

# SN.pos
for (i in c(25:ncol(SN.pos))) {
  SN.pos[, i] <- as.numeric(as.character(SN.pos[, i]))
}

afterSN.pos <- lapply(SN.pos, class)

# RT.pos
for (i in c(25:ncol(RT.pos))) {
  RT.pos[, i] <- as.numeric(as.character(RT.pos[, i]))
}

afterRT.pos <- lapply(RT.pos, class)

# Area.pos
for (i in c(25:ncol(Area.pos))) {
  Area.pos[, i] <- as.numeric(as.character(Area.pos[, i]))
}

afterArea.pos <- lapply(Area.pos, class)

# MZ.pos
for (i in c(25:ncol(MZ.pos))) {
  MZ.pos[, i] <- as.numeric(as.character(MZ.pos[, i]))
}

afterMZ.pos <- lapply(MZ.pos, class)


# SN.neg
for (i in c(25:ncol(SN.neg))) {
  SN.neg[, i] <- as.numeric(as.character(SN.neg[, i]))
}

afterSN.neg <- lapply(SN.neg, class)


# RT.neg
for (i in c(25:ncol(RT.neg))) {
  RT.neg[, i] <- as.numeric(as.character(RT.neg[, i]))
}

afterRT.neg <- lapply(RT.neg, class)

# Area.neg
for (i in c(25:ncol(Area.neg))) {
  Area.neg[, i] <- as.numeric(as.character(Area.neg[, i]))
}

afterArea.neg <- lapply(Area.neg, class)

# MZ.neg
for (i in c(25:ncol(MZ.neg))) {
  MZ.neg[, i] <- as.numeric(as.character(MZ.neg[, i]))
}

afterMZ.neg<- lapply(MZ.neg, class)

# Rearrange datasets ------------------------------------------------------
testSN.pos <- SN.pos %>%
  tidyr::gather(
    key = "ReplicateName",
    value = "SNValue",
    starts_with("X")) %>%
  select(ReplicateName, SNValue, everything())

testRT.pos <- RT.pos %>%
  tidyr::gather(
    key = "ReplicateName",
    value = "RTValue",
    starts_with("X")) %>%
  select(ReplicateName, RTValue, everything())

testArea.pos <- Area.pos %>%
  tidyr::gather(
    key = "ReplicateName",
    value = "AreaValue",
    starts_with("X")) %>%
  select(ReplicateName, AreaValue, everything())

testMZ.pos <- MZ.pos %>%
  tidyr::gather(
    key = "ReplicateName",
    value = "MZValue",
    starts_with("X")) %>%
  select(ReplicateName, MZValue, everything())


testSN.neg <- SN.neg %>%
  tidyr::gather(
    key = "ReplicateName",
    value = "SNValue",
    starts_with("X")) %>%
  select(ReplicateName, SNValue, everything())

testRT.neg <- RT.neg %>%
  tidyr::gather(
    key = "ReplicateName",
    value = "RTValue",
    starts_with("X")) %>%
  select(ReplicateName, RTValue, everything())

testArea.neg <- Area.neg %>%
  tidyr::gather(
    key = "ReplicateName",
    value = "AreaValue",
    starts_with("X")) %>%
  select(ReplicateName, AreaValue, everything())

testMZ.neg <- MZ.neg %>%
  tidyr::gather(
    key = "ReplicateName",
    value = "MZValue",
    starts_with("X")) %>%
  select(ReplicateName, MZValue, everything())

# Combine to one dataset --------------------------------------------------
combined.pos <- testArea.pos %>%
  left_join(testMZ.pos) %>%
  left_join(testSN.pos) %>%
  left_join(testRT.pos) %>%
  mutate(Column = "HILICPos") %>%
  select(ReplicateName, Ordered.ID, Column, AreaValue, MZValue, RTValue, SNValue, everything())

combined.neg <- testArea.neg %>%
  left_join(testMZ.neg) %>%
  left_join(testSN.neg) %>%
  left_join(testRT.neg) %>%
  mutate(Column = "HILICNeg") %>%
  select(ReplicateName, Ordered.ID, Column, AreaValue, MZValue, RTValue, SNValue, everything())

combined <- rbind(combined.pos, combined.neg)

# This is a comment

rm(list=setdiff(ls(), "combined"))
