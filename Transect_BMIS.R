## Information ##
# This is an example of a Best-Matched Internal Standard (BMIS) script for a HILICPos + HILICNeg QE run.
# Written by K.Heal, A.Boysen, and R.Lionheart

library(ggplot2)
library(stringr)
library(tidyverse)
options(scipen=999)


# Things to Return --------------------------------------------------------

# IS_inspectPlot (plot to make sure there aren't any internal standards we should kick out)
# QuickReport (% that picked BMIS, with cut off values)
# ISTest_plot (plot to evaluate if you cut off is appropriate)
# BMIS_normalizedData (tibble with the info you actually want!)

# Set parameters -----------------------------------------------------------------

cut.off <- 0.3 # 30% decrease in RSD of pooled injections, aka improvement cutoff
cut.off2 <- 0.1 # RSD minimum

# Imports -----------------------------------------------------------------

SampKey.all <- read.csv("data_extras/Sample.Key.EddyTransect.csv") %>%
  mutate(Sample.Name = Sample.Name %>%
           str_replace("-","."))
SampKey.all[SampKey.all == "180821_Poo_MesoScopeQC_1a"] <- "180821_Poo_MesoScopeQC_1"

Internal.Standards <- read.csv("data_extras/Ingalls_Lab_Standards.csv") %>%
  filter(Column == "HILIC") %>%
  filter(Compound.Type == "Internal Standard")

trimws(Internal.Standards$Compound.Name, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")


# HILICPos and HILICNeg data
transect <- read.csv("data_processed/QC_Output.csv", header = TRUE) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value))


# Change class + adjust data. Set cutoff values -----------------------------------------------------------------
transect <- transect %>%
  filter(!str_detect(ReplicateName, "Blk")) %>%
  filter(!str_detect(ReplicateName, "Std")) %>%
  mutate(ReplicateName = as.character(ReplicateName)) %>%
  mutate(Metabolite.name = as.character(Metabolite.name)) %>%
  mutate(RTValue = as.numeric(RTValue)) %>%
  mutate(AreaValue = as.numeric(AreaValue)) %>%
  mutate(SNValue = as.numeric(SNValue)) %>%
  filter(!(Column == "HILICNeg" & Metabolite.name == "Inosine")) %>%
  filter(!(Column == "HILICNeg" & Metabolite.name == "Guanine")) %>%
  mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name))


# Match transect data with Internal Standards list -----------------------------------------------------------------

transect.withIS <- transect %>%
  filter(Metabolite.name %in% Internal.Standards$Compound.Name)

transect.NoIS <- transect %>%
  filter(!Metabolite.name %in% Internal.Standards$Compound.Name)


# Read in Internal Standard data -----------------------------------------------------------------

transect.IS.data <- transect.withIS %>%
  select(ReplicateName, Metabolite.name, Area.with.QC) %>%
  mutate(MassFeature = Metabolite.name) %>%
  select(-Metabolite.name) %>%
  filter(!MassFeature == "Guanosine Monophosphate, 15N5")

# Drop syntactically correct "X" at start of ReplicateName.
transect.IS.data$ReplicateName <- gsub("^.{0,1}", "", transect.IS.data$ReplicateName)

# Add injection volume -----------------------------------------------------------------

transect.SampKey <- SampKey.all %>%
  filter(Sample.Name %in% transect.IS.data$ReplicateName) %>% 
  select(Sample.Name, Bio.Normalization) %>%
  mutate(MassFeature = "Inj_vol",
         Area.with.QC = Bio.Normalization,
         ReplicateName = Sample.Name) %>%
  select(ReplicateName, Area.with.QC, MassFeature)

# Create Internal standard data identify problematic compounds/replicates-----------------------------------------------------------------

transect.IS.data <- rbind(transect.IS.data, transect.SampKey)
# THIS REMOVAL OF DDA SAMPLES IS ADDED AS A STOPGAP MEASURE- NEEDS TO BE FIXED!!! ##
transect.IS.data <- transect.IS.data %>%
  filter(!grepl("DDA", ReplicateName))
# THIS REMOVAL OF DDA SAMPLES IS ADDED AS A STOPGAP MEASURE- NEEDS TO BE FIXED!!! ##

# Identify internal standards without an Area, i.e. any NA values.
IS_Issues <- transect.IS.data[is.na(transect.IS.data$Area.with.QC),]


# Visualize raw areas of Internal Standards -----------------------------------------------------------------

IS_inspectPlot <- ggplot(transect.IS.data, aes(x = ReplicateName, y = Area.with.QC)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~MassFeature, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 5),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10))+
  ggtitle("IS Raw Areas")
#print(IS_inspectPlot)


# Edit data so names match-----------------------------------------------------------------
transect.IS.data <- transect.IS.data %>%
  mutate(ReplicateName = ReplicateName %>%
           str_replace("-",".")) %>%
  arrange(ReplicateName)

transect.long  <- transect.NoIS %>%
  rename(MassFeature = Metabolite.name) %>%
  select(ReplicateName, MassFeature, Area.with.QC) %>%
  arrange(ReplicateName)

# Drop syntactically valid "X" from ReplicateName.
transect.long$ReplicateName <- gsub("^.{0,1}", "", transect.long$ReplicateName)

# Test that names are equal
test_IS.data <- as.data.frame(unique(transect.IS.data$ReplicateName))
test_IS.data <- test_IS.data %>% 
  mutate(ReplicateName = as.character(unique(transect.IS.data$ReplicateName))) 

test_long.data <- as.data.frame(unique(transect.long$ReplicateName))
test_long.data <- test_long.data %>% 
  mutate(ReplicateName = as.character(unique(transect.long$ReplicateName))) %>%
  filter(!grepl("DDA", ReplicateName)) 

identical(test_IS.data$ReplicateName, test_long.data$ReplicateName)

# Caluclate mean values for each Internal Standard----------------------------------------------------------------
transect.IS.means <- transect.IS.data %>%
  filter(!grepl("_Blk_", ReplicateName)) %>%
  mutate(MassFeature = as.factor(MassFeature)) %>%
  group_by(MassFeature) %>%
  summarise(Average.Area = mean(as.numeric(Area.with.QC), na.rm = TRUE)) %>%
  mutate(MassFeature = as.character(MassFeature))

transect.IS.means[is.na(transect.IS.means)] <- NA


# Normalize to each internal Standard----------------------------------------------------------------
transect.binded <- rbind(transect.IS.data, transect.long) %>%
  arrange(MassFeature)

Split_Dat <- list()

for (i in 1:length(unique(transect.IS.data$MassFeature))) {
  Split_Dat[[i]] <- transect.binded %>%
    mutate(MIS = unique(transect.IS.data$MassFeature)[i]) %>%
    left_join(transect.IS.data %>%
                rename(MIS = MassFeature, IS_Area = Area.with.QC) %>%
                select(MIS, ReplicateName, IS_Area), by = c("ReplicateName", "MIS")) %>%
    left_join(transect.IS.means %>%
                rename(MIS = MassFeature), by = "MIS") %>%
    mutate(Adjusted_Area = Area.with.QC/IS_Area*Average.Area)
}

transect.area.norm <- do.call(rbind, Split_Dat) %>%
  select(-IS_Area, -Average.Area)

# Standardize name structure to: Date_type_ID_replicate_anythingextra) ----------------------------------------------------------------
transect.mydata_new <- transect.area.norm %>%
  separate(ReplicateName, c("runDate", "type", "SampID", "replicate"), "_") %>%
  mutate(Run.Cmpd = paste(transect.area.norm$ReplicateName, transect.area.norm$MassFeature))


# Find the B-MIS for each MassFeature----------------------------------------------------------------

# Look only at the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo),
# then choose which IS reduces the RSD the most (Poo.Picked.IS)
transect.poodat <- transect.mydata_new %>%
  filter(type == "Poo") %>%
  group_by(SampID, MassFeature, MIS) %>%
  summarise(RSD_ofPoo_IND = sd(Adjusted_Area, na.rm = TRUE) / mean(Adjusted_Area, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>%
  group_by(MassFeature, MIS) %>%
  summarise(RSD_ofPoo =  mean(RSD_ofPoo_IND, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo = ifelse(RSD_ofPoo == "NaN", NA, RSD_ofPoo)) 


transect.poodat <- transect.poodat %>%
  left_join(transect.poodat %>% group_by(MassFeature) %>%
              summarise(Poo.Picked.IS = unique(MIS)[which.min(RSD_ofPoo)] [1]))


# Get the original RSD, calculate RSD change, decide if MIS is acceptable----------------------------------------------------------------
transect.poodat <- left_join(transect.poodat, transect.poodat %>%
                                   filter(MIS == "Inj_vol" ) %>%
                                   mutate(Orig_RSD = RSD_ofPoo) %>%
                                   select(-RSD_ofPoo, -MIS)) %>%
  mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
  mutate(percentChange = del_RSD/Orig_RSD) %>%
  mutate(accept_MIS = (percentChange > cut.off & Orig_RSD > cut.off2))


# Change the BMIS to "Inj_vol" if the BMIS is not an acceptable----------------------------------------------------------------

# Adds a column that has the BMIS, not just Poo.Picked.IS
# Changes the FinalBMIS to inject_volume if its no good

transect.fixedpoodat <- transect.poodat %>%
  filter(MIS == Poo.Picked.IS) %>% 
  mutate(FinalBMIS = ifelse(accept_MIS == "FALSE", "Inj_vol", Poo.Picked.IS)) %>%
  mutate(FinalRSD = RSD_ofPoo)

newpoodat <- transect.poodat %>%
  left_join(transect.fixedpoodat %>% select(MassFeature, FinalBMIS)) %>%
  filter(MIS == FinalBMIS) %>%
  mutate(FinalRSD = RSD_ofPoo)

Try <- newpoodat %>%
  filter(FinalBMIS != "Inj_vol")

QuickReport <- print(paste("% of MFs that picked a BMIS",
                           length(Try$MassFeature) / length(newpoodat$MassFeature),
                           "RSD improvement cutoff", cut.off,
                           "RSD minimum cutoff", cut.off2,
                           sep = " "))


# Evaluate and visualize the results of your BMIS cutoff----------------------------------------------------------------
IS_toISdat <- transect.mydata_new %>%
  filter(MassFeature %in% transect.IS.data$MassFeature) %>%
  select(MassFeature, MIS, Adjusted_Area, type) %>%
  filter(type == "Smp") %>%
  group_by(MassFeature, MIS) %>%
  summarise(RSD_ofSmp = sd(Adjusted_Area, na.rm = TRUE)/mean(Adjusted_Area, na.rm = TRUE)) %>%
  left_join(transect.poodat %>% select(MassFeature, MIS, RSD_ofPoo, accept_MIS))

injectONlY_toPlot <- IS_toISdat %>%
  filter(MIS == "Inj_vol")


ISTest_plot <- ggplot() +
  geom_point(dat = IS_toISdat, shape = 21, color = "black", size = 2,aes(x = RSD_ofPoo, y = RSD_ofSmp, fill = accept_MIS)) +
  scale_fill_manual(values=c("white","dark gray")) +
  geom_point(dat = injectONlY_toPlot, aes(x = RSD_ofPoo, y = RSD_ofSmp), size = 3) +
  facet_wrap(~ MassFeature)
#print(ISTest_plot)


# Return data that is normalized via BMIS----------------------------------------------------------------

## original
transect.BMIS_normalizedData <- newpoodat %>% select(MassFeature, FinalBMIS, Orig_RSD, FinalRSD) %>%
  left_join(transect.mydata_new, by = "MassFeature") %>%
  filter(MIS == FinalBMIS) %>%
  unique()

write.csv(transect.BMIS_normalizedData, file = "data_processed/BMIS_Output.csv")

rm(list = ls())

