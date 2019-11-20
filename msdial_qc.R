## Information ##
# This is an example of a quality control filter for a HILICPos + HILICNeg QE run.
# Written by K.Heal and R.Lionheart

library(rlist)
library(tidyverse)

# Functions ---------------------------------------------------------------

IdentifyRunTypes <- function(msdial.file) {
  # Identify run typfes and return each unique value present in the Skyline output.
  #
  # Args
  #   msdial.file: Raw output file from Skyline.
  #
  # Returns
  #   run.type: list of labels identifying the run types, isolated from Replicate.Name.
  #   Options conssist of samples (smp), pooled (poo), standards (std), and blanks (blk).
  #
  run.type <- tolower(str_extract(msdial.file$ReplicateName, "(?<=_)[^_]+(?=_)"))
  print(paste("Your runtypes are:", toString(unique(run.type))))
}

# Define parameters -------------------------------------------------------

area.min   <- 1000
RT.flex    <- 0.4
blk.thresh <- 0.3
SN.min     <- 4

# Identify runtypes and select columns -------------------------------------------------------

msdial.runtypes <- IdentifyRunTypes(combined)
combined <- combined %>%
  select(ReplicateName:Alignment.ID, Metabolite.name) %>%
  mutate(Run.Type = (tolower(str_extract(combined$ReplicateName, "(?<=_)[^_]+(?=_)")))) 

# Create retention time comparison table -------------------------------------------------------

RT.table <- combined %>%
  filter(Run.Type == "std") %>%
  arrange(Metabolite.name) %>%
  group_by(Metabolite.name) %>%
  mutate(RT.min = min(RTValue, na.rm = TRUE)) %>%
  mutate(RT.max = max(RTValue, na.rm = TRUE)) %>%
  select(Metabolite.name:RT.max) %>%
  unique()

# Create blanks comparison table -------------------------------------------------------

blank.table <- combined %>%
  filter(Run.Type == "blk") %>%
  mutate(Blk.Area = AreaValue) %>%
  arrange(Metabolite.name) %>%
  group_by(Metabolite.name) %>%
  mutate(Blk.min = min(AreaValue)) %>%
  mutate(Blk.max = max(AreaValue)) %>%
  select(Metabolite.name:Blk.max) %>%
  select(-Blk.Area) %>%
  unique()

# Add Signal to Noise and Area flags -------------------------------------------------------

SN.Area.Flags <- combined %>%
  arrange(Metabolite.name) %>%
  mutate(SN.Flag       = ifelse(((SNValue) < SN.min), "SN.Flag", NA)) %>%
  mutate(area.min.Flag = ifelse((AreaValue < area.min), "area.min.Flag", NA))

# Add Retention time flags -------------------------------------------------------

add.RT.Flag <- SN.Area.Flags %>%
  group_by(Metabolite.name) %>%
  left_join(RT.table, by = c("Metabolite.name", "Run.Type")) %>%
  mutate(RT.Flag = ifelse((RTValue >= (RT.max + RT.flex) | RTValue <= (RT.min - RT.flex)), "RT.Flag", NA)) %>%
  select(-c("RT.max", "RT.min"))

# Add Blank flags -------------------------------------------------------

add.blk.Flag <- add.RT.Flag %>%
  left_join(blank.table, by = c("Metabolite.name", "Run.Type")) %>%
  mutate(blank.Flag = ifelse((AreaValue / Blk.max) < blk.thresh, "blank.Flag", NA)) %>%
  select(-c("Blk.min", "Blk.max"))

# Combine all flags and  standardize Metabolite.name variable -------------------------------------------------------

final.table <- add.blk.Flag %>%
  mutate(all.Flags      = paste(SN.Flag, area.min.Flag, RT.Flag, blank.Flag, sep = ", ")) %>%
  mutate(all.Flags      = as.character(all.Flags %>% str_remove_all("NA, ") %>% str_remove_all("NA"))) %>%
  mutate(Area.with.QC   = ifelse(is.na(area.min.Flag), AreaValue, NA)) %>%
  select(ReplicateName:AreaValue, Area.with.QC, everything()) %>%
  ungroup(Metabolite.name) %>%
  mutate(Metabolite.name = as.character(Metabolite.name)) %>%
  mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name))


# Add parameter comments and print to file -------------------------------------------------------

Description <- c("Hello! Welcome to the world of MSDIAL QE Quality Control! ",
                 "Minimum area for a real peak: ",
                 "RT flexibility: ",
                 "Blank can be this fraction of a sample: ",
                 "S/N ratio: " ,
                 "Processed on: ")
Value <- c(NA, area.min, RT.flex, blk.thresh, SN.min, Sys.time())

df <- data.frame(Description, Value)
final.table <- bind_rows(df, final.table)

rm(list=setdiff(ls(), "final.table"))

write.csv(final.table, "./data_processed/QC_Output.csv", row.names = FALSE)

