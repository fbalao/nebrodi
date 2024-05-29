#### This is the master script for the analysis of Abies nebrodensis manuscript
library(here)

# 1. SNP panel design ----
source(here("scripts","snparraydesign","nebrodi_SNPpanel_selection.R"))

# 2. Genetic diversity comparison ----
source(here("scripts","RADseq_comparative","abies_comparativeRADseq.R"))

# 3. SNP panel Genotyping ----
source(here("scripts","Genotyping","Analisis_OpenArray_nebrodi_adults_seedlings.R"))
source(here("scripts","Genotyping","Cleaning_genotypes_mothers.R"))
source(here("scripts","Genotyping","Aneb_adults_seedlings_colony.R"))

# 4. Adults genetic diversity and structure ----
source(here("scripts","Diversity","DiversityAnalysis_Aneb.R"))

# 5. Nursery seedlings genotyping -----
source(here("scripts","Genotyping","Analisis_OpenArray_nebrodi_nursery_seedlings.R"))

# 6. Hybrids analysis ----
source(here("scripts","hybrids","hybridPCA.R"))
source(here("scripts","hybrids","hybridindex.R"))
source(here("scripts","hybrids","hybridindexnursery.R"))
       
# 7. AGF ----
source(here("scripts","AGF","AGF.R"))