# Script to reproduce the results of:
#####
# Title: Contribution from the Silent Majority to Dominate Dengue Transmission
# Authors: Quirine A. ten Bosch, Hannah E. Clapham, Louis Lambrechts, Veasna Duong,
# Philippe Buchy, Benjamin M. Althouse, Alun L. Lloyd, Lance A. Waller,
# Amy C. Morrison, Uriel Kitron, Gonzalo M. Vazquez-Prokopec, Thomas W. Scott,
# T. Alex Perkins
# Journal: PLOS Pathogens (in Review)
# Year: 2018
#####

rm(list=ls())

# Set work directory ------------------------------------------------------
setwd(dirname(list.files(pattern='Main_AsymptomaticsPaper.R', recursive=TRUE, full.names=TRUE)))

# Run source files -------------------------------------------------------
source('Main_ViremiaToInfectiousness.R')
source('Main_ProportionImmuneHistory.R')
source('Main_FOICalculations.R')
save.image(file = "Workspace.lowerAsym.RData")
source('Supp_FOICalculations_4infections.R')
save.image(file = "Workspace.lowerAsym.RData")
source('Supp_TernaryCalculations.R') 

# Compile Outcomes for paper ------------------------------------------
source('Main_DeriveOutcomes.R')
source('Main_CompileOutcomes.R')
source('Supp_NetInfectiousnessComparison.R')   
save.image(file = "Workspace.lowerAsym.RData")

# Create Final Figures for paper ------------------------------------------
source('Main_CoreFigures.R')
source('Supp_Figures.R')
source('Supp_TernaryPlot.R') 

# Uncertainty analysis ----------------------------------------------------

rm(list=ls())
elim.uncertainty = 'viremia' 
source('Main_ViremiaToInfectiousness_onelessuncertainty.R')
save.image(file = "Workspace.Uncertainty.viremia.RData")

rm(list=ls())
elim.uncertainty = 'ratios' 
source('Main_ViremiaToInfectiousness_onelessuncertainty.R')
save.image(file = "Workspace.Uncertainty.ratios.RData")

rm(list=ls())
elim.uncertainty = 'infectiousness' 
source('Main_ViremiaToInfectiousness_onelessuncertainty.R')
save.image(file = "Workspace.Uncertainty.infectiousness.RData")

rm(list=ls())
elim.uncertainty = 'IIP' 
source('Main_ViremiaToInfectiousness_onelessuncertainty.R')
save.image(file = "Workspace.Uncertainty.iip.RData")

source('Main_VarianceContributions.R')
