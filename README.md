SilentMajority_2018
====================

This repository contains code used in the following paper.

Quirine A. ten Bosch, Hannah E. Clapham, Louis Lambrechts, Veasna Duong, Philippe Buchy, Benjamin M. Althouse, Alun L. Lloyd, Lance A. Waller,Amy C. Morrison, Uriel Kitron, Gonzalo M. Vazquez-Prokopec, Thomas W. Scott,T. Alex Perkins (2018)
**Contributions from the silent majority dominate dengue virus transmission**. *PLOS Pathogens* 

All code contained within this repository is released under the [CRAPL v0.1 License](http://matt.might.net/articles/crapl/).  Because of the large sizes of the output files, we have deposited those on Open Science Framework [https://osf.io/pjbhz/](https://osf.io/pjbhz/). Here you can also find data used in the analysis. 

====================

### Set up of code base: 

* `Scripts` contains code to run analyses and reproduce figures in the manuscript and supplements. Files starting with Main or Supp generate the analyses and figures for the manuscript and supplements respectively. Function-files contain models and functions for the Main and Supp files with the same extension.   
* `Data` contains population data used in `Scripts`

### Scripts folder

The scripts in this folder are used to run the several steps of the analysis described in the manuscript. 

To run the full set of analyses:

* `Main.Asymptomatics.Paper.R`, which calls the files below

to run the main analyses: 
* `Main_ViremiaToInfectiousness.R` 
* `Main_ProportionImmuneHistory.R`
* `Main_FOICalculations.R`

to run the supplementary analyses: 
* `Supp_FOICalculations_4infections.R` 
* `Supp_TernaryCalculations.R`

to derive and compile outcomes as presented in the manuscript:
* `Main_DeriveOutcomes.R`
* `Main_CompileOutcomes.R`
* `Supp_NetInfectiousnessComparison.R`

to create figures for the paper:
* `Main_CoreFigures.R`
* `Supp_Figures.R`
* `Supp_TernaryPlot.R`

The above code uses dose response curves generated with:
* `Supp_DoseResponseCurves.R`

And results from a meta-analysis on the proportion of apparent infections by immune history:
* `Main_Meta_PrimSec.R`
* `Main_Meta_PrimSec_BB.R`
* `Main_Meta_PostSec.R`

### Data folder

Contains age specific population data for Thailand and Brazil. For all other data, please refer to Open Science Framework project [https://osf.io/pjbhz/](https://osf.io/pjbhz/).


