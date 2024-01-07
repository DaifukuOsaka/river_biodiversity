# river_biodiversity
R code files supporting Qing Fu for his MSc thesis project "Model-based evaluation of optimal strategies embedding diverse datasets to assess biodiversity in rivers".
The most recent draft version of the thesis is available [here](https://www.overleaf.com/read/hbcbwdvfpskb#aac673). Feel free to [contact me](mailto:qing.fu@uzh.ch?subject=[Transparent%20Master%20Thesis]%20Inquiry) with questions, suggestions, etc.

## Navigating this Repo
This repo contains all the resources to replicate and analyze the simulation study and application with the exception of data from Castanho et al.'s (2020) paper, which can be downloaded from the authors' own repo [here](https://github.com/bcastanho/PRQ2019). All relevant scripts are included in the folder [/Rscripts](https://github.com/pitrieger/masterthesis/tree/main/Rscripts). The implementations of the detection methods for single- and multi-factor models are contained in the sub-folders [/simulation](https://github.com/pitrieger/masterthesis/tree/main/Rscripts/simulation) and [/application](https://github.com/pitrieger/masterthesis/tree/main/Rscripts/application), respectively, which corresponds with where they are used in the paper. Within these folders, the main scripts for running and analyzing the simulation study are [Simulation_final.R](https://github.com/pitrieger/masterthesis/blob/main/Rscripts/simulation/Simulation_final.R) and [Analysis_final.R](https://github.com/pitrieger/masterthesis/blob/main/Rscripts/simulation/Analysis_final.R), respectively. The main script for running the application is [application_main.R](https://github.com/pitrieger/masterthesis/blob/main/Rscripts/application/application_main.R). Simulation data is available in [/publicdata](https://github.com/pitrieger/masterthesis/blob/main/publicdata).

## Current work
Doing some final polishing and proofreading.

## Open Questions


## Log
- 16/08/2021 - Finished initial proposal
- 15/09/2021 - Registered thesis with study administration
- 05/10/2021 - Officially started writing
- 11/2021 - Finished draft for introduction of (confirmatory) factor analysis
- 11/2021 - Finished draft for introduction of measurement invariance
- 11/2021 - Implemented detection methods (existing and original contribution) for single-factor models
- 11/2021 - Implemented preliminary simulation study
- 12/2021 - Implemented detection methods for more general CFA models
- 12/2021 - Finished draft for introduction of detection methods
- 12/2021 - Finished draft for simulation study
- 01/2022 - Finished draft introduction
- 01/2022 - Finished draft for application
- 01/2022 - Programmed Shiny app for interactive version of simulation study
- 01/2022 - Formalized idea of novel method (to be added to draft)
- 02/2022 - Ran final version of simulation study
- 02/2022 - Added formalized novel method introduction
