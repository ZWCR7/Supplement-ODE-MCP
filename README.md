# Code-and-Data-for-``Structural Change Detection in Dynamic Systems''
Supplementary codes and data for ``Structural Change Detection in Dynamic Systems''

## Data and software availability
- The Italian COVID-19 dataset is available from Johns Hopkins University Center for Systems Science and Engineering (JHU CCSE) Coronavirus (https://systems.jhu.edu/research/public-health/ncov/) and the dataset can be derived through the R package ``coronavirus''. 
- The global temperature reconstruction dataset for the Holocene (the past 11,300 years) is compiled by [1].

## Code and Data for ``Structural Change Detection in Dynamic Systems''
### Overview

1. Directory *utils* contains the main functions for conducting structural change detection in dynamic systems. Specifically,
- The R script ***seededInterval*** contains the code for generating seeded intervals.
- The R script ***residualFunc*** contains the code for calculating the target optimization function.
- The R script ***StatGenerator*** contains the code for computing the test statistics, including the GLR statistic.
- The R script ***SNOT-Algorithm*** contains the code for conducting SNOT algorithm with different test statistics.
- The R script ***mopsFunc*** contains the code for the SPC FDR control procedure.

2. Directory *Simulations* contains the main functions for conducting simulations. Specifically,
- The R script ***Simulation-FN-General*** contains the code for conducting simulations for FN equation with general structural change.
- The R script ***Simulation-FN-Structural*** contains the code for conducting simulations for FN equation with systematic structural change.
- The R script ***Simulation-FN-Compare*** contains the code for comparing detection performances of different detection methods for general structural change with trajectory V(t) only.
- The R script ***Linear-Compare*** contains the code for conducting simulations for linear system.
- The R script ***Regression-Compare*** contains the code for conducting simulations for linear regression system.
- The R script ***Summarize-Simulation*** contains the code for summarizing simulation results.

3. Directory *Applications* contains the main functions and datasets for conducting analysis on the Italian COVID-19 data and the global temperature reconstruction dataset for the Holocene. Specifically, 
- The R script ***covid19Data_Generation*** cotains the code for generating the Italian COVID19 dataset from the original dataset through R package ``coronavirus''.
- The R script ***SIR*** contains the code for detecting structural changes in the Italian COVID-19 dataset.
- The R script ***Holocene Temperature*** contains the code for detecting structural changes in the global temperature reconstruction dataset for the Holocene.


### Workflows
Overall, please put the three folders in a same directory.
#### Simulation for comparing detection performance of different methods
1. Please create folders *Results_FN/* and *Results_Lin* for saving the simulation results of FN system and Linear (Linear Regression) system.
2. Run Rscripts ***Simulation-FN-General***, ***Simulation-FN-Structural***, ***Simulation-FN-Compare***, ***Linear-Compare*** and ***Regression-Compare*** to get the simulation results.
3. Run Rscript ***Summarize-Simulation*** to generate the simulation statistics and tables.

### Applications for COVID19 and Holocene Temperature
1. Run Rscript ***covid19Data_Generation*** to generate the Italian COVID19 dataset.
2. Run Rscript ***SIR*** to get the COVID19 application results and Figure 1.
3. Download the global temperature reconstruction dataset from [1] and put it in the folder Applications.
4. Run Rscript ***Holocene Temperature*** to get the Holocene Temperature application results and Figure 2.

## References
[1]. Marcott et al. A reconstruction of regional and global temperature for the past 11,300 years. Science, 339(6124):1198-1201, 2013.
