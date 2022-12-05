# Rcode_PIM_composite_outcome

The folders and the R code help to reproducing Tables 1-3 in the article "Probabilistic index models for composite outcomes" by Zhiqiang Cao and Fan Li (under review) 

For questions or comments about the code, please contact Zhiqiang Cao <zcaoae@connect.ust.hk>. 
You will need to change the directory to use the example code in script.  This folder includes the following functions:

1. pim.co.R is to fits a probabilistic index model for composite outcomes, it can be used to fit standard PIM for composite outcomes, or covariates have heterogeneity

2. CreateScoreFun.R, delta.data.co.R, extract.times.R, pim.co.fit.R are functions used in pim.co.R

3. simulation_tables1_and_2.R is to reproduce simulation results in Tables 1-2 of paper

4. simulation_tables3_homoscedastic.R is to reproduce simulation results of normal linear regression with Homoscedastic in Tables 3 of paper

5. simulation_tables3_heteroscedastic.R is to reproduce simulation results of normal linear regression with Heteroscedastic in Tables 3 of paper

