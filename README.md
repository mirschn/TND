# TND
Test-negative design simulation study for vaccine effectiveness

Simulation study for the paper "Estimands and estimation of COVID-19 vaccine effectiveness under the test-negative design: connections to causal inference" (available soon).

## Scenario 1
No effect modification, no interference

## Scenario 2
Effect modification: the covariate C interacts with V (vaccination) in the model for infection with SARS-CoV-2.
No interference.

## Scenario 3
Partial interference: patient recruited from 10 blocks. The % infected in the block acts as an effect modifier for the vaccine in the model for infection with SARS-CoV-2.
No other effect modification.

There are 5 files:

## datagenv2_pub
Data-generation function for the first and second scenarios, without interference.

## datagenv2_interference_pub
Data-generation function for the third scenario, with interference.

## truevals_pub
Demonstrates how to use the previous two functions to obtain the true values of the various parameters of interest described in the paper (i.e. the conditional risk ratio, the marginal risk ratio, and the direct effects under interference)

## simstudy_all_Aug2021v2_pub ; simstudy_int_Aug2021_pub
Files that run the simulation study. The first file has the scenarios without interference; the second file has the scenario with partial interference. The methods are logistic regression and inverse probability of treatment weighting where the propensity score is fit using 1) all subjects and 2) only control subjects.

## simstudy_int_Nov2021_byblock
Runs analyses on a single simulated data-set in order to demonstrate stratification by block under interference.
