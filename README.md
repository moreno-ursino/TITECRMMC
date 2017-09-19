In this folder, you can find the code to run the model and the simulation presented in “Dose-finding designs for cumulative toxicities using multiple constraints” paper in Biostatistics.

Here you are few guidelines to replicate the work or to use the TITECRMMC model.

scenarios_probability.R: this file helps to compute the probabilities of moderate toxicities and of DLTs for each scenario.

- genscen_titecrmmcLIK.R: you can find the code to simulate clinical trial data for each scenario presented in the paper and to run the estimation with the frequentist TITECRMMC. You can find the code also for the “probit” link function. To change scenario, change the transition matrices copying them from scenarios_probability.R

- sample2432.R: this file helps in running the estimation of MTD after 24 and 32 patients; you use the data from the previous clinical trial run for 40 patients. 

- bayes_empfix.R: to run the Bayesian TITECRMMC with empirical link function and fixed accrual.

- bayes_emppois.R: to run the Bayesian TITECRMMC with empirical link function and Poisson accrual.

- bayes_probfix.R: to run the Bayesian TITECRMMC with probit link function and fixed accrual.

- bayes_probpois.R: to run the Bayesian TITECRMMC with probit link function and Poisson accrual.

- titecrm_trial.R: to adapt the TITECRM to our generated datasets.

- bench_CRMMC.R: to obtain the results of the benchmark for the CRMMC.

The data used in the example of the paper correspond to the scenario 10 trial number 4.

The work was performed using R 3.2.1 with rstan  version 2.6.0. For any issue in updating the models for new version of rstan, please contact moreno.ursino@gmail.com.
