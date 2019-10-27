# Blocking for toy model in Matlab

Matlab implementation of the toy model (Andrieu, 2010)
with `sigma_x=10.0` and `sigma_y=1.0`.

Running the script "simulateBlock" will produce results for a run of M=10000 MCMC steps with N=50 particles for particle Gibbs with ancestor sampling, marginalized particle Gibbs with ancestor sampling and marginalized particle Gibbs with ancestor sampling + blocking. All parameters are set to the same values as in section 4.1 in the paper.
