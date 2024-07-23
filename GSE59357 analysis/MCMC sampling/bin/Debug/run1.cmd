model in Az_MCMC_ErbB_signaling_pathway.jags
data in jagsData.r
compile
initialize
update 10000
monitor a, thin(1)
update 20000
coda *, stem(CODA1)

