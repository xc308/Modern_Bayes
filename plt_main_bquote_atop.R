#============
# Plots main 
#============

plot(1: 1e4, alpha_post, type = "l",
     xlab = "MCMC Iterations", 
     ylab = expression(paste("posterior of ", alpha)),
     #main = bquote(mu ~ ":" ~ .(mu2) ~ mm),
     main = bquote(atop("Posterior of "~alpha, 
                        atop("lag-1 acf:"~ .(round(acf_a_lag1, 2)),
                             "N_eff:"~ .(round(N_eff, 2))))),
     cex = 0.55)
