#=============
# Chp 6 MCMC
#=============

# In low-dim, 
  # Important sampling and RJ sampling works pretty well

# but in high-dim
  # proposal worked in 2-D might not work in any dimension

# It's hard to capture high-dim space! 
  # turn to MCMC. 

# Markov Chain 
  # – where we go next only depends on our last state 
    # (the Markov property).
# Monte Carlo – just simulating data.

# Why MCMC?
  # (a) the region of high probability tends to be “connected”
    # get from one point to another without going through a low-probability region
  
  # (b) we tend to be interested in the expectations of functions 
    # that are relatively smooth and have lots of “symmetries”
    # That is, one only needs to evaluate them at a small number 
    # of representative points in order to get the general picture.
  

#---------------------------------  
# Advantages/Disadvantages of MCMC:
#---------------------------------  

## Advantanges:
  # applicable even when we can’t directly draw samples
  # works for complicated distributions in high-dimensional spaces,
    #even when we don’t know where the regions of high probability are
  # relatively easy to implement
  # fairly reliable


## Disadvantages:
  # slower than simple Monte Carlo or importance sampling 
    # (i.e.,requires more samples for the same level of accuracy

  # can be very difficult to assess accuracy and evaluate convergence,
    # even empirically


#---------------
# Ergodic Thrm
#---------------

# The ergodic theorem says: “if we start at a point xo and we
# keeping moving around in our high dimensional space, then we
# are guaranteed to eventually reach all points in that space with
# probability 1."


##-----------
# Metropolis Algorithm
##-----------

# Setup: Assume p.m.f.(mass) pi on Chi (countable)
          # i.e. measure sapce (Chi, sigma-algebra, pi).
        # f : Chi -> R

# Goal:
  # sample x from pi
  # apporximate E_pi[f(x)]

# supporting theory:
  # “If we take samples X = (X0, X1, . . . , ) 
    # then by the ergodic theorem,
      # they will eventually reach pi, 
        # which is known as the stationary distribution 
          # (the true pmf)."


# Approach : ergodic thrm
  # 1. If we run the Markov chain long enough, 
      # then the last state is approximately from pi.
  # 2. Under some regularity conditions,
      # 1/n * Sum_i f(x_i) --> E_pi[f(x)]














