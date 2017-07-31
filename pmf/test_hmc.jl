include("hmc.jl")

S  = [1/3 1/2; 1/2 1];
iS = [12 -6; -6 4];
mu = [0.5, 0.5]

loglik(x)     = -dot(iS*(x-mu), (x-mu))/2
gradloglik(x) = iS*(mu-x)

x0 = [0.2, 0.2]

samples = hmc(loglik, gradloglik, x0; steps=10000, stepsize=1e-1)
