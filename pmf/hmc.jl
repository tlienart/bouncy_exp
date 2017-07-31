
# following http://www.inference.org.uk/mackay/itprnn/ps/387.412.pdf

function hmc(loglik::Function, gradloglik::Function,
            x0::Vector{Float64}; steps=1000, burnin=100,
            lfsteps=2, stepsize=1e-3)

    @assert steps>burnin>0

    res = Vector{Vector{Float64}}(steps-burnin)

    x   = x0
    ll  = loglik(x)
    gll = gradloglik(x)


    for i in 1:steps

        p = randn(length(x0))   # initial momentum
        H = dot(p,p)/2 - ll     # hamiltonian

        xnew   = x
        gllnew = gll

        # Leapfrog
        for tau = 1:lfsteps
            p      += stepsize * gllnew / 2 # make half step in p
            xnew   += stepsize * p          # half step in x
            gllnew  = gradloglik(xnew)      # new gradient
            p      += stepsize * gllnew / 2 # half step in p
        end

        llnew = loglik(xnew)
        Hnew  = dot(p,p)/2 - llnew

        ΔH = Hnew-H

        accept = ΔH < 0 || (rand() < exp(-ΔH))
        if accept
            g  = gllnew
            x  = xnew
            ll = llnew
        end

        i > burnin ? res[i-burnin] = x  : nothing
    end
    res
end
