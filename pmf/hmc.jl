using ProgressMeter
# following http://www.inference.org.uk/mackay/itprnn/ps/387.412.pdf

function hmc(
            loglik::Function,
            gradloglik::Function,
            x0::Vector{Float64};
            steps::Int = 1000,
            burnin::Int = 100,
            lfsteps::Int = 2,
            stepsize::Float64 = 1e-3
            )

    @assert steps>burnin>0

    # Structure to store the results, each element = one sample
    res = Vector{Vector{Float64}}(steps-burnin)

    x   = x0
    dx  = length(x)
    ll  = loglik(x)
    gll = gradloglik(x)

    prog = Progress(steps, 1)

    for i in 1:steps

        p = randn(dx)          # initial momentum
        H = dot(p,p)/2 - ll    # hamiltonian

        xnew   = x
        gllnew = gll

        # Leapfrog
        for tau = 1:lfsteps
            p      += stepsize * gllnew / 2   # make half step in p
            xnew   += stepsize * p            # half step in x
            gllnew  = gradloglik(xnew)        # new gradient
            p      += stepsize * gllnew / 2   # half step in p
        end
        llnew = loglik(xnew)
        Hnew  = dot(p,p)/2 - llnew

        # Check new step
        ΔH = Hnew-H
        if (ΔH < 0) || (rand() < exp(-ΔH)) # accept
            g, x, ll  = gllnew, xnew, llnew
        end

        # Store as of burnin
        i > burnin ? res[i-burnin] = x  : nothing
        next!(prog)
    end
    res
end
