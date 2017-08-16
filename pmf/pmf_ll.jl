function pmf_ll(
            rows::Vector{Int},
            cols::Vector{Int},
            rates::Vector{Int},
            nrows::Int,
            ncols::Int,
            sr::Float64,
            su::Float64,
            sv::Float64,
            d::Int
            )::Tuple{Function,Function}

    n = length(rates)
    @assert length(rows)==length(cols)==n "inconsistent dimensions"

    function loglik(x) # neg log

        # x vector is a big column vector made of [u; v]
        # u is a big column vector made of [u1;u2;..;uN]
        # v is a big column vector made of [v1;v2;..;vP]

        u = x[1:(d*nrows)]
        v = x[(d*nrows+1):end]

        # accumulator for likelihood
        s = 0.0

        # computing everything as negative loglikelihood then returning neg
        for (k, rk) in enumerate(rates)
            # for a given rate r_ij, retrieve i,j then ui, uj
            i,j   = rows[k], cols[k]
            ui,vj = u[((i-1)*d+1):(i*d)], v[((j-1)*d+1):(j*d)]

            # then it's just a gaussian pdf and we add everything
            s += (rk - dot(ui,vj))^2/(2sr)
        end
        # add the likelihoods corresponding to spherical priors
        s += norm(u)^2/(2su) + norm(v)^2/(2sv)

        -s
    end
    function gradloglik(x)
        mu = 1:(d*nrows)             # mask for u in x
        mv = d*nrows + (1:(d*ncols)) # mask for v in x
        u = x[mu]
        v = x[mv]

        # accumulator for the gradient
        g = similar(x)

        # part associated with prior
        g[mu] += u/su
        g[mv] += v/sv

        # part associated with rates
        for (k, rk) in enumerate(rates)
            i,j   = rows[k], cols[k]
            mui   = ((i-1)*d+1):(i*d) # mask for ui in u
            mvj   = ((j-1)*d+1):(j*d) # mask for vj in v
            ui,vj = u[mui], v[mvj]
            # grad in ui
            g[mui] += (rk - dot(ui, vj)) * vj / (2sr)
            # grad in vj
            g[d*nrows+mvj] += (rk - dot(ui, vj)) * ui / (2sr)
        end
        
        g
    end
    (loglik,gradloglik)
end
