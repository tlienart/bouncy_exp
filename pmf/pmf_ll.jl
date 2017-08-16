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
        u = x[1:(d*nrows)]
        v = x[(d*nrows+1):end]
        s = 0.0
        for (k, rk) in enumerate(rates)
            i,j   = rows[k], cols[k]
            ui,vj = u[((i-1)*d+1):(i*d)], v[((j-1)*d+1):(j*d)]
            s    += (rk - dot(ui,vj))^2/(2sr)
        end
        s += norm(u)^2/(2su) + norm(v)^2/(2sv)
        -s
    end
    function gradloglik(x)
        mu = 1:(d*nrows)
        mv = d*nrows + (1:(d*ncols))
        u = x[mu]
        v = x[mv]
        g = similar(x)
        g[mu] += u/su
        g[mv] += v/sv
        for (k, rk) in enumerate(rates)
            i,j   = rows[k], cols[k]
            mui   = ((i-1)*d+1):(i*d)
            mvj   = ((j-1)*d+1):(j*d)
            ui,vj = u[mui], v[mvj]
            # grad in ui
            g[mui] += (dot(ui, vj) - rk) * vj / (2sr)
            # grad in vj
            g[d*nrows+mvj] += (dot(ui, vj) - rk) * ui / (2sr)
        end
        g
    end
    (loglik,gradloglik)
end
