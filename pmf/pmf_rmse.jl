function pmf_rmse(rows, cols, rates, nrows, ncols, d, x)
    u  = x[1:(d*nrows)]
    v  = x[(d*nrows+1):end]

    mir = minimum(rates)
    mar = maximum(rates)

    s = 0.0
    for (k, rk) in enumerate(rates)
        i, j  = rows[k], cols[k]
        mui   = ((i-1)*d+1):(i*d)
        mvj   = ((j-1)*d+1):(j*d)
        ui,vj = u[mui], v[mvj]

        cand = dot(ui,vj)
        cand = (cand > mar) ? mar : cand
        cand = (cand < mir) ? mir : cand

        s += (rk - cand)^2
    end
    sqrt(s/length(rates))
end

function pmf_rmse2(rows, cols, rates, nrows, ncols, d, listx)
    meanpred = similar(rates)
    for x in listx
        pred = similar(rates)

        u    = x[1:(d*nrows)]
        v    = x[(d*nrows+1):end]

        mir = minimum(rates)
        mar = maximum(rates)

        s = 0.0
        for (k, rk) in enumerate(rates)
            i, j  = rows[k], cols[k]
            mui   = ((i-1)*d+1):(i*d)
            mvj   = ((j-1)*d+1):(j*d)
            ui,vj = u[mui], v[mvj]

            cand = dot(ui,vj)
            cand = (cand > mar) ? mar : cand
            cand = (cand < mir) ? mir : cand

            pred[k] = cand
            #s += (rk - cand)^2
        end
        meanpred += pred
    end
    meanpred /= length(listx)
    sqrt(mean((rates-meanpred).^2))
end
