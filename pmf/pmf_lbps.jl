function pmf_lbps(data::Dict, params::Dict, verb::Bool=true)::Dict

    cprint("Starting LBPS...", verb) ; tb = time()

    rows  = data["ROWS"]
    cols  = data["COLS"]
    rates = data["RATES"]

    d  = params["LATENT_D"]
    sU = params["SIGMA_U"]
    sV = params["SIGMA_V"]
    sR = params["SIGMA_R"]
    λr = params["LAMBDAREF"]

    expname = params["EXPNAME"]
    maxN    = params["MAXNEVENTS"]
    maxT    = params["MAXT"]

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cprint("preparing the graph... ", verb) ; ta = time()

    # there may be discrepancy with lines missing etc
    nU = maximum(rows)
    nV = maximum(cols)

    # create N factors for the users
    mvgU             = MvDiagonalGaussian(zeros(d), sU*ones(d))
    gllU(x)          = gradloglik(mvgU, x)
    nexteventU(x, v) = nextevent_bps(mvgU, x, v)
    factorU(k)       = Factor(nexteventU, gllU, k)

    allfactors = [factorU(k) for k in 1:nU]

    # factors: create M factors for the movies,
    mvgV             = MvDiagonalGaussian(zeros(d), sV*ones(d))
    gllV(x)          = gradloglik(mvgV, x)
    nexteventV(x, v) = nextevent_bps(mvgV, x, v)
    factorV(k)       = Factor(nexteventV, gllV, nU+k)

    # stack all the factors
    allV = [factorV(k) for k in 1:nV]
    push!(allfactors, allV...)

    # Structure for factors: a vector formed of one part for U one for V
    maskU(x) = x[1:d]
    maskV(x) = x[d+1:end]

    # structure of the graph: each factor is connected to its own var
    structure = [[k] for k in 1:(nU+nV)]

    # now we need to add the factors corresponding to a specific rate
    # factors are stacked (1 -> nU), (nU+1 -> nU+nV), (nU+nV+1 -> nU+nV+nR)
    # a factor r_ij therefore corresponds to a connection [i to nU+j] in list

    for k in 1:length(rates)
        i, j, rij = rows[k], cols[k], rates[k]
        # the likelihood
        gij = PMFGaussian(rates[k], sR, d)
        # now we can build the factor
        fij = Factor( (x,w) -> nextevent_bps(gij, x, w),
                       x    -> gradloglik(gij, x),
                       nU+nV+k )
        # add the factor to the list
        push!(allfactors, fij)
        # and add the corresponding connection structure
        push!(structure, [i, nU+j])
    end

    # final step to define the graph
    fg = FactorGraph(structure, allfactors)

    cprintln("[done in $(round(time()-ta,1))s]", verb)

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cprint("initialising the simulation... ", verb) ; ta = time()

    nvars = nU + nV

    x0 = [sU*randn(d) for i in 1:nU]
    append!(x0, [sV*randn(d) for i in 1:nV])

    v0 = [randn(d) for i in 1:nvars]
    v0 = map(v->v/norm(v), v0)

    lsim = LocalSimulation(fg, x0, v0, maxT, maxN, λr)

    cprintln("[done in $(round(time()-ta,1))s]", verb)

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cprintln("starting the simulation... ", verb) ; ta = time()
    (all_evlist, details) = simulate(lsim)

    cprintln("... simulation finished ($(round(time()-ta,1))s)", verb)
    cprintln("... LBPS finished ($(round(time()-tb,1))s)", verb)

    # --------------------------
    # Returning relevant results

    results = Dict(
        "ALL_EVLIST"  => all_evlist,
        "SIM_DETAILS" => details,
        "RATES"       => rates       # the centred, scaled rates, for comp
        )
end
