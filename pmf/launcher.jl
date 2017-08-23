print("loading PDMP... "); ta = time()
using PDMP
println("[done in $(round(time()-ta,1))s]")

print("loading other packages... "); ta = time()
using JLD
include("pmf_lbps.jl")
println("[done in $(round(time()-ta,1))s]")

cprint(s, b)   = b ? print(s)   : nothing
cprintln(s, b) = b ? println(s) : nothing

########################################################
# Idential preparation of the data for all experiments #
########################################################
verb = true

cprint("reading and preparing data... ", verb) ; ta = time()

rows  = vec(readdlm("data/rows.csv",  Int))
cols  = vec(readdlm("data/cols.csv",  Int))
rates = vec(readdlm("data/rates.csv", Float64))

# centre and scale the rates
range  = maximum(rates)-minimum(rates)
rates -= mean(rates)
rates /= range

data = Dict(
    "ROWS"  => rows,
    "COLS"  => cols,
    "RATES" => rates
)

base_sigma_r = 0.5 # Salakhutdinov & Mni
# https://pymc-devs.github.io/pymc3/notebooks/pmf-pymc.html
base_sigma_u = mean(var(rates, 2))
base_sigma_v = mean(var(rates, 1))

cprintln("[done in $(round(time()-ta,1))s]", verb)
########################################################

#################
# PMF with LBPS #
#################

lbpsparams = Dict(
    "EXPNAME"    => "A",   # name of the experiment
    "LATENT_D"   => 10,    # dimension of latent space
    "SIGMA_U"    => base_sigma_u,
    "SIGMA_V"    => base_sigma_v,
    "SIGMA_R"    => base_sigma_r,
    "LAMBDAREF"  => 0.01,  # refreshment rate
    "MAXNEVENTS" => 10,    # maximum number of events to generate
    "MAXT"       => Inf,   # maximum time
)
results = pmf_lbps(data, params)

save("exp_$(EXPNAME).jld",
        "RESULTS",     results,
        "SIM_PARAMS",  params)
