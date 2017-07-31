print("loading PDMP... "); ta = time()
using PDMP
println("[done in $(round(time()-ta,1))s]")

print("loading other packages... "); ta = time()
using JLD
include("pmf_lbps.jl")
println("[done in $(round(time()-ta,1))s]")

# Basic experiment format (this could be repeated in a for loop etc)

params = Dict(
    "EXPNAME"    => "A",   # name of the experiment
    "LATENT_D"   => 30,    # dimension of latent space
    "SIGMA_U"    => 20.0,
    "SIGMA_V"    => 20.0,
    "SIGMA_R"    =>  1.0,
    "LAMBDAREF"  =>  0.01, # refreshment rate
    "MAXNEVENTS" => 10,    # maximum number of events to generate
    "MAXT"       => Inf,   # maximum time
    "FROWS"      => "data/rows.csv",
    "FCOLS"      => "data/cols.csv",
    "FRATES"     => "data/rates.csv"
)

results = pmf_lbps(params)

save("exp_$(EXPNAME).jld",
        "RESULTS",     results,
        "SIM_PARAMS",  params)
