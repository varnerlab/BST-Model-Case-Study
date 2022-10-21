
# setup paths -
const _PROJECT_ROOT = pwd()
const _PATH_TO_SRC = joinpath(_PROJECT_ROOT, "src")
const _PATH_TO_MODELS = joinpath(_PROJECT_ROOT, "models")

# load my codes -
include(joinpath(_PATH_TO_SRC, "Utility.jl"))
include(joinpath(_PATH_TO_SRC, "Factory.jl"))
include(joinpath(_PATH_TO_SRC, "Kinetics.jl"))
include(joinpath(_PATH_TO_SRC, "Balances.jl"))
include(joinpath(_PATH_TO_SRC, "Simulation.jl"))
