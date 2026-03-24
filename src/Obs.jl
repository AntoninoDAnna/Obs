module Obs
using ObsIO
using BitFlags
using ADerrors

AbstractCorr = ObsIO.AbstractCorr
include("enum.jl")
include("utilities.jl")
include("improvements.jl")
include("observables.jl")


end # module Obs
