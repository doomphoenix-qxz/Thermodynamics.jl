module Thermodynamics

include("abstracts.jl")
include("IdealGasBackend.jl")

export phasestate, spec_state!, mixture, thermocalc, thermocalc!, makethermodata, makethermoplot

end # module Thermodynamics
