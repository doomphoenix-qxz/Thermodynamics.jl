module Thermodynamics

using ThermodynamicsBase
using ForwardDiff: derivative
using Reexport
using NLSolve: nlsolve
using CSV

include("backends/idealgas.jl")
include("backends/multiparameterhelmholtz.jl")
include("backends/pengrobinson.jl")

export phasestate, spec_state!, mixture, thermocalc, thermocalc!, makethermodata, makethermoplot

end # module Thermodynamics
