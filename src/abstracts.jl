using Enums

abstract type AbstractSubstance end
abstract type AbstractPhaseState end
abstract type AbstractMixture end
abstract type AbstractThermoBackend end
abstract type AbstractEquationOfState end
abstract type AbstractProperty end

@enum ThermoProp pressure=1 volume=2 temperature=3 entropy=4 internalenergy=5 enthalpy=6 helmholtzenergy=7  gibbsenergy=8 

function calculate(s::AbstractSubstance, st::AbstractPhaseState, b::AbstractThermoBackend, p::AbstractProperty) end
