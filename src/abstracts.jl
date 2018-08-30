abstract type AbstractSubstance end
abstract type AbstractPhaseState end
abstract type AbstractMixture end
abstract type AbstractThermoBackend end
abstract type AbstractEquationOfState end
abstract type AbstractProperty end

function calculate(s::AbstractSubstance, st::AbstractPhaseState, b::AbstractThermoBackend, p::AbstractProperty) end
