using Enums

abstract type AbstractSubstance end
abstract type AbstractPhaseState end
abstract type AbstractMixture end
abstract type AbstractThermoBackend end
abstract type AbstractEquationOfState<:AbstractThermoBackend end
abstract type AbstractMixingRules<:AbstractThermoBackend end
abstract type AbstractProperty end

@enum ThermoProp pressure=1 volume=2 temperature=3 entropy=4 internalenergy=5 enthalpy=6 helmholtzenergy=7  gibbsenergy=8 heatcapacity=9 vaporpressure=10 chemicalpotential=11 fugacity=12 fugacitycoefficient=13 activity=14 activitycoefficient=15
@enum MixingRules none=0 wongsandler=1

function thermocalc(st::AbstractPhaseState, eos::AbstractEquationOfState, p::AbstractProperty)
println("We're sorry, the property $show(p) isn't implemented for $show(st) , using the $show(eos) equation of state.") end
function thermocalc!(st::AbstractPhaseState, eos::AbstractEquationOfState, p::AbstractProperty) end

