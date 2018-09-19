using Enums

abstract type AbstractSubstance end
abstract type AbstractPhaseState end
abstract type AbstractMixture end
abstract type AbstractThermoBackend end
abstract type AbstractEquationOfState end
abstract type AbstractProperty end

@enum ThermoProp pressure=1 volume=2 temperature=3 entropy=4 internalenergy=5 enthalpy=6 helmholtzenergy=7  gibbsenergy=8 heatcapacity=9 vaporpressure=10 chemicalpotential=11 fugacity=12 fugacitycoefficient=13 activity=14 activitycoefficient=15
# @enum MixingRules ...to implement later

function calculate(s::AbstractSubstance, st::AbstractPhaseState, b::AbstractThermoBackend, p::AbstractProperty) end
