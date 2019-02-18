#module IdealGasBackend

# One important thing to do here is to make sure that we can handle cases
# with missing data. For Ideal Gas, the possible missing pieces of data are
# Heat Capacity (Cp) data and transport properties.
# The constructors for state structs that depend on heat capacity data (such
# as IdealGasPotentialsState) are responsible to throw an error if heat capacity
# data is missing.

struct IdealGasData{T_} where T <: AbstractFloat
    CompoundName::String
    MW::T_
    CpCoefficients::Union{Vector{T_}, Missing}
    CpEquationNumber::Union{Int, Missing}
end

struct IdealGasMinimalState{T_} where T_ <: AbstractFloat
    P::T_
    V::T_
    ρ::T_
    T::T_
    Cp::Union{T_, Missing}
end

struct IdealGasPotentialsState{T_} where T_ <: AbstractFloat
    S::T_
    U::T_
    H::T_
    A::T_
    G::T_
end

struct IdealGasChemicalState{T_} where T_ <: AbstractFloat
    xᵢ::T_
    f::T_
    ϕ::T_
end

struct IdealGasTransportState{T_} where T_ <: AbstractFloat
    λ::Union{T_, Missing}
    μ::Union{T_, Missing}
    ν::Union{T_, Missing}
    α::Union{T_, Missing}
    Pr::Union{T_, Missing}
end

struct IdealGasFullState{T_} where T_ <: AbstractFloat
    data::IdealGasData{T_}
    minstate::IdealGasMinimalState{T_}
    potentials::IdealGasPotentialsState{T_}
    chemstate::IdealGasChemicalState{T_}
    transtate::IdealGasTransportState{T_}
end





#end # module  IdealGasBackend
