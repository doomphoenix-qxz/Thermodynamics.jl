# module PengRobinson

using PolynomialRoots

struct PengRobinsonData{T} <: AbstractCubicEOS where T <: AbstractFloat
    CompoundName::String
    MW::T
    Pc::T
    Tc::T
    a::T
    b::T
    ω::T
    κ::T
end



function α(eos::PengRobinsonData, T)
    Tr = T/eos.Tc
    return (1+eos.κ*(1-√(Tr)))^2
end
function A(eos::PengRobinsonData, P, T)
    return eos.a * α(eos, T) * P / (R() ^2 * T^2)
end
function B(eos::PengRobinsonData, P, T)
    return eos.b*P/(T*R())
end
function c1(eos::PengRobinsonData, P, T)
    return (B(eos, P, T) - 1)
end
function c2(eos::PengRobinsonData, P, T)
    return A(eos, P, T) - 2B(eos, P,T) - 3B(eos, P,T)^2
end
function c3(eos::PengRobinson, p, T)
    return B(eos, P,T)^3 + B(eos, P,T)^2 - A(eos, P,T)*B(eos, P,T)
end


function PengRobinsonData{T}(CompoundName::String, MW::T, Pc::T, Tc::T, ω::T) where T<: AbstractFloat
    a = 0.45724 * R() ^2 * Tc^2 /pc
    b = 0.0778*R() * Tc /pc
    κ = 0.37464 + 1.54226*ω - 0.26992*ω^2
    return PengRobinsonData(CompoundName, MW, Pc, Tc, a,b,ω,κ)
end

function get_the_real(vec::Vector{Complex{T}}) where {T<:Number}
    for i in vec
        if abs(imag(i)) < 2*eps(T)
            return real(i)
        end
    end
end

function volume_roots(eos::PengRobinsonData, P, T)
    c1_ = c1(eos, P, T)
    c2_ = c2(eos, P, T)
    c3_ = c3(eos, P, T)
    roots_ = roots([c3_, c2_, c1_, 1.0])
    return roots_
end

# end # ends module PengRobinson
