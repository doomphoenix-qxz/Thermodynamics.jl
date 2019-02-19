# module PengRobinson

using PolynomialRoots

struct PengRobinsonData{T} <: AbstractEquationOfState where T <: AbstractFloat
    CompoundName::String
    MW::T
    Pc::T
    Tc::T
    a::T
    b::T
    ω::T
    κ::T
end



function α(p_::PengRobinsonData, T)
    Tr = T/p_.Tc
    return (1+p_.κ*(1-√(Tr)))^2
end
function A(p_::PengRobinsonData, P, T)
    return p_.a * α(p_, T) * P / (R() ^2 * T^2)
end
function B(p_::PengRobinsonData, P, T)
    return p_.b*P/(T*R())
end
function c1(p_::PengRobinsonData, P, T)
    return (B(p_, P, T) - 1)
end
function c2(p_::PengRobinsonData, P, T)
    return A(p_, P, T) - 2B(p_, P,T) - 3B(p_, P,T)^2
end
function c3(p_::PengRobinson, p, T)
    return B(p_, P,T)^3 + B(p_, P,T)^2 - A(p_, P,T)*B(p_, P,T)
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

function volume_roots(p_::PengRobinsonData, P, T)
    c1_ = c1(p_, P, T)
    c2_ = c2(p_, P, T)
    c3_ = c3(p_, P, T)
    roots_ = roots([c3_, c2_, c1_, 1.0])
    return roots_
end

# end # ends module PengRobinson
