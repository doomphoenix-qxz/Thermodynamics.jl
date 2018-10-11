module PengRobinson

using PolynomialRoots

struct PengRobinson{T<:Real}
    pc::T
    Tc::T
    a::T
    b::T
    ω::T
    κ::T
end



function α(p_::PengRobinson, T)
    Tr = T/p_.Tc
    return (1+p_.κ*(1-√(Tr)))^2
end
function A(p_::PengRobinson, p, T)
    return p_.a * α(p_, T) * p / (R() ^2 * T^2)
end
function B(p_::PengRobinson, p, T)
    return p_.b*p/(T*R())
end
function c1(p_::PengRobinson, p, T)
    return (B(p_, p, T) - 1)
end
function c2(p_::PengRobinson, p, T)
    return A(p_, p, T) - 2B(p_, p,T) - 3B(p_, p,T)^2
end
function c3(p_::PengRobinson, p, T)
    return B(p_, p,T)^3 + B(p_, p,T)^2 - A(p_, p,T)*B(p_, p,T)
end


function pengRobinson(pc, Tc, ω)
    a = 0.45724 * R() ^2 * Tc^2 /pc
    b = 0.0778*R() * Tc /pc
    κ = 0.37464 + 1.54226*ω - 0.26992*ω^2
    return pengRobinson(pc, Tc, a,b,ω,κ)
end

function get_the_real(vec::Vector{Complex{T}}) where {T<:Number}
    for i in vec
        if abs(imag(i)) < 2*eps(T)
            return real(i)
        end
    end
end

function volume_roots(p_::PengRobinson, p, T)
    c1_ = c1(p_, p, T)
    c2_ = c2(p_, p, T)
    c3_ = c3(p_, p, T)
    roots_ = roots([c3_, c2_, c1_, 1.0])
    return roots_
end

end # ends module PengRobinson
