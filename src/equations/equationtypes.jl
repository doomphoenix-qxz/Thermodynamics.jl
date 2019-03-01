struct IdealGasCpEquation1{T<:Number}
    vals::Array{T,1}
end

function thermocalc(eq::IdealGasCpEquation1, T)
    A,B,C,D,E = eq.vals
    return A+B*T+C*T^2+D*T^3+E/(T^2)
end
function thermocalc(eq::IdealGasCpEquation1, P, T)
    return thermocalc(eq, T)
end
