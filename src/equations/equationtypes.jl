struct CpEquation1{T<:Number}
    vals::Array{T,1}
end

function thermocalc(eq::CpEquation1, Temp)
    A,B,C,D,E = vals
    return A+B*Temp+C*Temp^2+D*Temp^3+E/(Temp^2)
end
