struct parameterizedeos
  Tc::Float64
  ρc::Float64
  R::Float64

  #parameters for ideal gas portion
  n₀::Vector{Float64}
  γ₀::Vector{Float64}
end

h2o_n₀ = [-8.3204464837497, 6.6832105275932, 3.00632, 0.012436, 0.97315, 1.27950,
          0.96956, 0.24873]
h2o_γ₀ = [1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105]

function Σ(expr)
    return sum(eval(@. expr))
end

function ig(eos, δ, τ)
  end_ = Σ(eos.n₀[4:8]*log(1-exp(-eos.γ₀)*τ))
  return log(δ) + eos.n₀[1] + eos.n₀[2]*τ + eos.n₀[3]*log(τ) + end_
end

Tc = 647.096
ρc = 322
R = 0.46151805

eos = parameterizedeos(Tc,ρc,R,h2o_n₀,h2o_γ₀)
δ₁ = 838.025/ρc
τ₁ = Tc/500
print(ig(eos,δ₁,τ₁))
