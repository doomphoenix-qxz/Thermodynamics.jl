#parameterized_eos.jl
################################################################################
# Note: The following is based on the nomenclature in the paper
# "Revised Release on the IAPWS Formulation 1995 for the Thermodynamic
# Properties of Ordinary Water Substance for General and Scientific Use"
# by IAPWS, the International Association for the Properties of Water and Steam.
#
# It is more or less as follows:
#
# ϕ : Dimensionless Helmholtz energy of water as a fucntion of density and temperature
# δ : Dimensionless density of Water
# τ : Dimensionless temperature of water
# n, d, t, α, β, γ, ϵ, A, B, C, D : these are parameters used in fitting the
# equation of state to data. As far as we're concerned, they are constants with
# no real physical significance.
#
# Note that although the code was written by consulting the IAPWS paper, it is
# applicable to any
################################################################################


struct MultiparameterEOS{T_} where T_ <: AbstractFloat
    Tc::T_
    ρc::T_
    pc::T_
    ω::T_
    R::T_

    #parameters for ideal gas portion
    n₀::Vector{T_}
    γ₀::Vector{T_}

    # parameters for residual term 1
    n₁::Vector{T_}
    d₁::Vector{T_}
    t₁::Vector{T_}

    # parameters for residual term 2
    n₂::Vector{T_}
    d₂::Vector{T_}
    t₂::Vector{T_}
    c₂::Vector{T_}

    # parameters for residual term 3
    n₃::Vector{T_}
    d₃::Vector{T_}
    t₃::Vector{T_}
    α₃::Vector{T_}
    β₃::Vector{T_}
    γ₃::Vector{T_}
    ϵ₃::Vector{T_}

    # parameters for residual term 4
    n₄::Vector{T_}
    a₄::Vector{T_}
    b₄::Vector{T_}
    β₄::Vector{T_}
    A₄::Vector{T_}
    B₄::Vector{T_}
    C₄::Vector{T_}
    D₄::Vector{T_}

end

function ϕ₀(eos::MultiparameterEOS, δ, τ)
    first = log(δ) + eos.n₀[1] + eos.n₀[2]*τ + eos.n₀[3]*log(τ)
    end_ = sum(eos.n₀[4:8].*log.(1-exp.(-eos.γ₀.*τ)))
    return first + end_
end

function res_1(eos::MultiparameterEOS, δ, τ)
    return sum(eos.n₁ .* (δ .^ eos.d₁) .* (τ .^ eos.t₁))
end

function res_2(eos::MultiparameterEOS, δ, τ)
    return sum(eos.n₂ .* (δ .^ eos.d₂) .* (τ .^ eos.t₂) .* exp.(-δ .^ eos.c₂))
end

function res_3(eos::MultiparameterEOS, δ, τ)
    return sum(eos.n₃ .* (δ .^ eos.d₃) .* (τ .^ eos.t₃) .* exp.(-eos.α₃ .* (δ - eos.ϵ₃).^2 - eos.β₃ .* (τ - eos.γ₃).^2))
end

function θ(eos::MultiparameterEOS, δ, τ)
    return sum((1-τ) + eos.A₄.*((δ - 1).^2).^(1 ./(2 .* eos.β₄)))
end

function Δ(eos::MultiparameterEOS, δ, τ)
    return sum(θ(eos, δ, τ).^2 + eos.B₄ .* ((δ - 1).^2) .* eos.a₄)
end

function Ψ(eos::MultiparameterEOS, δ, τ)
    return sum(exp(-eos.C₄ .* (δ - 1)^2 - eos.D₄ .* (τ - 1)^2))
end

function res_4(eos::MultiparameterEOS, δ, τ)
    return sum(eos.n₄ .* Δ(eos,δ,τ).^eos.b₄ .* δ * Ψ(eos,δ,τ))
end

function ϕr(eos::MultiparameterEOS, δ, τ)
    return res_1(eos, δ, τ)+res_2(eos, δ, τ)+res_3(eos, δ, τ)+res_4(eos, δ, τ)
end

function ϕ₀_δ(eos::MultiparameterEOS, δ, τ)
    f(δ_) = ϕ₀(eos, δ_, τ)
    return ForwardDiff.derivative(f, δ)
end
function ϕ₀_δδ(eos::MultiparameterEOS, δ, τ)
    f(δ_) = ϕ₀_δ(eos, δ_, τ)
    return ForwardDiff.derivative(f, δ)
end
function ϕ₀_τ(eos::MultiparameterEOS, δ, τ)
    f(τ_) = ϕ₀(eos, δ, τ_)
    return ForwardDiff.derivative(f, τ)
end
function ϕ₀_ττ(eos::MultiparameterEOS, δ, τ)
    f(τ_) = ϕ₀_τ(eos, δ, τ_)
    return ForwardDiff.derivative(f, τ)
end
function ϕ₀_δτ(eos::MultiparameterEOS, δ, τ)
    return 0.0
end

function ϕr_δ(eos::MultiparameterEOS, δ, τ)
    f(δ_) = ϕr(eos, δ_, τ)
    return ForwardDiff.derivative(f, δ)
end
function ϕr_δδ(eos::MultiparameterEOS, δ, τ)
    f(δ_) = ϕr_δ(eos, δ_, τ)
    return ForwardDiff.derivative(f, δ)
end
function ϕr_τ(eos::MultiparameterEOS, δ, τ)
    f(τ_) = ϕr(eos, δ, τ_)
    return ForwardDiff.derivative(f, τ)
end
function ϕr_ττ(eos::MultiparameterEOS, δ, τ)
    f(τ_) = ϕr_τ(eos, δ, τ_)
    return ForwardDiff.derivative(f, τ)
end
function ϕr_δτ(eos::MultiparameterEOS, δ, τ)
    f(τ_) = ϕr_δ(eos, δ, τ_)
    return ForwardDiff.derivative(f, τ)
end



function p(eos::MultiparameterEOS, ρ, T)
    δ = ρ/eos.ρc
    τ = eos.Tc/T
    #@printf "δ = %f\n" δ
    #@printf "τ = %f\n" τ
    #@printf "derivative is %f\n" ϕr_δ(eos, δ, τ)
    p_nondim = 1 + δ*ϕr_δ(eos, δ, τ)
    return p_nondim * ρ*eos.R*T
end

function ρ_from_p(eos::MultiparameterEOS, p_, T)
    function solveit(ρ_guess)
        return abs(p_ - p(eos, ρ_guess, T))
    end
    a = optimize(solveit, 0.00001, 1500.0)
    print(a)
end

function u(eos::MultiparameterEOS, p_, T)
    ρ = ρ_from_p(eos, p_, T)
    δ = ρ/eos.ρc
    τ = eos.Tc/T
    u_nondim = τ * (ϕ₀_τ(eos, δ, τ) + ϕr_τ(eos, δ, τ))
    return u_nondim * T*eos.R
end

function s(eos::MultiparameterEOS, p_, T)
    ρ = ρ_from_p(eos, p_, T)
    δ = ρ/eos.ρc
    τ = eos.Tc/T
    s_nondim = τ * (ϕ₀_τ(eos, δ, τ) + ϕr_τ(eos, δ, τ)) - (ϕ₀(eos, δ, τ) + ϕr(eos, δ, τ))
    return s_nondim * eos.R
end

function H(eos::MultiparameterEOS, p_, T)
    ρ = ρ_from_p(eos, p_, T)
    δ = ρ/eos.ρc
    τ = eos.Tc/T
    H_nondim = 1 + τ * (ϕ₀_τ(eos, δ, τ) + ϕr_τ(eos, δ, τ)) + δ*ϕr_δ(eos, δ, τ)
    return H_nondim * T*eos.R
end

function Cp(eos::MultiparameterEOS, p_, T)
    ρ = ρ_from_p(eos, p_, T)
    δ = ρ/eos.ρc
    τ = eos.Tc/T
    part1 = -τ^2 * (ϕ₀_ττ(eos, δ, τ) + ϕr_ττ(eos, δ, τ))
    p2_top = (1 + δ*ϕr_δ(eos, δ, τ) - δ*τ*ϕr_δτ(eos, δ, τ))^2
    p2_bot = 1 + 2δ*ϕr_δ(eos, δ, τ) + δ^2*ϕr_δδ(eos, δ, τ)
    Cp_nondim = part1 + p2_top/p2_bot
    return Cp_nondim * eos.R
end


module IAPWS
using parameterized_eos
function gen_iapws()
    h2o_n₀ = [-8.3204464837497, 6.6832105275932, 3.00632, 0.012436, 0.97315, 1.27950,
              0.96956, 0.24873]
    h2o_γ₀ = [1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105]

    h2o_n₁ = [0.12533547935523e-1,
              0.78957634722828e1,
              -0.87803203303561e1,
              0.31802509345418,
              -0.26145533859358,
              -0.78199751687981e-2,
              0.88089493102134e-2]
    h2o_d₁ = [1, 1, 1, 2, 2, 3, 4]
    h2o_t₁ = [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1]

    h2o_n₂ = [-0.66856572307965,
              0.20433810950965,
              -0.66212605039687e-4,
              -0.19232721156002,
              -0.25709043003438,
              0.16074868486251,
              -0.40092828925807e-1,
              0.39343422603254e-6,
              -0.75941377088144e-5,
              0.56250979351888e-3,
              -0.15608652257135e-4,
              0.11537996422951e-8,
              0.36582165144204e-6,
              -0.13251180074668e-11,
              -0.62639586912454e-9,
              -0.10793600908932,
              0.17611491008752e-1,
              0.22132295167546,
              -0.40247669763528,
              0.58083399985759,
              0.49969146990806e-2,
              -0.31358700712549e-1,
              -0.74315929710341,
              0.47807329915480,
              0.20527940895948e-1,
              -0.13636435110343,
              0.14180634400617e-1,
              0.83326504880713e-2,
              -0.29052336009585e-1,
              0.38615085574206e-1,
              -0.20393486513704e-1,
              -0.16554050063734e-2,
              0.19955571979541e-2,
              0.15870308324157e-3,
              -0.16388568342530e-4,
              0.43613615723811e-1,
              0.34994005463765e-1,
              -0.76788197844621e-1,
              0.22446277332006e-1,
              -0.62689710414685e-4,
              -0.55711118565645e-9,
              -0.19905718354408,
              0.31777497330738,
              -0.11841182425981]
    h2o_c₂ = [1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 2,2,2,2,2,2,2,2, 2,2,2,2,2,2,2,2,2,
              2,2,2, 3,3,3,3, 4,6,6,6,6]
    h2o_d₂ = [1,1,1,2,2,3,4,4,5,7,9,10,11,13,15,1,2,2,2,3,4,4,4,5,6,6,7,9,9,9,9,9,
              10,10,12,3,4,4,5,14,3,6,6,6]
    h2o_t₂ = [4,6,12,1,5,4,2,13,9,3,4,11,4,13,1,7,1,9,10,10,3,7,10,10,6,10,10,1,2,3,
              4,8,6,9,8,16,22,23,23,10,50,44,46,50]

    h2o_d₃ = [3, 3, 3]
    h2o_t₃ = [0,1,4]
    h2o_n₃ = [-0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4]
    h2o_α₃ = [20,20,20]
    h2o_β₃ = [150,150,150]
    h2o_γ₃ = [1.21,1.21,1.25]
    h2o_ϵ₃ = [1,1,1]

    h2o_a₄ = [3.5,3.5]
    h2o_b₄ = [0.85,0.95]
    h2o_β₄ = [0.3,0.3]
    h2o_n₄ = [-0.14874640856724, 0.31806110878444]
    h2o_A₄ = [0.32,0.32]
    h2o_B₄ = [0.2,0.2]
    h2o_C₄ = [28,32]
    h2o_D₄ = [700,800]

    Tc = 647.096
    ρc = 322
    R = 0.46151805*1000

    return MultiparameterEOS(Tc,ρc,R,h2o_n₀,h2o_γ₀,h2o_n₁,h2o_d₁,h2o_t₁,
      h2o_n₂,h2o_d₂,h2o_t₂,h2o_c₂,h2o_n₃,h2o_d₃,h2o_t₃,h2o_α₃,h2o_β₃,h2o_γ₃,
      h2o_ϵ₃,h2o_n₄,h2o_a₄,h2o_b₄,h2o_β₄,h2o_A₄,h2o_B₄,h2o_C₄,h2o_D₄)
end
end # ends module IAPWS

module co2
import parameterized_eos
function gen_co2()
    co2_n₀= [ 8.37304456, -3.70454304,  2.5       ,  1.99427042,  0.62105248,
        0.41195293,  1.04028922,  0.08327678]

    co2_γ₀ = [3.15163,   6.1119 ,   6.77708,
            11.32384,  27.08792]

    co2_n₁ = [ 0.38856823,  2.93854759, -5.58671885, -0.767532  ,  0.31729004,
            0.54803316,  0.12279411]
    co2_d₁ = [ 1.,  1.,  1.,  1.,  2.,  2.,  3.]
    co2_t₁ = [ 0.  ,  0.75,  1.  ,  2.  ,  0.75,  2.  ,  0.75]
    co2_n₂ = [  2.16589615e+00,   1.58417351e+00,  -2.31327054e-01,
             5.81169164e-02,  -5.53691372e-01,   4.89466159e-01,
            -2.42757398e-02,   6.24947905e-02,  -1.21758602e-01,
            -3.70556853e-01,  -1.67758797e-02,  -1.19607366e-01,
            -4.56193625e-02,   3.56127893e-02,  -7.44277271e-03,
            -1.73957049e-03,  -2.18101213e-02,   2.43321666e-02,
            -3.47401334e-02,   1.43387158e-01,  -1.34919691e-01,
            -2.31512251e-02,   1.23631255e-02,   2.10583220e-03,
            -3.39585190e-04,   5.59936518e-03,  -3.03351181e-04]

    co2_d₂ = [  1.,   2.,   4.,   5.,   5.,   5.,   6.,   6.,   6.,   1.,   1.,
             4.,   4.,   4.,   7.,   8.,   2.,   3.,   3.,   5.,   5.,   6.,
             7.,   8.,  10.,   4.,   8.]

    co2_c₂ = [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,  2.,  2.,
            2.,  2.,  2.,  3.,  3.,  3.,  4.,  4.,  4.,  4.,  4.,  4.,  5.,  6.]

    co2_t₂ = [  1.5,   1.5,   2.5,   0. ,   1.5,   2. ,   0. ,   1. ,   2. ,
             3. ,   6. ,   3. ,   6. ,   8. ,   6. ,   0. ,   7. ,  12. ,
            16. ,  22. ,  24. ,  16. ,  24. ,   8. ,   2. ,  28. ,  14. ]

    co2_n₃ = [  -213.65488688,  26641.56914927, -24027.21220456,  -283.41603424,
              212.472844  ]
    co2_n₄ = [-0.66642277,  0.72608632,  0.05506867]
    co2_d₃ = [2.,2.,2.,3.,3.]
    co2_t₃ = [1.,0.,1.,3.,3.]
    co2_α₃ = [25.,25.,25.,15.,20.]
    co2_β₃ = [325.,300.,300.,275.,275.]
    co2_γ₃ = [1.16,1.19,1.19,1.25,1.22]
    co2_ϵ₃ = [1.,1.,1.,1.,1.]

    co2_a₄ = [3.5,3.5,3.0]
    co2_b₄ = [0.875,0.925,0.875]
    co2_β₄ = [0.3,0.3,0.3]
    co2_A₄ = [0.7,0.7,0.7]
    co2_B₄ = [0.3,0.3,1.0]
    co2_C₄ = [10.0,10.0,12.5]
    co2_D₄ = [275,275,275]

    Tc = 304.1282
    ρc = 467.6
    pc = 7.38e6
    R = 188.9241
    ω = 0.228
    return MultiparameterEOS(Tc,ρc,pc,ω,R,co2_n₀,co2_γ₀,co2_n₁,co2_d₁,co2_t₁,
      co2_n₂,co2_d₂,co2_t₂,co2_c₂,co2_n₃,co2_d₃,co2_t₃,co2_α₃,co2_β₃,co2_γ₃,
      co2_ϵ₃,co2_n₄,co2_a₄,co2_b₄,co2_β₄,co2_A₄,co2_B₄,co2_C₄,co2_D₄)
end
#eos = gen_co2()
#Tc = 304.1282
#ρc = 467.6
#δ₁ = 838.025/ρc
#τ₁ = Tc/500
#test_peos(eos, δ₁, τ₁)
end #ends module co2
