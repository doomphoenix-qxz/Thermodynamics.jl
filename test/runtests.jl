using Thermodynamics
using Test

function test_peos(eos::Thermodynamics.parameterizedeos, δ, τ)
    print("Testing module parameterized_eos:\n")
    print("δ = "*string(δ))
    print("\nτ = "*string(τ))
    print("\nϕ₀ = "*string(ϕ₀(eos, δ, τ)))
    print("\nϕ₀_δ = "*string(ϕ₀_δ(eos, δ, τ)))
    print("\nϕ₀_δδ = "*string(ϕ₀_δδ(eos, δ, τ)))
    print("\nϕ₀_τ = "*string(ϕ₀_τ(eos, δ, τ)))
    print("\nϕ₀_ττ = "*string(ϕ₀_ττ(eos, δ, τ)))
    print("\nϕ₀_δτ = "*string(ϕ₀_δτ(eos, δ, τ)))
    print("\nϕr = "*string(ϕr(eos, δ, τ)))
    print("\nϕr_δ = "*string(ϕr_δ(eos, δ, τ)))
    print("\nϕr_δδ = "*string(ϕr_δδ(eos, δ, τ)))
    print("\nϕr_τ = "*string(ϕr_τ(eos, δ, τ)))
    print("\nϕr_ττ = "*string(ϕr_ττ(eos, δ, τ)))
    print("\nϕr_δτ = "*string(ϕr_δτ(eos, δ, τ)))
    print("\n\nNew test conditions:")
    ρ = 467.5483
    T = 305
    @printf "\nT = %f K, ρ = %f kg/m³" T ρ

    p_ = 7.5e6
    print("\np = "*string(p_))
    p_ = p(eos, ρ, T)
    print("\np = "*string(p_))
    print("\nSolving in reverse, with p = "string(p_)*" and T="*string(T))
    #print(ρ_from_p(eos, p_, 300))
    ρ_from_p(eos, p_, T)
end
