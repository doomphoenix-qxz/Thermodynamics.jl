# Thermodynamics.jl

## What is this package?

This is a Julia package that aims to become a thermodynamics engine and property database for as many chemicals as possible. I have grandiose hopes that it will one day become the de-facto standard for thermodynamic calculations all over the world. But the realities of opensource software development are what they are, and I can't make any guarantees.

## Why make this package when we have REFPROP/CoolProp/DIPPR/Aspen Properties/Insert my favorite thermo engine here?

Because I want a thermodynamic package that a) has data and properties for a large number of chemicals, b) is Open Source, and c) is written in Julia. I don't know of any other package satisfying these criteria.

## Why write in Julia instead of C/C++?

Because I like Julia. Because I don't have to write my own [nonlinear numerical methods](https://github.com/JuliaNLSolvers/NLsolve.jl). Because [adding physical units won't be hard](https://github.com/ajkeller34/Unitful.jl). And did I mention that I like Julia?

## Who are you?

I'm a chemical engineer. Physicists and physical chemists are arguably better thermodynamicists than chemical engineers, but a chemical engineer will do fine for equations of state and phase equilibria unless you really need statistical mechanics. (Contributions from anyone who knows statistical mechanics welcome!) Mechanical engineers do more of the cycle analysis and exergy analysis, so if any MechE's out there would like to contribute, that would be great too!

## What kind of thermodynamics functionality do you want to provide?

Just about anything you can think of. A database of thermodynamics and transport properties for as much as possible, several equations of state including reference equations for the fluids that have lots of data, phase equilibria, analysis of various cycles, perhaps even exergy.
