# Thermodynamics.jl

## What is this package?

This is a Julia package that aims to become a thermodynamics engine and property database for as many chemicals as possible. I have absurd, grandiose hopes that it will one day become the de-facto standard for thermodynamic calculations all over the world. But I make no guarantees that it will. In fact, I can't even guarantee that I'll have time to work on it past next week, although I very much hope that I will.

## Why make this package when we have REFPROP/CoolProp/DIPPR/Aspen Properties/Insert my favorite thermo engine here?

Because I want a thermodynamic package that a) has data and properties for a large number of chemicals, b) is Open Source, and c) is written in Julia. If your favorite thermo engine happens to satisfy all of these criteria, then by all means let me know, but I'm pretty sure I've seen most of the thermodynamics code written for Julia. There's not a whole lot. 

## Why write in Julia instead of C/C++?

Because I like Julia. Because I don't have to write my own [nonlinear numerical methods](https://github.com/JuliaNLSolvers/NLsolve.jl). Because [adding physical units won't be hard](https://github.com/ajkeller34/Unitful.jl). And did I mention that I like Julia?

## Who are you?

All you need to know is that I'm a chemical engineer. Physicists and physical chemists are arguably better thermodynamicists than chemical engineers, but a chemical engineer will do fine unless you really need statistical mechanics. (Contributions from anyone who knows statistical mechanics welcome!)

## What kind of thermodynamics functionality do you want to provide?

Just about anything you can think of. Thermodynamics and transport properties for anything I can, phase equilibria, analysis of various cycles, perhaps even exergy, though I'd need a mechanical engineer to help with that.
