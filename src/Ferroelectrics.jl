module Ferroelectrics

using Roots
using LinearAlgebra
using Makie, CairoMakie, LaTeXStrings
using FileIO
using StaticArrays
using FFTW

mutable struct Parameter
	PV::Vector{Float64}
	V0::Float64
	A2::Matrix{Float64}
	A4::Array{Float64, 4}
	G2::Matrix{Float64}
	A2t::Matrix{Float64}
	A4t::Array{Float64, 4}
	A2s::Matrix{Float64}
	A4s::Array{Float64, 4}
	rescaling::Float64
	length_scale::Float64
end


mutable struct PolarisationField1D
    field::Matrix{Float64}
    lp::Int64
    ls::Float64
    x::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    parameters::Parameter
	is_electrostatic::Bool
	b::Float64
end

export PolarisationField1D, xyzfsrt, srtfxyz, inbfsrt, srtfinb

include("Parameters.jl")
export set_parameters_lithium!, set_parameters_barium!

include("Plots.jl")
export plot_field

include("Derivatives.jl")

include("Properties.jl")
export changeAs, energy

include("Diff.jl")
export gradient_flow!

PolarisationField1D(lp,ls,R2,material,b) = PolarisationField1D( zeros(3,lp), lp, ls, -ls*(lp - 1)/2.0 : ls : ls*(lp - 1)/2.0, set_parameters_lithium!(R2,material) , false, b)


function setgrid(N,dx)
	
	println("grid has ", N, " points and goes from ", -0.5*dx*(N-1), " to ", 0.5*dx*(N-1))	
	return -0.5*dx*(N-1):dx:0.5*dx*(N-1)
	
end

function xyzfsrt(P,R2)
	return R2'*P
end

function srtfxyz(P,R2)
	return R2*P
end

function inbfsrt(P,R2)
	return R2*P
end

function srtfinb(P,R2)
	return R2*P
end


	
end