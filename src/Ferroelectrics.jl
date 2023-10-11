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
	is_electrostatic::Bool
end

mutable struct Grid
	lp::Int64
    ls::Float64
    x::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
end

mutable struct PolarisationField1D
    field::Matrix{Float64}
	grid::Grid
    parameters::Parameter
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


PolarisationField1D(lp,ls,R2,material;is_electrostatic=false) = PolarisationField1D( zeros(3,lp), setgrid(lp, ls), set_parameters_lithium!(R2,material, is_electrostatic))

function setgrid(N,dx)
	
	println("grid has ", N, " points and goes from ", -0.5*dx*(N-1), " to ", 0.5*dx*(N-1))	
	return Grid(N, dx, -0.5*dx*(N-1):dx:0.5*dx*(N-1) )
	
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