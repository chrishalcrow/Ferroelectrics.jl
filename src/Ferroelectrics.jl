module Ferroelectrics

using Roots, ForwardDiff
using LinearAlgebra
using Makie, GLMakie, CairoMakie, LaTeXStrings
using DifferentialEquations, DiffEqCallbacks
using FileIO

export dpot

export export_arrow, SLf1, SLf2, SLden, stringfig, stringfig2, SLEXT


export imusingnotebook, imusingterminal, ConnectVacThroughAnotherVac, setgrid, InterpolateWalls
export xyzfsrt, srtfxyz, inbfsrt, srtfinb, getR2, getRvac, orderedtanhs, flow_s!

include("Parameters.jl")
export set_parameters!, get_As

include("Plots.jl")
export simple_plot, plot_field, plot_eng, energy

include("Derivatives.jl")
export dxDsing, d2xDsing

include("Properties.jl")
export baryon, PotEng

include("DiffRAK.jl")
export flow2DRAK!, flow2DRAKeng!, flow2DRAKno!, flow2DRAKnoN!

include("Diff.jl")
export flowRAK!, flowRAKf!


include("Modes.jl")
export getfreqs, getvecs, getH

include("2D.jl")
export flow2Del!, energy2Del, flow2D!, makeSkyrmeLine, getDDu!, getDDP!, dEdP2Dsf!, setgrid2D, energy2D, dEdP2Ds!, updateYbdy!
export energy2Del_noED, flow2D

include("String2D.jl")
export makeString, dist, string_length, string_positions, resample_string


include("Plots2D.jl")
export arrowplot, one_surface, two_surfaces, plot_surfaces, plot_strain_and_QPP, arrow_and_contour_plot

export arrow_and_contour_and_surface_plot, arrow_and_contour_plot2, arrow_and_contour_plot2D

include("Elastic.jl")
export getnbasis, getcs, get_C2m, GramSchmidtC, getALLcs, makeϵbf!, getc00, DOT_C, getQPP!, get_ϵbestL!, dotC, DOT_NO_C, get_ϵbestL_no_compat!
export check_compat_pt, get_string_energy


include("Figs/Fig5.jl")
include("Figs/Fig4.jl")
export getPts
include("Figs/Fig6.jl")
include("Figs/Fig8.jl")
include("Figs/Fig9.jl")
include("Figs/Fig10.jl")
include("Figs/Fig11.jl")
export fig5, getf4data, makefig5, getfig5data, getfig6data, makefig6,  getfig7data, makefig7, makefig8, makefig10
export getfig10data, getfig8data, flowRAKfig9!, makefig9, makefig11, makefig7, makefig5B





function InterpolateWalls(solT,solB,x,y)
    
    Lx = Lx = x[end-2] - x[3] 
    Nx = size(x)[1]
    Ny = size(y)[1]

    P0 = zeros(3,Nx,Ny)

    for a in 1:3, j in 1: Ny
    
        P0[a,1,j] = solT[a,j]
        P0[a,2,j] = solT[a,j]
   
        for i in 3:Nx-2
            P0[a,i,j] = (x[i]/Lx + 0.5)*solB[a,j] + (-x[i]/Lx + 0.5)*solT[a,j]
        end
    
        P0[a,Nx-1,j] = solB[a,j]
        P0[a,Nx,j] = solB[a,j]
        
    end

    return P0

end




function CtimesQ(C2,Q)

	qtemp = zeros(6,3,3)

	for α in 1:6, β in 1:6, b in 1:3, c in 1:3

		qtemp[α,b,c] += C2[α,β]*Q[β,b,c]

	end

	return qtemp

end




function one_surface(BD,x,y;scale=3.0)

    fig = Figure(resolution = (400, 400).*scale, fontsize = 10)
    ax = Axis3(fig[1, 1])
    sm = surface!(ax, x, y, BD)
    return fig

end

function two_surfaces(BD1,BD2,x,y;scale=3.0)

    fig = Figure(resolution = (800, 400).*scale, fontsize = 10)
    
    ax = Axis3(fig[1, 1])
    sm = surface!(ax, x, y, BD1)

    ax2 = Axis3(fig[1, 2])
    sm = surface!(ax2, x, y, BD2)

    return fig

end



function getRvac(ss, PV)

	st = ss - PV*dot(ss,PV)/dot(PV,PV);

	if sum(abs.(cross(st,PV))) == 0.0
		println("UH OH! PV and s are in the same direction!")
	end

	rt = cross(st,PV)

	return [ normalize!(rt) normalize!(st)  normalize!(PV) ]

end


function getR2(ss, rr)
	
	ss /= sqrt(sum(ss.*ss))
	rr /= sqrt(sum(rr.*rr))
	tt = cross(ss,rr)
	
	R2 = zeros(3,3)
	R2 .= hcat(ss, rr, tt)'


	
	return R2
	
end

function imusingnotebook()
	CairoMakie.activate!()
end

function imusingterminal()
	GLMakie.activate!()
end


function flow!(P,DP,N,dx,G2,A2,A4,R2,a111,a112,a123,V0,dt,n; engstops = 100000)
	
	ED = zeros(3,N)
	
	Em = 0.0
	E0 = 0.0
	Ep = 0.0
	counter = 1
	
	for t in 1:n
		
		DP = dEdP(P,N,dx,G2,A2,A4,R2,a111,a112,a123)
		
		for i in 3:N-2, a in 1:3
			P[a,i] -= dt*DP[a,i]
		end

		if t % engstops == 0
			
			println(energy(P, ED, dx, N, R2, G2, A2, A4, a111, a112, a123, V0) )
			#Em = E0
			#E0 = Ep
			#Ep = energy(P, ED, dx, N, R2, G2, A2, A4, a111, a112, a123, V0) 
			
			#println(Ep)
			
			#println( (Ep - E0)^2/( Ep*Em - E0^2 ) )
			
			#anton_criteria()
		end
		
	end
	
	
end


function flow_s!(P,DP,N,dx,G2,A2s,A4s,R2,a111,a112,a123,V0,dt,n; engstops = 100000)
	
	ED = zeros(3,N)
	
	Em = 0.0
	E0 = 0.0
	Ep = 0.0
	counter = 1

	dedp = zeros(3,N)
	
	for t in 1:n
		
		dEdPs!(dedp,P,N,dx,G2,A2s,A4s,R2,a111,a112,a123)
		
		for i in 3:N-2, a in 1:3
			P[a,i] -= dt*dedp[a,i]
		end

		if t % engstops == 0
			
			println(energy(P, ED, dx, N, R2, G2, A2s, A4s, a111, a112, a123, V0) )
			#Em = E0
			#E0 = Ep
			#Ep = energy(P, ED, dx, N, R2, G2, A2, A4, a111, a112, a123, V0) 
			
			#println(Ep)
			
			#println( (Ep - E0)^2/( Ep*Em - E0^2 ) )
			
			#anton_criteria()
		end
		
	end
	
	
end




function flow_tensor!(P,DP,N,dx,G2,A2,A4,R2,A6,V0,dt,n; engstops = 100000)
	
	ED = zeros(3,N)
	
	Em = 0.0
	E0 = 0.0
	Ep = 0.0
	counter = 1
	
	for t in 1:n
		
		DP = dEdP_tensor(P,N,dx,G2,A2,A4,R2,A6)
		
		for i in 3:N-2, a in 1:3
			P[a,i] -= dt*DP[a,i]
		end

		if t % engstops == 0
			
			println(energy_tensor(P, ED, dx, N, R2, G2, A2, A4, A6, V0) )

			#Em = E0
			#E0 = Ep
			#Ep = energy(P, ED, dx, N, R2, G2, A2, A4, a111, a112, a123, V0) 
			
			#println(Ep)
			
			#println( (Ep - E0)^2/( Ep*Em - E0^2 ) )
			
			#anton_criteria()
		end
		
	end
	
	
end




function orderedtanhs(P0,x,a,b,c,sep,R2)

	P = P0.*ones(3, size(x)[1] )
	X = zeros(3)

	if a == 1
		X[1] = -sep
	elseif a == 2
		X[2] = -sep
	else 
		X[3] = -sep
	end

	if b == 1
		X[1] = 0.0
	elseif b == 2
		X[2] = 0.0
	else 
		X[3] = 0.0
	end

	if c == 1
		X[1] = sep
	elseif c == 2
		X[2] = sep
	else 
		X[3] = sep
	end

	P[1,:] .= P0.*tanh.( x .- X[1])
	P[2,:] .= P0.*tanh.( x .- X[2])
	P[3,:] .= P0.*tanh.( x .- X[3])
	
	return srtfxyz(P,R2)

end


function makethreetanh(x,X)
	
	P = zeros(3, size(x)[1])
	
	P[1,:] .= tanh.( x .- X[1])
	P[2,:] .= tanh.( x .- X[2])
	P[3,:] .= tanh.( x .- X[3])
	
	return P
	
end

function ConnectVacThroughAnotherVac(x, PV, PV2)
	
	P = zeros(3,size(x)[1])
	
	for i in 1:size(x)[1], a in 1:3
		P[a,i] = PV[a]*tanh(x[i])
		P[a,i] += PV2[a]*sech(x[i])^2
	end
	
	return P
	
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



function setgrid(N,dx)
	
	println("grid has ", N, " points and goes from ", -0.5*dx*(N-1), " to ", 0.5*dx*(N-1))	
	return -0.5*dx*(N-1):dx:0.5*dx*(N-1)
	
end

function d1V6(Ppt, a111, a112, a123)
	@inbounds @fastmath Ppt[1]*( a111*6.0*Ppt[1]^4 + a112*(4.0*Ppt[1]^2*( Ppt[2]^2 + Ppt[3]^2 ) + 2.0*(Ppt[2]^4 +Ppt[3]^4) ) + 2.0*a123*Ppt[2]^2*Ppt[3]^2 )
end

function d2V6(Ppt, a111, a112, a123)
	@inbounds @fastmath Ppt[2]*(a111*6.0*Ppt[2]^4 + a112*(4.0*Ppt[2]^2*( Ppt[3]^2 + Ppt[1]^2 ) + 2.0*(Ppt[3]^4 +Ppt[1]^4) ) + 2.0*a123*Ppt[3]^2*Ppt[1]^2)
end

function d3V6(Ppt, a111, a112, a123)
	@inbounds @fastmath Ppt[3]*(a111*6.0*Ppt[3]^4 + a112*(4.0*Ppt[3]^2*( Ppt[1]^2 + Ppt[2]^2 ) + 2.0*(Ppt[1]^4 +Ppt[2]^4) ) + 2.0*a123*Ppt[1]^2*Ppt[2]^2)
end

function d11V6(Ppt, a111, a112, a123)
	a111*30.0*Ppt[1]^4 + a112*(12.0*Ppt[1]^2*( Ppt[2]^2 + Ppt[3]^2 ) + 2.0*(Ppt[2]^4 +Ppt[3]^4) ) + 2.0*a123*Ppt[2]^2*Ppt[3]^2
end

function d22V6(Ppt, a111, a112, a123)
	a111*30.0*Ppt[2]^4 + a112*(12.0*Ppt[2]^2*( Ppt[3]^2 + Ppt[1]^2 ) + 2.0*(Ppt[3]^4 +Ppt[1]^4) ) + 2.0*a123*Ppt[3]^2*Ppt[1]^2
end

function d33V6(Ppt, a111, a112, a123)
	a111*30.0*Ppt[3]^4 + a112*(12.0*Ppt[3]^2*( Ppt[1]^2 + Ppt[2]^2 ) + 2.0*(Ppt[1]^4 +Ppt[2]^4) ) + 2.0*a123*Ppt[1]^2*Ppt[2]^2
end

function d12V6(Ppt, a111, a112, a123)
	 a112*(4.0*Ppt[1]^3*( 2.0*Ppt[2] ) + 2.0*Ppt[1]*(4.0*Ppt[2]^3 ) ) + 4.0*a123*Ppt[1]*Ppt[2]*Ppt[3]^2
end

function d21V6(Ppt, a111, a112, a123)
	 a112*(4.0*Ppt[1]^3*( 2.0*Ppt[2] ) + 2.0*Ppt[1]*(4.0*Ppt[2]^3 ) ) + 4.0*a123*Ppt[1]*Ppt[2]*Ppt[3]^2
end

function d13V6(Ppt, a111, a112, a123)
	 a112*(4.0*Ppt[1]^3*( 2.0*Ppt[3] ) + 2.0*Ppt[1]*(4.0*Ppt[3]^3 ) ) + 4.0*a123*Ppt[1]*Ppt[2]^2*Ppt[3]
end

function d31V6(Ppt, a111, a112, a123)
	 a112*(4.0*Ppt[1]^3*( 2.0*Ppt[3] ) + 2.0*Ppt[1]*(4.0*Ppt[3]^3 ) ) + 4.0*a123*Ppt[1]*Ppt[2]^2*Ppt[3]
end

function d23V6(Ppt, a111, a112, a123)
	 a112*(4.0*Ppt[2]^3*( 2.0*Ppt[3] ) + 2.0*Ppt[2]*(4.0*Ppt[3]^3 ) ) + 4.0*a123*Ppt[2]*Ppt[3]*Ppt[1]^2
end

function d32V6(Ppt, a111, a112, a123)
	 a112*(4.0*Ppt[2]^3*( 2.0*Ppt[3] ) + 2.0*Ppt[2]*(4.0*Ppt[3]^3 ) ) + 4.0*a123*Ppt[2]*Ppt[3]*Ppt[1]^2
end


function dV6(Ppt, a111, a112, a123)
	[ d1V6(Ppt, a111, a112, a123), d2V6(Ppt, a111, a112, a123), d3V6(Ppt, a111, a112, a123) ]
end




function dEdP(P,N,dx,G2,A2,A4,R2, a111, a112, a123)

	dedp = zeros(3,N)
	
	DDPpt = zeros(3)
	Ppt = zeros(3)
	
	
	for i in 3:N-2
		
		Ppt = getP(P,i)
		getDDP!(DDPpt,P,i,dx)
		
		#println(DDPpt)
	
		for a in 1:3, b in 1:3
		
			dedp[a,i] += -G2[a,b]*DDPpt[b]
			
			dedp[a,i] += (A2[a,b] + A2[b,a])*Ppt[b]
			
			for c in 1:3, d in 1:3
				dedp[a,i] += (A4[a,b,c,d] + A4[b,a,c,d] + A4[b,c,a,d] + A4[b,c,d,a])*Ppt[b]*Ppt[c]*Ppt[d]
			end
			
		end
			
		dedp[:,i] .+= R2*dV6(R2'*Ppt,a111, a112, a123)
		
	end
	
	return dedp
	
end


function dEdPs(P,N,dx,G2,A2s,A4s,R2, a111, a112, a123)

	dedp = zeros(3,N)
	
	DDPpt = zeros(3)
	Ppt = zeros(3)

	dD6 = zeros(3)
	R2Pt = zeros(3)
	
	
	for i in 3:N-2


		
		#Ppt = getP(P,i)

		Ppt[1] = P[1,i] 
		Ppt[2] = P[2,i] 
		Ppt[3] = P[3,i] 

		getDDP!(DDPpt,P,i,dx)
		
		#println(DDPpt)
	
		for a in 1:3

			R2Pt[a] = 0.0 


			for b in 1:3

				R2Pt[a] += R2[b,a]*Ppt[b]
			
				dedp[a,i] += -G2[a,b]*DDPpt[b]
				
				dedp[a,i] += 2.0*A2s[a,b]*Ppt[b]
				
				for c in 1:3, d in 1:3
					dedp[a,i] += 4.0*A4s[a,b,c,d]*Ppt[b]*Ppt[c]*Ppt[d]
				end

			end
			
		end

		dD6[1] = d1V6(Ppt, a111, a112, a123)
		dD6[2] = d2V6(Ppt, a111, a112, a123)
		dD6[3] = d3V6(Ppt, a111, a112, a123)

		for a in 1:3, b in 1:3

			dedp[a,i] += R2[a,b]*dD6[b]
			
		end
		
	end
	
	return dedp
	
end


function dEdPs!(dedp,P,N,dx,G2,A2s,A4s,R2, a111, a112, a123)

	
	DDPpt = zeros(3)
	Ppt = zeros(3)

	dD6 = zeros(3)
	R2Pt = zeros(3)
	
	
	for i in 3:N-2
		
		Ppt[1] = P[1,i] 
		Ppt[2] = P[2,i] 
		Ppt[3] = P[3,i] 

		getDDP!(DDPpt,P,i,dx)
		
		#println(DDPpt)
	
		for a in 1:3

			dedp[a,i] = 0.0

			R2Pt[a] = 0.0 


			for b in 1:3

				R2Pt[a] += R2[b,a]*Ppt[b]
			
				dedp[a,i] += -G2[a,b]*DDPpt[b]
				
				dedp[a,i] += A2s[a,b]*Ppt[b]
				
				for c in 1:3, d in 1:3
					dedp[a,i] += A4s[a,b,c,d]*Ppt[b]*Ppt[c]*Ppt[d]
				end

			end
			
		end

		#dD6 = dV6(R2Pt,a111,a112,a123)
		dD6[1] = d1V6(Ppt, a111, a112, a123)
		dD6[2] = d2V6(Ppt, a111, a112, a123)
		dD6[3] = d3V6(Ppt, a111, a112, a123)

		for a in 1:3, b in 1:3

			dedp[a,i] += R2[a,b]*dD6[b]
			
		end
		
	end
	
	
end
			
		


function dEdP_tensor(P,N,dx,G2,A2,A4,R2, A6)

	dedp = zeros(3,N)
	
	DDPpt = zeros(3)
	Ppt = zeros(3)
	
	
	for i in 3:N-2
		
		Ppt = getP(P,i)
		getDDP!(DDPpt,P,i,dx)
		
		#println(DDPpt)
	
		for a in 1:3, b in 1:3
		
			dedp[a,i] += -G2[a,b]*DDPpt[b]
			
			dedp[a,i] += (A2[a,b] + A2[b,a])*Ppt[b]
			
			for c in 1:3, d in 1:3
				dedp[a,i] += (A4[a,b,c,d] + A4[b,a,c,d] + A4[b,c,a,d] + A4[b,c,d,a])*Ppt[b]*Ppt[c]*Ppt[d]
			

				for e in 1:3, f in 1:3
					dedp[a,i] += (A6[a,b,c,d,e,f] + A6[b,a,c,d,e,f] + A6[b,c,a,d,e,f] + A6[b,c,d,a,e,f] + A6[b,c,d,e,a,f] + A6[b,c,d,e,f,a] )*Ppt[b]*Ppt[c]*Ppt[d]*Ppt[e]*Ppt[f]
				end

			end
			
		end
			
		#dedp[:,i] .+= R2*dV6(R2'*Ppt,a111, a112, a123)
		
	end
	
	return dedp
	
end
			








#=
ss = [ 1.0, -1.0, 0.0 ]
rr = [ 1.0,1.0,0.0 ]
tt = cross(ss,rr)


tt /= sqrt(sum(tt.*tt))
ss /= sqrt(sum(ss.*ss))
rr /= sqrt(sum(rr.*rr))

	R2 = hcat(ss,rr,tt)
	=#
	
end