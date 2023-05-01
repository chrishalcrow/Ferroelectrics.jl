






function dudf2Del!(dP,P,p,t)

	pars,elasticstuff, ϵ = p
	Nx,Ny,dx,dy,G4,A2s,A4s, a111, a112, a123, C2,q3,Q3,R2 = pars

	Ppt = zeros(3)
	ept = zeros(6)
	DDPpt = zeros(3,2,2)

	dD6 = zeros(3)
	R2Pt = zeros(3)
	
	
	@inbounds for i in 3:Nx-2, j in 3:Ny-2
		
		for a in 1:3 
			Ppt[a] = P[a,i,j]
		end 
		for a in 1:6
			ept[a] = ϵ[a,i,j]
		end

		getDDP!(DDPpt,P,i,j,dx,dy)
			
		for a in 1:3

			dP[a,i,j] = 0.0

			R2Pt[a] = 0.0 


			for b in 1:3

				for α in 1:6
					dP[a,i,j] -= (q3[α,a,b]+q3[α,b,a])*ept[α]*Ppt[b]
				end

				for c in 1:2, d in 1:2
					dP[a,i,j] += -0.5*(G4[c,a,d,b]+G4[c,b,d,a])*DDPpt[b,c,d]
				end

				R2Pt[a] += R2[b,a]*Ppt[b]
			
				dP[a,i,j] += 2.0*A2s[a,b]*Ppt[b]
				
				for c in 1:3, d in 1:3
					dP[a,i,j] += 4.0*A4s[a,b,c,d]*Ppt[b]*Ppt[c]*Ppt[d]
				end

			end
			
		end

		dD6[1] = d1V6(R2Pt, a111, a112, a123)
		dD6[2] = d2V6(R2Pt, a111, a112, a123)
		dD6[3] = d3V6(R2Pt, a111, a112, a123)

		for a in 1:3, b in 1:3

			dP[a,i,j] += R2[a,b]*dD6[b]
			
		end

		for a in 1:3
			dP[a,i,j] *= -1.0
		end
		
	end



	
	
end







function flow2DRAK!(P0,ϵ0,tf,pars,elasticstuff,QPP,ϵbasis)
    
    EList  = zeros(1000)

	counter = 1

	p = [pars,elasticstuff, ϵ0, QPP, ϵbasis, counter, EList]


    tspan = (0.0,tf)
	
	cb = DiscreteCallback(condition,affect!, save_positions=(false,false))

    prob = ODEProblem(dudf2Del!,P0,tspan,p)
	sol = solve(prob, Tsit5(), save_everystep=false, callback=cb)


    return sol[end], EList
    
end




function condition(u, t, integrator)
	true 
end 

function affect!(integrator)




	#Nx = integrator.p[1][1]
	#Ny = integrator.p[1][2]
	#print(integrator.t, ",   " )

	#integrator.p[6] += 1.0

	cc00,CL,ncsb,Lx,Ly,pmax,qmax,C2m = integrator.p[2]
	Nx,Ny,dx,dy,G4,A2s,A4s, a111, a112, a123, C2,q3,Q3,R2 = integrator.p[1]
	QPP = integrator.p[4]
	ϵbasis = integrator.p[5]
	#ϵ0 = integrator.p[3]
	#integrator.p[6] += 1.0



#counter += 1.0

	getQPP!(QPP, Q3, integrator.u, Nx, Ny)

	updateYbdy!(integrator.u, Nx, Ny)

   	#getQPP!(QPP,p[1][13],integrator.u,Nx,Ny)

   		#if compat == true
			get_ϵbestL!(integrator.p[3],ϵbasis,QPP,cc00,CL,ncsb,Nx,Ny,dx,dy,Lx,Ly,pmax,qmax,C2m)
		#else
			#get_ϵbestL_no_compat!(integrator.p[3],QPP,Nx,Ny)
		#end

	integrator.p[7][integrator.p[6]] = energy2Del_noED(integrator.u, integrator.p[3], G4, A2s, A4s, C2, q3, a111,a112,a123,0.0,R2,Nx,Ny,dx,dy)  
	
	integrator.p[6] += 1
	#println(integrator.p[6])

	#println( energy2Del_noED(integrator.u, integrator.p[3], G4, A2s, A4s, C2, q3, a111,a112,a123,0.0,R2,Nx,Ny,dx,dy) ) 
	


end




function conditionno(u, t, integrator)
	true 
end 

function affectno!(integrator)

	Nx,Ny,dx,dy,G4,A2s,A4s, a111, a112, a123, C2,q3,Q3,R2, Eext = integrator.p[1]

	updateYbdy!(integrator.u, integrator.p[1][1], integrator.p[1][2])

	if integrator.p[2][1] % 10 == 9
		integrator.p[3][ (integrator.p[2][1]+10) ÷ 10 ] =  energy2D_noED(integrator.u, G4, A2s, A4s,a111,a112,a123,0.0,R2,Nx,Ny,dx,dy) 
	end

	integrator.p[2][1] += 1

end

function flow2DRAKno!(P0,tf,pars)
    

	counter = 1
	EList = zeros(1000)

	p = (pars,[counter],EList)

    tspan = (0.0,tf)
	
	cb = DiscreteCallback(conditionno,affectno!, save_positions=(false,false))

    prob = ODEProblem(dudf2Dno!,P0,tspan,p)
	sol = solve(prob, Tsit5(), save_everystep=false, callback=cb)


    P0 .= sol[end]

    return EList[1:p[2][1] ÷ 10]
    
end





function dudf2Dno!(dP,P,p,t)

	
	Nx,Ny,dx,dy,G4,A2s,A4s, a111, a112, a123, C2,q3,Q3,R2, Eext = p[1]

	Ppt = zeros(3)
	DDPpt = zeros(3,2,2)

	dD6 = zeros(3)
	R2Pt = zeros(3)
	
	
	@inbounds for i in 3:Nx-2, j in 3:Ny-2
		
		Ppt[1] = P[1,i,j] 
		Ppt[2] = P[2,i,j] 
		Ppt[3] = P[3,i,j] 

		getDDP!(DDPpt,P,i,j,dx,dy)
			
		for a in 1:3

			dP[a,i,j] = -Eext[a,i,j]

			R2Pt[a] = 0.0 

			
			for b in 1:3

				for c in 1:2, d in 1:2
					dP[a,i,j] -= -0.5*(G4[c,a,d,b]+G4[c,b,d,a])*DDPpt[b,c,d]
				end

				R2Pt[a] += R2[b,a]*Ppt[b]
			
				dP[a,i,j] -= 2.0*A2s[a,b]*Ppt[b]
				
				for c in 1:3, d in 1:3
					dP[a,i,j] -= 4.0*A4s[a,b,c,d]*Ppt[b]*Ppt[c]*Ppt[d]
				end

			end
			
		end

		dD6[1] = d1V6(R2Pt, a111, a112, a123)
		dD6[2] = d2V6(R2Pt, a111, a112, a123)
		dD6[3] = d3V6(R2Pt, a111, a112, a123)

		for a in 1:3, b in 1:3

			dP[a,i,j] -= R2[a,b]*dD6[b]
			
		end
		
	end
	
	
end





























function affectnoN!(integrator)

	Nx,Ny,dx,dy,G4,A2s,A4s, a111, a112, a123, C2,q3,Q3,R2, Eext = integrator.p[1]

	updateYbdy!(integrator.u, integrator.p[1][1], integrator.p[1][2])

	integrator.p[2][1] += 1

	#if integrator.p[2][1] % 1 == 9
	integrator.p[3][ integrator.p[2][1]  ] =  energy2D_noED(integrator.u, G4, A2s, A4s,a111,a112,a123,0.0,R2,Nx,Ny,dx,dy) 
	#end



	if integrator.p[3][ integrator.p[2][1]  ] > integrator.p[3][ integrator.p[2][1] - 1 ]
		for a in 1:3, i in 1:Nx, j in 1:Ny
			integrator.u[a+3,i,j] = 0.0
		end
	end

	

end

function flow2DRAKnoN!(P0,Pd0,tf,pars)
    

	counter = 1
	EList = zeros(10000)
	EList[1] = 1000000

	p = (pars,[counter],EList)

	u0 = vcat(P0,Pd0)

    tspan = (0.0,tf)
	
	cb = DiscreteCallback(conditionno,affectnoN!, save_positions=(false,false))

    prob = ODEProblem(dudf2DnoN!,u0,tspan,p)
	sol = solve(prob, Tsit5(), save_everystep=false, callback=cb)

	P0 .= sol[end][1:3,:,:]

    #u0 .= sol[end]

    return sol, EList[1:p[2][1]]
    
end





function dudf2DnoN!(dP,P,p,t)

	
	Nx,Ny,dx,dy,G4,A2s,A4s, a111, a112, a123, C2,q3,Q3,R2, Eext = p[1]

	Ppt = zeros(3)
	DDPpt = zeros(3,2,2)

	dD6 = zeros(3)
	R2Pt = zeros(3)
	
	
	@inbounds for i in 3:Nx-2, j in 3:Ny-2
		
		Ppt[1] = P[1,i,j] 
		Ppt[2] = P[2,i,j] 
		Ppt[3] = P[3,i,j] 

		getDDP!(DDPpt,P,i,j,dx,dy)

		for a in 1:6
			dP[a,i,j] = 0.0
		end
			
		for a in 1:3

			# THE OTHER ONE!!!
			dP[a,i,j] = P[a+3,i,j]

			dP[a+3,i,j] = -Eext[a,i,j]

			R2Pt[a] = 0.0 

			
			for b in 1:3

				for c in 1:2, d in 1:2
					dP[a+3,i,j] += 0.5*(G4[c,a,d,b]+G4[c,b,d,a])*DDPpt[b,c,d]
				end

				R2Pt[a] += R2[b,a]*Ppt[b]
			
				dP[a+3,i,j] -= 2.0*A2s[a,b]*Ppt[b]
				
				for c in 1:3, d in 1:3
					dP[a+3,i,j] -= 4.0*A4s[a,b,c,d]*Ppt[b]*Ppt[c]*Ppt[d]
				end

			end
			
		end

		dD6[1] = d1V6(R2Pt, a111, a112, a123)
		dD6[2] = d2V6(R2Pt, a111, a112, a123)
		dD6[3] = d3V6(R2Pt, a111, a112, a123)

		for a in 1:3, b in 1:3

			dP[a+3,i,j] -= R2[a,b]*dD6[b]
			
		end
		
	end
end