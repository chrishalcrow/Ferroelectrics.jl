function dEdu!(du,u,p,t)
	
	# p = (N,dx,G2,A2,A4,R2, a111, a112, a123)
	
	#dedp = zeros(3,p[1])
	
	
	DDPpt = zeros(3)
	Ppt = zeros(3)
	
    for i in 3:p[1]-2
		for a in 1:3
			du[a,i] = 0.0
		end
		
		Ppt = getP(u,i)
		getDDP!(DDPpt,u,i,p[2])
		
		for a in 1:3, b in 1:3
			
			du[a,i] -= -p[3][a,b]*DDPpt[b]
			
			du[a,i] -= (p[4][a,b] + p[4][b,a])*Ppt[b]
			
			for c in 1:3, d in 1:3
				du[a,i] -= (p[5][a,b,c,d] + p[5][b,a,c,d] + p[5][b,c,a,d] + p[5][b,c,d,a])*Ppt[b]*Ppt[c]*Ppt[d]
			end
			
		end
			
		du[:,i] .-= p[6]*dV6(p[6]'*Ppt,p[7], p[8], p[9])

    end
end

function dEduf!(du,u,p,t)
	
	#N,dx,G2,A2,A4,R2, a111, a112, a123 = p
	
	#dedp = zeros(3,p[1])
	
	
	DDPpt = zeros(3)
	Ppt = zeros(3)
	dD6 = zeros(3)
	R2Pt = zeros(3)
	
    @inbounds for i in 3:p[1]-2
		
		Ppt[1] = u[1,i] 
		Ppt[2] = u[2,i] 
		Ppt[3] = u[3,i] 

		getDDP!(DDPpt,u,i,p[2])
		
		for a in 1:3
			
			du[a,i] = 0.0
			R2Pt[a] = 0.0 

			for  b in 1:3

				R2Pt[a] += p[6][b,a]*Ppt[b]

				du[a,i] -= -p[3][a,b]*DDPpt[b]
				
				du[a,i] -= 2.0*(p[4][a,b])*Ppt[b]
				
				for c in 1:3, d in 1:3
					du[a,i] -= 4.0*p[5][a,b,c,d]*Ppt[b]*Ppt[c]*Ppt[d]
				end

			end
			
		end

		dD6[1] = d1V6(R2Pt, p[7], p[8], p[9])
		dD6[2] = d2V6(R2Pt, p[7], p[8], p[9])
		dD6[3] = d3V6(R2Pt, p[7], p[8], p[9])

		for a in 1:3, b in 1:3

			du[a,i] -= p[6][a,b]*dD6[b]
			
		end
			

    end
end


function dEduff!(du,u,p,t)
	
	N,dx,G2,A2,A4,R2, a111, a112, a123,tf,EList,V0, DDPpt, Ppt, dD6, R2Pt = p
	

	
	
	#DDPpt = zeros(3)
	#Ppt = zeros(3)
	#dD6 = zeros(3)
	#R2Pt = zeros(3)
	
    for i in 3:N-2
		
		Ppt[1] = u[1,i] 
		Ppt[2] = u[2,i] 
		Ppt[3] = u[3,i] 

		getDDP!(DDPpt,u,i,p[2])
		
		for a in 1:3
			
			du[a,i] = 0.0
			R2Pt[a] = 0.0 

			for  b in 1:3

				R2Pt[a] += R2[b,a]*Ppt[b]

				du[a,i] -= -G2[a,b]*DDPpt[b]
				
				du[a,i] -= 2.0*(A2[a,b])*Ppt[b]
				
				for c in 1:3, d in 1:3
					du[a,i] -= 4.0*A4[a,b,c,d]*Ppt[b]*Ppt[c]*Ppt[d]
				end

			end
			
		end

		dD6[1] = d1V6(R2Pt, a111, a112, a123)
		dD6[2] = d2V6(R2Pt, a111, a112, a123)
		dD6[3] = d3V6(R2Pt, a111, a112, a123)

		for a in 1:3, b in 1:3

			du[a,i] -= R2[a,b]*dD6[b]
			
		end
			

    end
end



	


function affect_checkeng!(integrator)
    
    n = floor(Int,round(integrator.t))

    ED = zeros(integrator.p[1])
    
    integrator.p[11][n] = energy(integrator.u, ED, integrator.p[2], integrator.p[1], integrator.p[6], integrator.p[3], integrator.p[4], integrator.p[5], integrator.p[7], integrator.p[8], integrator.p[9], integrator.p[12]) 
    

    #println("energy at t = ", integrator.t, " is ", integrator.p[11][n] )

    if n > 3 && abs(integrator.p[11][n] - integrator.p[11][n-1]) < 10.0^(-8)
        terminate!(integrator)
        integrator.p[11][(n+1):integrator.p[10]] .= integrator.p[11][n]
    end

end



function flowRAK!(phi,p)
    
	tf = p[10]
    EList = zeros(tf)

	# p = (N,dx,G2,A2,A4,R2, a111, a112, a123,tf,EList,V0)

    tspan = (0.0, tf)
    engtimes = [ 1.0*i for i=1:tf ]

    cb = PresetTimeCallback(engtimes, affect_checkeng!, save_positions=(false,false))

    prob = ODEProblem(dEduf!,phi,tspan,p)
	sol = solve(prob, Tsit5(),reltol=1e-8, abstol=1e-8, save_everystep=false, callback=cb )
    #sol = solve(prob, Tsit5(),reltol=1e-8, abstol=1e-8, save_everystep=true, callback=cb, tstops = engtimes )

    #if sol.t[end] == tf
    #   println("AGGGG");
    #else
    #    println(sol.t[end])
    #end

    return sol[end]
    
end


function flowRAKf!(phi,q)
    
	#ED = zeros(integrator.p[1])

	
    #EList = zeros(tf)

    DDPpt = zeros(3)
    Ppt = zeros(3)
    dD6 = zeros(3)
    R2Pt = zeros(3)

	p = (q[1], q[2], q[3], q[4], q[5], q[6], q[7], q[8], q[9], q[10], q[11], q[12], DDPpt, Ppt, dD6, R2Pt)

	tf = p[10]

	# p = (N,dx,G2,A2,A4,R2, a111, a112, a123,tf,EList,V0, ED)

    tspan = (0.0, tf)
    engtimes = [ 1.0*i for i=1:tf ]

    cb = PresetTimeCallback(engtimes, affect_checkengf!, save_positions=(false,false))

    prob = ODEProblem(dEduff!,phi,tspan,p)
	sol = solve(prob, Tsit5(),reltol=1e-8, abstol=1e-8, save_everystep=false, callback=cb )
    #sol = solve(prob, Tsit5(),reltol=1e-8, abstol=1e-8, save_everystep=true, callback=cb, tstops = engtimes )

    if sol.t[end] == tf
       println("AGGGG");
    else
        println(sol.t[end])
    end

    return sol[end]
    
end



function affect_checkengf!(integrator)
    
    n = floor(Int,round(integrator.t))

    integrator.p[11][n] = energyNOED(integrator.u, integrator.p[2], integrator.p[1], integrator.p[6], integrator.p[3], integrator.p[4], integrator.p[5], integrator.p[7], integrator.p[8], integrator.p[9], integrator.p[12]) 


    if n > 3 && abs(integrator.p[11][n] - integrator.p[11][n-1]) < 10.0^(-8)
        terminate!(integrator)
        integrator.p[11][(n+1):integrator.p[10]] .= integrator.p[11][n]
    end

end



