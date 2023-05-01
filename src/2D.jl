



function dEdP2Dsfel!(dedp,P,ϵ,Nx,Ny,dx,dy,G4,A2s,A4s,q3,R2, a111, a112, a123)

	
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

			dedp[a,i,j] = 0.0

			R2Pt[a] = 0.0 


			for b in 1:3

				for α in 1:6
					dedp[a,i,j] -= (q3[α,a,b]+q3[α,b,a])*ept[α]*Ppt[b]
				end

				for c in 1:2, d in 1:2
					dedp[a,i,j] += -0.5*(G4[c,a,d,b]+G4[c,b,d,a])*DDPpt[b,c,d]
				end

				R2Pt[a] += R2[b,a]*Ppt[b]
			
				dedp[a,i,j] += 2.0*A2s[a,b]*Ppt[b]
				
				for c in 1:3, d in 1:3
					dedp[a,i,j] += 4.0*A4s[a,b,c,d]*Ppt[b]*Ppt[c]*Ppt[d]
				end

			end
			
		end

		dD6[1] = d1V6(R2Pt, a111, a112, a123)
		dD6[2] = d2V6(R2Pt, a111, a112, a123)
		dD6[3] = d3V6(R2Pt, a111, a112, a123)

		for a in 1:3, b in 1:3

			dedp[a,i,j] += R2[a,b]*dD6[b]
			
		end
		
	end
	
	
end




function flow2Del!(P0,ϵ0,Nloops,dt,pars,elasticstuff;testenergy=false,compat=true,Eplot=false)


	if Eplot == true
		fig =  Figure(resolution = (300, 300))
		display(fig)
		ax = Axis(fig[1,1])
		t0 = 0.0
	end

	cc00,CL,ncsb,Lx,Ly,pmax,qmax,Cabm = elasticstuff
	Nx,Ny,dx,dy,G4,A2s,A4s, a111, a112, a123, C2,q3,Q3,R2 = pars

	QPP = zeros(6,Nx,Ny)
	getQPP!(QPP,Q3,P0,Nx,Ny)
	dedp = zeros(3,Nx,Ny)
	ϵbasis = zeros(8,3,Nx,Ny)

	if compat == true
		get_ϵbestL!(ϵ0,ϵbasis,QPP,cc00,CL,ncsb,Nx,Ny,dx,dy,Lx,Ly,pmax,qmax,Cabm)
	else
		get_ϵbestL_no_compat!(ϵ0,QPP,Nx,Ny)
	end

	if testenergy == true
		ED = zeros(Nx,Ny)
		println("initial energy is ", energy2Del(P0, ϵ0, ED, G4, A2s, A4s, C2, q3, a111,a112,a123,0.0,R2,Nx, Ny, dx, dy) )
	end

	for c in 1:Nloops


    	dEdP2Dsfel!(dedp,P0,ϵ0,Nx,Ny,dx,dy,G4,A2s,A4s,q3,R2, a111, a112, a123)
   		
   		for a in 1:3, i in 1:Nx, j in 1:Ny
   			P0[a,i,j] -= dt*dedp[a,i,j]
   		end

   		#P0 .-= dt.*dedp

   		updateYbdy!(P0,Nx,Ny)

   		getQPP!(QPP,Q3,P0,Nx,Ny)

   		if compat == true
			get_ϵbestL!(ϵ0,ϵbasis,QPP,cc00,CL,ncsb,Nx,Ny,dx,dy,Lx,Ly,pmax,qmax,Cabm)
		else
			get_ϵbestL_no_compat!(ϵ0,QPP,Nx,Ny)
		end

		if Eplot == true && c % 10 == 1
			scatter!(ax,Float64(c),energy2Del_noED(P0, ϵ0, G4, A2s, A4s, C2, q3, a111,a112,a123,0.0,R2,Nx, Ny, dx, dy) )
			display(fig)
		end
	

    	

    end

    if testenergy == true
    	ED = zeros(Nx,Ny)
		println("final   energy is ",  energy2Del(P0, ϵ0, ED, G4, A2s, A4s, C2, q3, a111,a112,a123,0.0,R2,Nx, Ny, dx, dy) )
	end

end





function flow2D!(P0,Nloops,dt,pars,Eext;testenergy=false,Eplot=false, boundary="periodic")


	if Eplot == true
		fig =  Figure(resolution = (300, 300))
		display(fig)
		ax = Axis(fig[1,1])
		t0 = 0.0
	end

	if boundary == "periodic"
    	println("you're using periodic boundary conditions")
	elseif boundary == "neumann"
		println("you're using Neumann boundary conditions")
	end

	Nx,Ny,dx,dy,G4,A2s,A4s,R2, a111, a112, a123 = pars

		dedp = zeros(3,Nx,Ny);

	if testenergy == true
		ED = zeros(Nx,Ny)
		println("initial energy is ", energy2D(P0, ED, G4, A2s, A4s, a111,a112,a123,0.0,R2,Nx, Ny, dx, dy) )
	end

	for c in 1:Nloops

    	dEdP2Dsf!(dedp,P0,Nx,Ny,dx,dy,G4,A2s,A4s,R2, a111, a112, a123,Eext)
    	P0 .-= dt.*dedp

    	if boundary == "periodic"
    		updateYbdy!(P0,Nx,Ny)
    	elseif boundary == "neumann"
    		updateYneu!(P0,Nx,Ny)
    	end


		if Eplot == true && c % 1000 == 1
			scatter!(ax,Float64(c),energy2D(P0, ED, G4, A2s, A4s, a111,a112,a123,0.0,R2,Nx, Ny, dx, dy) )
			display(fig)
		end

    end

    if testenergy == true
    	ED = zeros(Nx,Ny)
		println("final energy is ", energy2D(P0, ED, G4, A2s, A4s, a111,a112,a123,0.0,R2,Nx, Ny, dx, dy) )
	end

end




function makeSkyrmeLine(Nx,Ny,dx,dy,x,y,Lx,Rvac,PV)

	P0 = zeros(3,Nx,Ny)
	Pin = zeros(3,Nx,Ny)

	for a in 1:3, i in 1:Nx, j in 1:Ny
    	P0[a,i,j] = PV[a]
	end

	f(x) = pi*( tanh(x) + 1.0 )/2.0

	for i in 1:Nx, j in 1:Ny
    
    	P0[1,i,j] = 0.0
    	P0[2,i,j] = 0.0
    	P0[3,i,j] = 0.0
    
    	Pin[1,i,j] = sin(f(y[j]))*cos(x[i]*2.0*pi/Lx)
    	Pin[2,i,j] = sin(f(y[j]))*sin(x[i]*2.0*pi/Lx)
    	Pin[3,i,j] = cos(f(y[j]))
        
    	for a in 1:3, b in 1:3
        	P0[a,i,j] += Rvac[a,b]*Pin[b,i,j]
    	end
    
	end

	updateYbdy!(P0,Nx,Ny)

	return P0

end


function getlevic()

	levic = zeros(3,3,3)

	levic[1,2,3] = 1.0
	levic[3,1,2] = 1.0
	levic[2,3,1] = 1.0
	levic[1,3,2] = -1.0
	levic[2,1,3] = -1.0
	levic[3,2,1] = -1.0

	return levic

end



function baryon(P0, x,y ; BD = zeros(size(x)[1],size(y)[1]) )

	dx = step(x)
	dy = step(y)

	Nx = size(x)[1]
	Ny = size(y)[1]

	levic = getlevic()

	Ppt = zeros(3)
	DPpt = zeros(3,2)

	unitP = zeros(3,Nx,Ny)
	for i in 3:Nx-2, j in 3:Ny-2
		unitP[:,i,j] = normalize( [P0[1,i,j], P0[2,i,j], P0[3,i,j] ])
	end

	
    for i in 3:Nx-2, j in 3:Ny-2
        
        BD[i-2,j-2] = 0.0
        
        Ppt =  [unitP[1,i,j], unitP[2,i,j], unitP[3,i,j] ]
        getDP2D!(DPpt,unitP,i,j,dx,dy)

        for a in 1:3, b in 1:3, c in 1:3

        	BD[i-2,j-2] += levic[a,b,c]*Ppt[a]*DPpt[b,1]*DPpt[c,2]

        end

    end

    return sum(BD)*dx*dy/(4.0*pi)

end



function setgrid2D(Nx,Ny,dx,dy)

	#println("grid has ", Nx, "x", Ny" points. \n The x coord and goes from ", -0.5*dx*(Nx-1), " to ", 0.5*dx*(Nx-1), " while y goes from ", -0.5*dy*(Ny-1), " to ", 0.5*dy*(Ny-1) )	
	return -0.5*dx*(Nx-1):dx:0.5*dx*(Nx-1), -0.5*dy*(Ny-1):dy:0.5*dx*(Ny-1)
	
end


function getDDu!(DDupt,u,i,j,dx,dy)

	for a in 1:2, b in 1:2
		DDupt[a,b,1,1] = d2xDu(u,a,b,i,j,dx)
		DDupt[a,b,2,2] = d2yDu(u,a,b,i,j,dy)
		DDupt[a,b,1,2] = (dxdydiffDu(u, a, b, i, j, dx, dy) - DDupt[a,b,1,1] - DDupt[a,b,2,2])/2.0
		DDupt[a,b,2,1] = DDupt[a,b,1,2]
	end

end




function engpt2D(Ppt,DPpt,R2Tpt,R2,G4,A2,A4,a111,a112,a123,V0)
	
	engt = 0.0

	for a in 1:3, b in 1:3

		engt += A2[a,b]*Ppt[a]*Ppt[b]
		
		for c in 1:3, d in 1:3
			engt += A4[a,b,c,d]*Ppt[a]*Ppt[b]*Ppt[c]*Ppt[d]
		end
		
		for c in 1:2, d in 1:2

			engt += 0.5*G4[c,a,d,b]*DPpt[a,c]*DPpt[b,d]

		end

	end

	engt += V6o(R2'*Ppt,a111,a112,a123)
	engt -= V0

	return engt
	
end




function energy2D_noED(P, G4, A2, A4,a111,a112,a123,V0,R2,Nx,Ny,dx,dy)


	Ppt = zeros(3)
	DPpt = zeros(3,2)
	R2Tpt = zeros(3)
	engtot = 0.0
	
	@inbounds for i in 3:Nx-2, j in 3:Ny-2
		
			#Ppt = getP2D(P,i,j)

			Ppt[1] = P[1,i,j]
			Ppt[2] = P[2,i,j]
			Ppt[3] = P[3,i,j]

			R2Tpt[1] = 0.0
			R2Tpt[2] = 0.0
			R2Tpt[3] = 0.0

			getDP2D!(DPpt,P,i,j,dx,dy)
			
			engtot += engpt2D(Ppt,DPpt,R2Tpt,R2,G4,A2,A4,a111,a112,a123,V0)
	end
	
	return engtot*dx*dy
	
end



function energy2D(P, ED, G4, A2, A4,a111,a112,a123,V0,R2,Nx,Ny,dx,dy)


	Ppt = zeros(3)
	DPpt = zeros(3,2)
	R2Tpt = zeros(3)
	
	@inbounds for i in 3:Nx-2, j in 3:Ny-2
		
			#Ppt = getP2D(P,i,j)

			Ppt[1] = P[1,i,j]
			Ppt[2] = P[2,i,j]
			Ppt[3] = P[3,i,j]

			R2Tpt[1] = 0.0
			R2Tpt[2] = 0.0
			R2Tpt[3] = 0.0

			getDP2D!(DPpt,P,i,j,dx,dy)
			
			ED[i,j] = engpt2D(Ppt,DPpt,R2Tpt,R2,G4,A2,A4,a111,a112,a123,V0)
	end
	
	return sum(ED)*dx*dy
	
end




function engpt2Del(Ppt,ept,DPpt,R2Tpt,R2,G4,A2,A4,C2,q3,a111,a112,a123,V0)
	
	engpt = 0.0

	for a in 1:3, b in 1:3

		engpt += A2[a,b]*Ppt[a]*Ppt[b]
		
		for c in 1:3, d in 1:3
			engpt += A4[a,b,c,d]*Ppt[a]*Ppt[b]*Ppt[c]*Ppt[d]
		end
		
		for c in 1:2, d in 1:2

			engpt += 0.5*G4[c,a,d,b]*DPpt[a,c]*DPpt[b,d]

		end

	end

	for a in 1:6, b in 1:6

		engpt += 0.5*C2[a,b]*ept[a]*ept[b]

	end

	for a in 1:6

		for b in 1:6
			engpt += 0.5*C2[a,b]*ept[a]*ept[b]
		end

		for b in 1:3, c in 1:3
			engpt -= q3[a,b,c]*ept[a]*Ppt[b]*Ppt[c]
		end

	end

	engpt += V6o(R2'*Ppt,a111,a112,a123)
	engpt -= V0

	return engpt
	
end


function energy2Del(P, ϵ, ED, G4, A2, A4, C2, q3, a111,a112,a123,V0,R2,Nx,Ny,dx,dy)


	Ppt = zeros(3)
	ept = zeros(6)
	DPpt = zeros(3,2)
	R2Tpt = zeros(3)
	
	@inbounds for i in 3:Nx-2
		for j in 3:Ny-2
		
			#Ppt = getP2D(P,i,j)

			for a in 1:3
				Ppt[a] = P[a,i,j]
				R2Tpt[a] = 0.0
			end
			for a in 1:6
				ept[a] = ϵ[a,i,j]
			end

			getDP2D!(DPpt,P,i,j,dx,dy)
			
			ED[i,j] = engpt2Del(Ppt,ept,DPpt,R2Tpt,R2,G4,A2,A4,C2,q3,a111,a112,a123,V0)

		end
	end
	
	return sum(ED)*dx*dy
	
end


function energy2Del_noED(P, ϵ, G4, A2, A4, C2, q3, a111,a112,a123,V0,R2,Nx,Ny,dx,dy)


	Ppt = zeros(3)
	ept = zeros(6)
	DPpt = zeros(3,2)
	R2Tpt = zeros(3)

	etot = 0.0
	
	@inbounds for i in 3:Nx-2
		for j in 3:Ny-2
		
			#Ppt = getP2D(P,i,j)

			for a in 1:3
				Ppt[a] = P[a,i,j]
				R2Tpt[a] = 0.0
			end
			for a in 1:6
				ept[a] = ϵ[a,i,j]
			end

			getDP2D!(DPpt,P,i,j,dx,dy)
			
			etot += engpt2Del(Ppt,ept,DPpt,R2Tpt,R2,G4,A2,A4,C2,q3,a111,a112,a123,V0)

		end
	end
	
	return etot*dx*dy
	
end










function dEdP2Ds!(dedp,P,Nx,Ny,dx,dy,G4,A2s,A4s,R2, a111, a112, a123)

	
	Ppt = zeros(3)
	DDPpt = zeros(3,2,2)

	dD6 = zeros(3)
	R2Pt = zeros(3)
	
	
	for i in 3:Nx-2, j in 3:Ny-2
		
		Ppt[1] = P[1,i,j] 
		Ppt[2] = P[2,i,j] 
		Ppt[3] = P[3,i,j] 

		getDDP!(DDPpt,P,i,j,dx,dy)
			
		for a in 1:3

			dedp[a,i,j] = 0.0

			R2Pt[a] = 0.0 

			
			for b in 1:3

				for c in 1:2, d in 1:2
					dedp[a,i,j] += -0.5*(G4[c,a,d,b] + G4[c,b,d,a])*DDPpt[b]
				end

				R2Pt[a] += R2[b,a]*Ppt[b]
			
				
				
				dedp[a,i,j] += A2s[a,b]*Ppt[b]
				
				for c in 1:3, d in 1:3
					dedp[a,i,j] += A4s[a,b,c,d]*Ppt[b]*Ppt[c]*Ppt[d]
				end

			end
			
		end

		#dD6 = dV6(R2Pt,a111,a112,a123)
		dD6[1] = d1V6(R2Pt, a111, a112, a123)
		dD6[2] = d2V6(R2Pt, a111, a112, a123)
		dD6[3] = d3V6(R2Pt, a111, a112, a123)

		for a in 1:3, b in 1:3

			dedp[a,i,j] += R2[a,b]*dD6[b]
			
		end
		
	end
	
	
end




function dEdP2Dsf!(dedp,P,Nx,Ny,dx,dy,G4,A2s,A4s,R2, a111, a112, a123,Eext)

	
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

			dedp[a,i,j] = Eext[a,i,j]

			R2Pt[a] = 0.0 

			
			for b in 1:3

				for c in 1:2, d in 1:2
					dedp[a,i,j] += -0.5*(G4[c,a,d,b]+G4[c,b,d,a])*DDPpt[b,c,d]
				end

				R2Pt[a] += R2[b,a]*Ppt[b]
			
				dedp[a,i,j] += 2.0*A2s[a,b]*Ppt[b]
				
				for c in 1:3, d in 1:3
					dedp[a,i,j] += 4.0*A4s[a,b,c,d]*Ppt[b]*Ppt[c]*Ppt[d]
				end

			end
			
		end

		dD6[1] = d1V6(R2Pt, a111, a112, a123)
		dD6[2] = d2V6(R2Pt, a111, a112, a123)
		dD6[3] = d3V6(R2Pt, a111, a112, a123)

		for a in 1:3, b in 1:3

			dedp[a,i,j] += R2[a,b]*dD6[b]
			
		end
		
	end
	
	
end


function updateYbdy!(P0,Nx,Ny)

    @inbounds for a in 1:3, j in 1:Ny

        P0[a,1,j] = P0[a,Nx-3,j]
        P0[a,2,j] = P0[a,Nx-2,j]

        P0[a,Nx-1,j] = P0[a,3,j]
        P0[a,Nx,j] = P0[a,4,j]

    end
			
end


function updateYneu!(P0,Nx,Ny)

    @inbounds for a in 1:3, j in 1:Ny

        P0[a,1,j] = P0[a,3,j]
        P0[a,2,j] = P0[a,3,j]

        P0[a,Nx-1,j] = P0[a,Nx-2,j]
        P0[a,Nx,j] = P0[a,Nx-2,j]

    end
			
end






