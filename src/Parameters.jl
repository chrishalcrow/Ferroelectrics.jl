function δ2(i,j)
    if i == j
        return 1.0
    else
        return 0.0
    end
end

function δ4(i,j,k,l)
    if i == j && j ==k && k==l
        return 1.0
    else
        return 0.0
    end
end

function δ6(i,j,k,l,m,n)
    if i == j && j ==k && k==l && l==m && m==n
        return 1.0
    else
        return 0.0
    end
end

function get_As(A2,A4)

	A2s = zeros(3,3)
	A4s = zeros(3,3,3,3)

	for a in 1:3, b in 1:3
	    A2s[a,b] = 0.5*(A2[a,b] + A2[b,a])
	    for c in 1:3, d in 1:3
	        A4s[a,b,c,d] = 0.25*(A4[a,b,c,d] + A4[b,a,c,d] + A4[b,c,a,d] + A4[b,c,d,a])
	    end
	end

	return (A2,A4)

end


function set_parameters!(G2,A2t,A4t,T,R2,phase;  im2D=false)
	
	#GLMakie.activate()

	# VOIGT TRANSFORM 
	TT = zeros(6,3,3)
	PV = zeros(3)

	TT[1,1,1] = 1.0
	TT[1+1,1+1,1+1] = 1.0
	TT[2+1,2+1,2+1] = 1.0

	TT[3+1,1+1,2+1] = 0.5
	TT[3+1,2+1,1+1] = 0.5

	TT[4+1,0+1,2+1] = 0.5
	TT[4+1,2+1,0+1] = 0.5

	TT[5+1,1+1,0+1] = 0.5
	TT[5+1,0+1,1+1] = 0.5

	a1 = 334000*(-381 + T)
	a11 = -1.2275565521327014e9 + 4.69e6*T
	a12 = -3.4411355085840464e8
	a111 = 2.76e9 - 5.52e7*(-393 + T)
	a112 = 4.47e9
	a123 = 4.91e9

	C11 = 2.75e11
	C12 = 1.79e11
	C44 = 5.43e10

	q11 = 1.42e10
	q12 = -7.4e8
	q44 = 1.57e9

	G11 = 5.1e-10
	G12 = -2.e-11
	G44 = 2.e-11

	A0t = 0.0

	A2 = zeros(3,3)
	A2o = zeros(3,3)
	#A2t = zeros(3,3)
	
	A4 = zeros(3,3,3,3)
	A4o = zeros(3,3,3,3)
	A4ot = zeros(3,3,3,3)
	A4b = zeros(3,3,3,3)
	#A4t = zeros(3,3,3,3)
	
	Q4 = zeros(3,3,3,3)
	Q4o = zeros(3,3,3,3)
	Q3 = zeros(6,3,3)
	Q3o = zeros(6,3,3)
	
	q4 = zeros(3,3,3,3)
	q4o = zeros(3,3,3,3)
	q3 = zeros(6,3,3)
	
	G4o = zeros(3,3,3,3)
	G4 = zeros(3,3,3,3)

	
	C4 = zeros(3,3,3,3)
	C4o = zeros(3,3,3,3)
	C2 = zeros(6,6)
	
	for a in 1:3, b in 1:3
	    A2o[a,b] = a1*δ2(a,b)
	    for c in 1:3, d in 1:3
	        A4o[a,b,c,d] = (a11 - 0.5*a12)*δ4(a,b,c,d) + a12/6.0*(δ2(a,b)*δ2(c,d) + δ2(a,c)*δ2(b,d) + δ2(a,d)*δ2(b,c))
	        C4o[a,b,c,d] = (C11 - C12 - 2.0*C44)*δ4(a,b,c,d) + C12*(δ2(a,b)*δ2(c,d)) + C44*(δ2(a,c)*δ2(b,d) + δ2(a,d)*δ2(b,c))
			G4o[a,b,c,d] = (G11 - G12 - 2.0*G44)*δ4(a,b,c,d) + G12*(δ2(a,b)*δ2(c,d)) + G44*(δ2(a,c)*δ2(b,d) + δ2(a,d)*δ2(b,c))
	        q4o[a,b,c,d] = (q11 - q12 - q44)*δ4(a,b,c,d) + q12*(δ2(a,b)*δ2(c,d)) + q44/2.0*(δ2(a,c)*δ2(b,d) + δ2(a,d)*δ2(b,c))
	    end
	end

	C2o = [ C11 C12 C12 0 0 0 ; 
	        C12 C11 C12 0 0 0 ;
	        C12 C12 C11 0 0 0 ;
	        0 0 0 C44 0 0 ;
	        0 0 0 0 C44 0 ;
	        0 0 0 0 0 C44 ]
			
	C2oi = inv(C2o)

	q3o = zeros(6,3,3)

	for a in 1:6, b in 1:3, c in 1:3, d in 1:3, e in 1:3
	    q3o[a,b,c] += TT[a,d,e]*q4o[d,e,b,c]
	end

	for a in 1:6, b in 1:3, c in 1:3, d in 1:6
	    Q3o[a,b,c] += C2oi[a,d]*q3o[d,b,c]
	end

	for e in 1:3, f in 1:3, g in 1:3, h in 1:3
	    A4ot[e,f,g,h] = A4o[e,f,g,h]
	    for α in 1:6, β in 1:6 
	        A4ot[e,f,g,h] -= 0.5*C2o[α,β]*Q3o[α,e,f]*Q3o[β,g,h]
	    end
	end
	

	if phase == "rho"
		fr(p) = Vo([p,p,p]./sqrt(3.0),A2o,A4ot,a111,a112,a123)
		Dr(fr) = p -> ForwardDiff.derivative(fr,float(p))
		p0 = find_zero( Dr(fr), (0.1,1.5) )
		PV .= [p0,p0,p0]./sqrt(3.0)
	elseif phase == "oct"
		fo(p) = Vo([p,p,0.0]./sqrt(2.0),A2o,A4ot,a111,a112,a123)
		Do(fo) = p -> ForwardDiff.derivative(fo,float(p))
		p0 = find_zero( Do(fo), (0.1,1.5) )
		PV .= [p0,p0,0.0]./sqrt(2.0)
	elseif phase == "tet"
		ft(p) = Vo([p,0.0,0.0],A2o,A4ot,a111,a112,a123)
		Dt(ft) = p -> ForwardDiff.derivative(ft,float(p))
		p0 = find_zero( Dt(ft), (0.1,1.5) )
		PV .= [p0,0.0,0.0]
	else
		println(" NO ALLOWED PHASE GIVEN ERROR ERROR. ")
	end
	
	
	println("The ground state is PV = ", PV)
	
	A2o ./= abs(a1)
	A4o ./= abs(a1)/p0^2
	A4ot ./= abs(a1)/p0^2
	q4o ./= abs(a1)
	G4o /= G11
	Q4o .*= p0^2
	C4o ./= abs(a1)*p0^2
	
	PV ./= norm(PV)
	
	a111 /= abs(a1)/p0^4
	a112 /= abs(a1)/p0^4
	a123 /= abs(a1)/p0^4
	
	A2 = zeros(3,3)
	
	PV = R2*PV
	 
	for a in 1:3, b in 1:3, c in 1:3, d in 1:3
	    A2[a,b] += R2[a,c]*R2[b,d]*A2o[c,d]
	    A4b[a,b,c,d] = 0.0
	    for e in 1:3, f in 1:3, g in 1:3, h in 1:3
	    	A4b[a,b,c,d] += R2[a,e]*R2[b,f]*R2[c,g]*R2[d,h]*A4ot[e,f,g,h]
	        A4[a,b,c,d] += R2[a,e]*R2[b,f]*R2[c,g]*R2[d,h]*A4o[e,f,g,h]
	        C4[a,b,c,d] += R2[a,e]*R2[b,f]*R2[c,g]*R2[d,h]*C4o[e,f,g,h]
			G4[a,b,c,d] += R2[a,e]*R2[b,f]*R2[c,g]*R2[d,h]*G4o[e,f,g,h]
	        Q4[a,b,c,d] += R2[a,e]*R2[b,f]*R2[c,g]*R2[d,h]*Q4o[e,f,g,h]
	        q4[a,b,c,d] += R2[a,e]*R2[b,f]*R2[c,g]*R2[d,h]*q4o[e,f,g,h]
	    end
	end
	
	
	
	for a in 1:3, b in 1:3
		G2[a,b] = G4[1,a,1,b]
	end


	for a in 1:6, b in 1:3, c in 1:3, d in 1:3, e in 1:3
	    q3[a,b,c] += TT[a,d,e]*q4[d,e,b,c]
	    Q3[a,b,c] += TT[a,d,e]*Q4[d,e,b,c]
	end

	for a in 1:6, b in 1:6, c in 1:3, d in 1:3, e in 1:3, f in 1:3
	    C2[a,b] += TT[a,c,d]*TT[b,e,f]*C4[c,d,e,f]
	end



	for α in 1:6, β in 1:6, b in 1:3, c in 1:3
	    Q3[α,b,c] += inv(C2)[α,β]*q3[β,b,c]
	end

	qll = zeros(6,3,3)
	qperp = zeros(6,3,3)
	
	Qll = zeros(6,3,3)
	Qperp = zeros(6,3,3)
	
	Cll = zeros(6,6)
	Cperp = zeros(6,6)
	Cm = zeros(6,6)

	qperp[2,:,:] .= q3[2,:,:]
	qperp[3,:,:] .= q3[3,:,:]
	qperp[4,:,:] .= q3[4,:,:]

	qll[1,:,:] .= q3[1,:,:]
	qll[5,:,:] .= q3[5,:,:]
	qll[6,:,:] .= q3[6,:,:]

	Qperp[2,:,:] .= Q3[2,:,:]
	Qperp[3,:,:] .= Q3[3,:,:]
	Qperp[4,:,:] .= Q3[4,:,:]

	Qll[1,:,:] .= Q3[1,:,:]
	Qll[5,:,:] .= Q3[5,:,:]
	Qll[6,:,:] .= Q3[6,:,:]

	for a in 1:6, b in 1:6
    
	    if a in [1,5,6] && b in [1,5,6]
	        Cll[a,b] = C2[a,b]
	    elseif a in [2,3,4] && b in [2,3,4]
	        Cperp[a,b] = C2[a,b]
	    else
	        Cm[a,b] = C2[a,b]
	    end
		
	end

	Cll33 = [ Cll[1,1] Cll[1,5] Cll[1,6] ;
		      Cll[5,1] Cll[5,5] Cll[5,6] ;
		      Cll[6,1] Cll[6,5] Cll[6,6] ]

	Cll33inv = inv(Cll33)

	Cll66inv=[ Cll33inv[1,1] 0 0 0 Cll33inv[1,2] Cll33inv[1,3] ;
				0 0 0 0 0 0 ;
				0 0 0 0 0 0 ;
				0 0 0 0 0 0 ;
			    Cll33inv[2,1] 0 0 0 Cll33inv[2,2] Cll33inv[2,3] ;
			    Cll33inv[3,1]  0 0 0 Cll33inv[3,2] Cll33inv[3,3] ]

	ϵV = zeros(6)

	for α in 1:6, b in 1:3, c in 1:3
	    ϵV[α] += Q3[α,b,c]*PV[b]*PV[c]
	end

	H = zeros(6,3,3)
	K = zeros(6)
	
	

	for α in [1,5,6], b in 1:3, c in 1:3
	    H[α,b,c] = Qll[α,b,c]
	    for β in 1:6, γ in 1:6
	        H[α,b,c] += Cll66inv[α,β]*Cm[β,γ]*Qperp[γ,b,c]
	    end
	end


	for α in [1,5,6], β in 1:6, γ in 1:6
	    K[α] += -Cll66inv[α,β]*Cm[β,γ]*ϵV[γ]
	end
	
	for α in [2,3,4]
	    K[α] = ϵV[α]
	end
		
	for a in 1:3, b in 1:3, c in 1:3, d in 1:3
		A4t[a,b,c,d] = A4[a,b,c,d]
		for alp in 1:6
			A4t[a,b,c,d] -= 1.0*q3[alp,a,b]*H[alp,c,d]
			for bet in 1:6
				A4t[a,b,c,d] += 0.5*H[bet,a,b]*C2[bet,alp]*H[alp,c,d]
			end
		end
	end
	

	
	for a in 1:3, b in 1:3
		A2t[a,b] = A2[a,b]
		for alp in 1:6
			A2t[a,b] -= q3[alp,a,b]*K[alp]
			for bet in 1:6
				A2t[a,b] += 0.5*(C2[alp,bet]+C2[bet,alp])*H[bet,a,b]*K[alp]
			end
		end
	end
	
	for alp in 1:6, bet in 1:6
		A0t += 0.5*K[alp]*C2[alp,bet]*K[bet]
	end
	

	
	println(PV)
	
	V0 = PotEng(PV,R2,G2,A2t,A4t,a111,a112,a123,0.0)
	

	
	println("the zero point energy is: ", PotEng(PV,R2,G2,A2t,A4t,a111,a112,a123,V0))
	
	if im2D == true 
		return PV, A0t, a111, a112, a123, V0, G4, A4b, C2, Q3, q3, A4
	else
		return PV, A0t, a111, a112, a123, V0
	end
	
end

function outputA6(a111,a112,a123,R2)

	A6o = zeros(3,3,3,3,3,3)
	A6 = zeros(3,3,3,3,3,3)


	A6o[1,1,1,1,1,1] = a111;
	A6o[1,1,1,1,2,2] = a112/15.;
	A6o[1,1,1,1,3,3] = a112/15.;
	A6o[1,1,1,2,1,2] = a112/15.;
	A6o[1,1,1,2,2,1] = a112/15.;
	A6o[1,1,1,3,1,3] = a112/15.;
	A6o[1,1,1,3,3,1] = a112/15.;
	A6o[1,1,2,1,1,2] = a112/15.;
	A6o[1,1,2,1,2,1] = a112/15.;
	A6o[1,1,2,2,1,1] = a112/15.;
	A6o[1,1,2,2,2,2] = a112/15.;
	A6o[1,1,2,2,3,3] = a123/90.;
	A6o[1,1,2,3,2,3] = a123/90.;
	A6o[1,1,2,3,3,2] = a123/90.;
	A6o[1,1,3,1,1,3] = a112/15.;
	A6o[1,1,3,1,3,1] = a112/15.;
	A6o[1,1,3,2,2,3] = a123/90.;
	A6o[1,1,3,2,3,2] = a123/90.;
	A6o[1,1,3,3,1,1] = a112/15.;
	A6o[1,1,3,3,2,2] = a123/90.;
	A6o[1,1,3,3,3,3] = a112/15.;
	A6o[1,2,1,1,1,2] = a112/15.;
	A6o[1,2,1,1,2,1] = a112/15.;
	A6o[1,2,1,2,1,1] = a112/15.;
	A6o[1,2,1,2,2,2] = a112/15.;
	A6o[1,2,1,2,3,3] = a123/90.;
	A6o[1,2,1,3,2,3] = a123/90.;
	A6o[1,2,1,3,3,2] = a123/90.;
	A6o[1,2,2,1,1,1] = a112/15.;
	A6o[1,2,2,1,2,2] = a112/15.;
	A6o[1,2,2,1,3,3] = a123/90.;
	A6o[1,2,2,2,1,2] = a112/15.;
	A6o[1,2,2,2,2,1] = a112/15.;
	A6o[1,2,2,3,1,3] = a123/90.;
	A6o[1,2,2,3,3,1] = a123/90.;
	A6o[1,2,3,1,2,3] = a123/90.;
	A6o[1,2,3,1,3,2] = a123/90.;
	A6o[1,2,3,2,1,3] = a123/90.;
	A6o[1,2,3,2,3,1] = a123/90.;
	A6o[1,2,3,3,1,2] = a123/90.;
	A6o[1,2,3,3,2,1] = a123/90.;
	A6o[1,3,1,1,1,3] = a112/15.;
	A6o[1,3,1,1,3,1] = a112/15.;
	A6o[1,3,1,2,2,3] = a123/90.;
	A6o[1,3,1,2,3,2] = a123/90.;
	A6o[1,3,1,3,1,1] = a112/15.;
	A6o[1,3,1,3,2,2] = a123/90.;
	A6o[1,3,1,3,3,3] = a112/15.;
	A6o[1,3,2,1,2,3] = a123/90.;
	A6o[1,3,2,1,3,2] = a123/90.;
	A6o[1,3,2,2,1,3] = a123/90.;
	A6o[1,3,2,2,3,1] = a123/90.;
	A6o[1,3,2,3,1,2] = a123/90.;
	A6o[1,3,2,3,2,1] = a123/90.;
	A6o[1,3,3,1,1,1] = a112/15.;
	A6o[1,3,3,1,2,2] = a123/90.;
	A6o[1,3,3,1,3,3] = a112/15.;
	A6o[1,3,3,2,1,2] = a123/90.;
	A6o[1,3,3,2,2,1] = a123/90.;
	A6o[1,3,3,3,1,3] = a112/15.;
	A6o[1,3,3,3,3,1] = a112/15.;
	A6o[2,1,1,1,1,2] = a112/15.;
	A6o[2,1,1,1,2,1] = a112/15.;
	A6o[2,1,1,2,1,1] = a112/15.;
	A6o[2,1,1,2,2,2] = a112/15.;
	A6o[2,1,1,2,3,3] = a123/90.;
	A6o[2,1,1,3,2,3] = a123/90.;
	A6o[2,1,1,3,3,2] = a123/90.;
	A6o[2,1,2,1,1,1] = a112/15.;
	A6o[2,1,2,1,2,2] = a112/15.;
	A6o[2,1,2,1,3,3] = a123/90.;
	A6o[2,1,2,2,1,2] = a112/15.;
	A6o[2,1,2,2,2,1] = a112/15.;
	A6o[2,1,2,3,1,3] = a123/90.;
	A6o[2,1,2,3,3,1] = a123/90.;
	A6o[2,1,3,1,2,3] = a123/90.;
	A6o[2,1,3,1,3,2] = a123/90.;
	A6o[2,1,3,2,1,3] = a123/90.;
	A6o[2,1,3,2,3,1] = a123/90.;
	A6o[2,1,3,3,1,2] = a123/90.;
	A6o[2,1,3,3,2,1] = a123/90.;
	A6o[2,2,1,1,1,1] = a112/15.;
	A6o[2,2,1,1,2,2] = a112/15.;
	A6o[2,2,1,1,3,3] = a123/90.;
	A6o[2,2,1,2,1,2] = a112/15.;
	A6o[2,2,1,2,2,1] = a112/15.;
	A6o[2,2,1,3,1,3] = a123/90.;
	A6o[2,2,1,3,3,1] = a123/90.;
	A6o[2,2,2,1,1,2] = a112/15.;
	A6o[2,2,2,1,2,1] = a112/15.;
	A6o[2,2,2,2,1,1] = a112/15.;
	A6o[2,2,2,2,2,2] = a111;
	A6o[2,2,2,2,3,3] = a112/15.;
	A6o[2,2,2,3,2,3] = a112/15.;
	A6o[2,2,2,3,3,2] = a112/15.;
	A6o[2,2,3,1,1,3] = a123/90.;
	A6o[2,2,3,1,3,1] = a123/90.;
	A6o[2,2,3,2,2,3] = a112/15.;
	A6o[2,2,3,2,3,2] = a112/15.;
	A6o[2,2,3,3,1,1] = a123/90.;
	A6o[2,2,3,3,2,2] = a112/15.;
	A6o[2,2,3,3,3,3] = a112/15.;
	A6o[2,3,1,1,2,3] = a123/90.;
	A6o[2,3,1,1,3,2] = a123/90.;
	A6o[2,3,1,2,1,3] = a123/90.;
	A6o[2,3,1,2,3,1] = a123/90.;
	A6o[2,3,1,3,1,2] = a123/90.;
	A6o[2,3,1,3,2,1] = a123/90.;
	A6o[2,3,2,1,1,3] = a123/90.;
	A6o[2,3,2,1,3,1] = a123/90.;
	A6o[2,3,2,2,2,3] = a112/15.;
	A6o[2,3,2,2,3,2] = a112/15.;
	A6o[2,3,2,3,1,1] = a123/90.;
	A6o[2,3,2,3,2,2] = a112/15.;
	A6o[2,3,2,3,3,3] = a112/15.;
	A6o[2,3,3,1,1,2] = a123/90.;
	A6o[2,3,3,1,2,1] = a123/90.;
	A6o[2,3,3,2,1,1] = a123/90.;
	A6o[2,3,3,2,2,2] = a112/15.;
	A6o[2,3,3,2,3,3] = a112/15.;
	A6o[2,3,3,3,2,3] = a112/15.;
	A6o[2,3,3,3,3,2] = a112/15.;
	A6o[3,1,1,1,1,3] = a112/15.;
	A6o[3,1,1,1,3,1] = a112/15.;
	A6o[3,1,1,2,2,3] = a123/90.;
	A6o[3,1,1,2,3,2] = a123/90.;
	A6o[3,1,1,3,1,1] = a112/15.;
	A6o[3,1,1,3,2,2] = a123/90.;
	A6o[3,1,1,3,3,3] = a112/15.;
	A6o[3,1,2,1,2,3] = a123/90.;
	A6o[3,1,2,1,3,2] = a123/90.;
	A6o[3,1,2,2,1,3] = a123/90.;
	A6o[3,1,2,2,3,1] = a123/90.;
	A6o[3,1,2,3,1,2] = a123/90.;
	A6o[3,1,2,3,2,1] = a123/90.;
	A6o[3,1,3,1,1,1] = a112/15.;
	A6o[3,1,3,1,2,2] = a123/90.;
	A6o[3,1,3,1,3,3] = a112/15.;
	A6o[3,1,3,2,1,2] = a123/90.;
	A6o[3,1,3,2,2,1] = a123/90.;
	A6o[3,1,3,3,1,3] = a112/15.;
	A6o[3,1,3,3,3,1] = a112/15.;
	A6o[3,2,1,1,2,3] = a123/90.;
	A6o[3,2,1,1,3,2] = a123/90.;
	A6o[3,2,1,2,1,3] = a123/90.;
	A6o[3,2,1,2,3,1] = a123/90.;
	A6o[3,2,1,3,1,2] = a123/90.;
	A6o[3,2,1,3,2,1] = a123/90.;
	A6o[3,2,2,1,1,3] = a123/90.;
	A6o[3,2,2,1,3,1] = a123/90.;
	A6o[3,2,2,2,2,3] = a112/15.;
	A6o[3,2,2,2,3,2] = a112/15.;
	A6o[3,2,2,3,1,1] = a123/90.;
	A6o[3,2,2,3,2,2] = a112/15.;
	A6o[3,2,2,3,3,3] = a112/15.;
	A6o[3,2,3,1,1,2] = a123/90.;
	A6o[3,2,3,1,2,1] = a123/90.;
	A6o[3,2,3,2,1,1] = a123/90.;
	A6o[3,2,3,2,2,2] = a112/15.;
	A6o[3,2,3,2,3,3] = a112/15.;
	A6o[3,2,3,3,2,3] = a112/15.;
	A6o[3,2,3,3,3,2] = a112/15.;
	A6o[3,3,1,1,1,1] = a112/15.;
	A6o[3,3,1,1,2,2] = a123/90.;
	A6o[3,3,1,1,3,3] = a112/15.;
	A6o[3,3,1,2,1,2] = a123/90.;
	A6o[3,3,1,2,2,1] = a123/90.;
	A6o[3,3,1,3,1,3] = a112/15.;
	A6o[3,3,1,3,3,1] = a112/15.;
	A6o[3,3,2,1,1,2] = a123/90.;
	A6o[3,3,2,1,2,1] = a123/90.;
	A6o[3,3,2,2,1,1] = a123/90.;
	A6o[3,3,2,2,2,2] = a112/15.;
	A6o[3,3,2,2,3,3] = a112/15.;
	A6o[3,3,2,3,2,3] = a112/15.;
	A6o[3,3,2,3,3,2] = a112/15.;
	A6o[3,3,3,1,1,3] = a112/15.;
	A6o[3,3,3,1,3,1] = a112/15.;
	A6o[3,3,3,2,2,3] = a112/15.;
	A6o[3,3,3,2,3,2] = a112/15.;
	A6o[3,3,3,3,1,1] = a112/15.;
	A6o[3,3,3,3,2,2] = a112/15.;
	A6o[3,3,3,3,3,3] = a111;


	for a in 1:3, b in 1:3, c in 1:3, d in 1:3
	    for e in 1:3, f in 1:3, g in 1:3, h in 1:3
	    	for i in 1:3, j in 1:3, k in 1:3, l in 1:3

	    		A6[a,b,c,d,e,f] += R2[a,g]*R2[b,h]*R2[c,i]*R2[d,j]*R2[e,k]*R2[f,l]*A6o[g,h,i,j,k,l]

	    	end
	    end
	end


	return A6

end





