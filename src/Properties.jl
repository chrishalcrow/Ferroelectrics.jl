


function V6o(P,a111,a112,a123)

    a111*(P[1]^6 + P[2]^6 + P[3]^6) + a112*(P[1]^4*(P[2]^2 + P[3]^2) + P[2]^4*(P[1]^2 + P[3]^2) + P[3]^4*(P[2]^2 + P[1]^2)) + a123*P[1]^2*P[2]^2*P[3]^2
    
end

function Vo(P, A2, A4, a111, a112, a123)
    
    v0 = 0.0
    
    for a in 1:3, b in 1:3
        v0 += A2[a,b]*P[a]*P[b]
        for c in 1:3, d in 1:3
            v0 += A4[a,b,c,d]*P[a]*P[b]*P[c]*P[d]
        end
    end
    
    v0 += V6o(P,a111,a112,a123)
    
    return v0
end

function V6_tensor(Ppt,A6)

	engt = 0.0

	for a in 1:3, b in 1:3, c in 1:3, d in 1:3, e in 1:3, f in 1:3
		engt += A6[a,b,c,d,e,f]*Ppt[a]*Ppt[b]*Ppt[c]*Ppt[d]*Ppt[e]*Ppt[f];
	end

	return engt

end


function DerEng(DPpt,G2)
	
	derengpt = 0.0
	
	for a in 1:3, b in 1:3
		derengpt += 0.5*G2[a,b]*DPpt[a]*DPpt[b]
	end
	
	return derengpt
	
end


function PotEng(Ppt,R2,G2,A2,A4,a111,a112,a123,V0)

	engpt = 0.0

	for a in 1:3, b in 1:3
		
		engpt += A2[a,b]*Ppt[a]*Ppt[b]
	
		for c in 1:3, d in 1:3
			engpt += A4[a,b,c,d]*Ppt[a]*Ppt[b]*Ppt[c]*Ppt[d]
		end
	end
	
	engpt += V6o(R2'*Ppt,a111,a112,a123)
	
	engpt -= V0
	
	return engpt
	
end


function engpt(Ppt,DPpt,R2,G2,A2,A4,a111,a112,a123,V0)
	
	return PotEng(Ppt,R2,G2,A2,A4,a111,a112,a123,V0) + DerEng(DPpt,G2)
	
end

function energy(P,ED,dx,N,R2,G2,A2,A4,a111,a112,a123,V0)
	
	Ppt = zeros(3)
	DPpt = zeros(3,3)
	
	for i in 3:N-2
		
		Ppt = getP(P,i)
		getDP!(DPpt,P,i,dx)
		
		#println(Ppt)
		#println(DPpt)

		ED[i] = engpt(Ppt,DPpt,R2,G2,A2,A4,a111,a112,a123,V0)
		
	end
	
	return sum(ED)*dx
	
end

function energyNOED(P,dx,N,R2,G2,A2,A4,a111,a112,a123,V0)
	
	Ppt = zeros(3)
	DPpt = zeros(3,3)

	engt = 0.0
	
	for i in 3:N-2
		
		Ppt = getP(P,i)
		getDP!(DPpt,P,i,dx)

		engt += engpt(Ppt,DPpt,R2,G2,A2,A4,a111,a112,a123,V0)
		
		#println(Ppt)
		#println(DPpt)

		
		
	end
	
	return engt*dx
	
end


function PotEng_tensor(Ppt,R2,G2,A2,A4,A6,V0)

	engpt = 0.0

	for a in 1:3, b in 1:3
		
		engpt += A2[a,b]*Ppt[a]*Ppt[b]
	
		for c in 1:3, d in 1:3
			engpt += A4[a,b,c,d]*Ppt[a]*Ppt[b]*Ppt[c]*Ppt[d]
		
			for e in 1:3, f in 1:3
				engpt += A6[a,b,c,d,e,f]*Ppt[a]*Ppt[b]*Ppt[c]*Ppt[d]*Ppt[e]*Ppt[f]

			end

		end


	end

	engpt -= V0
	
	return engpt
	
end


function engpt_tensor(Ppt,DPpt,R2,G2,A2,A4,A6,V0)
	
	return PotEng_tensor(Ppt,R2,G2,A2,A4,A6,V0) + DerEng(DPpt,G2)
	
end

function energy_tensor(P,ED,dx,N,R2,G2,A2,A4,A6,V0)
	
	Ppt = zeros(3)
	DPpt = zeros(3,3)
	
	for i in 3:N-2
		
		Ppt = getP(P,i)
		getDP!(DPpt,P,i,dx)
		
		ED[i] = engpt_tensor(Ppt,DPpt,R2,G2,A2,A4,A6,V0)
		
	end
	
	return sum(ED)*dx
	
end