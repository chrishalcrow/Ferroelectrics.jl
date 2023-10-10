



function FastCyclicConv(f,g,n) 
	return (n > 1 ? ifft(fft(f,2).*fft(g,2),2) : ifft(fft(f).*fft(g)))
end

function FastLinearConvolution(f,g,n)
    if n == 1 && size(f)[1] == length(f) && size(g)[1] == length(g)
        N = length(f)
        M = length(g)

        f_pad = [ f; zeros(M-1) ]     
        g_pad = [ g; zeros(N-1) ]
    else
        N = size(f)[2]
        M = size(g)[2]

        f_pad = [ f[1:n,:] zeros(n, M-1) ]     
        g_pad = [ g[1:n,:] zeros(n, N-1) ]     
    end
    return FastCyclicConv(f_pad, g_pad,n)
end

function DipoleEnergy(P, sums, type, shift=2)
	dipole_energy = 0.0
	N = P.lp

	if type == 1 # Approach #1: capacitor model. P.pars[10] should be abs(g1), energy calculation for the whole sample
		dipole_energy = sum(P.field[1,(1+shift):(N-shift)].^2)*P.ls/(2*8.8541878128e-12*P.parameters.rescaling)
	elseif type == 2 #Approach #2: interacting planes infinite in y direction and width b (P.b in code) in z direction. P.pars[10] should be abs(g1)
		for i in (1+shift):(N-shift)
			for j in i:(N-shift)
				dipole_energy += P.field[1,i]*P.field[1,j]*sums[1,abs(i-j)+1] + P.field[3,i]*P.field[3,j]*sums[2,abs(i-j)+1]
			end
		end
		dipole_energy *= P.b/(4*pi*8.8541878128e-12*P.parameters.rescaling)
	elseif type == 3
		#alternative way using FFT and FastLinearConvolution
		#dipole_energy = real(sum(P.field[1,:]		 .*(FastLinearConvolution(P.field,sums,1)[1:N]))) 					 #only component P[1], energy calculation for the whole sample
		#dipole_energy = real(sum(P.field[1,3:N-2]	 .*(FastLinearConvolution(P.field[1,3:N-2],sums[1,:],1)[1:N-4]))) 		 #component P[1], energy calculation for sites 3:P.lp-2
		#dipole_energy = real(sum(P.field[1:2:3,:]	 .*(FastLinearConvolution(P.field[1:2:3,:],sums,2)[:,1:N]))) 		 #components P[1], P[3], energy calculation for the whole sample
		
		dipole_energy = real(sum(P.field[1:2:3,(1+shift):(N-shift)].*(FastLinearConvolution(P.field[1:2:3,(1+shift):(N-shift)],sums[:,1:(N-2*shift)],2)[:,1:(N-2*shift)])))   #components P[1], P[3], energy calculation for sites 3:P.lp-2

		dipole_energy *= P.b/(4*pi*8.8541878128e-12*P.parameters.rescaling)
	end
	
	return dipole_energy
end

function energy(P, sums, type=1, shift=2)
	
	ED = zeros(P.lp)
	
	Threads.@threads for i in (1+shift):(P.lp-shift)

		Ppt = getP(P,i)
		DPpt = getDP(P,i)

		ED[i] = engpt(Ppt,DPpt, P.parameters)
	end

	total_energy = sum(ED)*P.ls

	if P.is_electrostatic
		total_energy += DipoleEnergy(P, sums, type, shift)
	end
	

	return total_energy

	
end

function changeAs(P,sums,type,shift=2)
	if P.is_electrostatic
		N = P.lp
	
		P.field = zeros(3,N)
		P.field[3,:] .= 1.0
		A2 = DipoleEnergy(P, sums, type, shift)/((P.lp-2*shift)*P.ls)
		P.parameters.A2t[3,3] -= A2   #correction to A2t
		P.parameters.A2s[3,3] -= 2*A2 #correction to A2s

		P.field = zeros(3,N)
		P.field[1,:] .= 1.0
		A2 = DipoleEnergy(P, sums, type, shift)/((P.lp-2*shift)*P.ls)
		P.parameters.A2t[1,1] -= A2   #correction to A2t
		P.parameters.A2s[1,1] -= 2*A2 #correction to A2s
	end
end

function EnergyANF(P, ED)
	
	Threads.@threads for i in 3:P.lp-2

		Ppt = getP(P,i)
		DPpt = getDP(P,i)

		ED[i] = engpt(Ppt,DPpt, P.parameters)
		
	end
	
	return sum(ED)*P.ls
	
end

function engpt(Ppt,DPpt,pars)
	
	return PotEng(Ppt, pars.A2t, pars.A4t, pars.V0) + DerEng(DPpt,pars.G2)
	
end


function DerEng(DPpt,G2)
	
	derengpt = 0.0
	
	for a in 1:3, b in 1:3
		derengpt += 0.5*G2[a,b]*DPpt[a]*DPpt[b]
	end
	
	return derengpt
	
end


function PotEng(Ppt, A2t, A4t, V0)

	engpt = 0.0

	for a in 1:3, b in 1:3
		
		engpt += A2t[a,b]*Ppt[a]*Ppt[b]
	
		for c in 1:3, d in 1:3
			engpt += A4t[a,b,c,d]*Ppt[a]*Ppt[b]*Ppt[c]*Ppt[d]
		end
	end

	engpt -= V0
	
	return engpt
	
end


function DWwidth(P)
	for i in 2:P.lp
		if P.field[3,i]*P.field[3,i-1] < 0
			print("DW thickness ", round(1e9* P.parameters.length_scale *P.ls * 2. / abs(P.field[3,i - 1] - P.field[3,i]),digits=5)," nm")
			break
		end
	end
end