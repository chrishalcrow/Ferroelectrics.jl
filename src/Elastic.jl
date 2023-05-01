function get_ϵbestL_no_compat!(ϵbest,QPP,Nx,Ny)

    for a in 1:6, i in 1:Nx, j in 1:Ny
        ϵbest[a,i,j] = QPP[a,i,j]
    end

end

function check_compat(ϵ,Nx,Ny,dx,dy)

    howcompat = zeros(Nx,Ny)

    for i in 3:Nx-2, j in 3:Ny-2
        howcompat[i,j] = check_compat_pt(ϵ,i,j,dx,dy)
    end

    return howcompat

end



function check_compat_pt(ϵ,i,j,dx,dy)

        return d2yD2D(ϵ,1,i,j,dy)  - 2.0*(dxdydiffD2D(ϵ, 6, i, j, dx, dy) - d2xD2D(ϵ,6,i,j,dx) - d2yD2D(ϵ,6,i,j,dy))/2.0 + d2xD2D(ϵ,2,i,j,dx) 
end



function get_ϵbestL!(ϵbest,ϵbasis,QPP,cc00,CL,ncsb,Nx,Ny,dx,dy,Lx,Ly,pmax,qmax,Cabm)

    qPPCeL = zeros(8)

    @inbounds for a in 1:6, i in 1:Nx, j in 1:Ny
    	ϵbest[a,i,j] = 0.0
    end

    @inbounds for p in 0:pmax, q in 0:qmax

        makeϵbf!(ϵbasis,p,q,cc00,CL,ncsb,Nx,Ny,Lx,Ly)
        DOT_C_L!(qPPCeL, QPP,ϵbasis,Cabm,Nx,Ny,dx,dy)

        for i in 1:Nx, j in 1:Ny

			ϵbest[3,i,j] = QPP[3,i,j]
	    	ϵbest[4,i,j] = QPP[4,i,j]
	    	ϵbest[5,i,j] = QPP[5,i,j]

	        for a in 1:8

	        	ϵbest[1,i,j] += qPPCeL[a]*ϵbasis[a,1,i,j]
	        	ϵbest[2,i,j] += qPPCeL[a]*ϵbasis[a,2,i,j]
	        	ϵbest[6,i,j] += qPPCeL[a]*ϵbasis[a,3,i,j]

	        end


        end

    end
    
end

function DOT_C(ϵ1, ϵ2, Cabm, Nx, Ny,dx,dy)
   
    eCe = 0.0
    
    for a in 1:3, b in 1:3, i in 3:Nx-2, j in 3:Ny-2
        
        eCe += ϵ1[a,i,j]*Cabm[a,b]*ϵ2[b,i,j]
        
    end
    
    return eCe*dx*dy
    
end

function DOT_C_2(ϵ1, ϵL, c, Cabm, Nx, Ny,dx,dy)
   
    eCe = 0.0
    
    @inbounds for a in 1:3, b in 1:3, i in 3:Nx-2, j in 3:Ny-2
        
        eCe += ϵ1[a,i,j]*Cabm[a,b]*ϵL[c,b,i,j]
        
    end
    
    return eCe*dx*dy
    
end

function threetosix(a)
	if a == 3
		return 6
	else
		return a
	end
end

function DOT_C_L!(qPPCe,ϵ1,ϵL,Cabm,Nx,Ny,dx,dy)
   
    @inbounds for c in 1:8
        
        qPPCe[c] = 0.0
    
        @inbounds for a in 1:3, b in 1:3, i in 3:Nx-2, j in 3:Ny-2
        
            qPPCe[c] += ϵ1[threetosix(a),i,j]*Cabm[a,b]*ϵL[c,b,i,j]
        
        end
        
        qPPCe[c] *= dx*dy
    
    end
    
end



function DOT_NO_C(ϵ1, ϵ2, Nx, Ny,dx,dy)
   
    eCe = 0.0
    
    for i in 3:Nx-2, j in 3:Ny-2
        
        eCe += ϵ1[i,j]*ϵ2[i,j]
        
    end
    
    return eCe*dx*dy
    
end


function getALLcs(pmax,qmax,Cabm,Lx,Ly)

	CL = zeros(pmax,qmax,8,2,3)

	for p in 1:pmax
	    for q in 1:qmax

	    	pp = 2.0*p*pi/Lx	    	
	    	qq = 2.0*pi*q/Ly


			e1 = [ pp 0.0 0.0 ; 0.0 0.0 -0.5*qq ]
			e2 = [ 0.0 qq 0.0 ; 0.0 0.0 -0.5*pp ]
			e3 = [ 0.0 0.0 -0.5*qq ; pp 0.0 0.0 ]
			e4 = [ 0.0 0.0 -0.5*pp ; 0.0 qq 0.0 ]

			e5 = [ pp 0.0 0.0 ; 0.0 0.0 0.5*qq ]
			e6 = [ 0.0 qq 0.0 ; 0.0 0.0 0.5*pp ]
			e7 = [ 0.0 0.0 0.5*qq ; pp 0.0 0.0 ]
			e8 = [ 0.0 0.0 0.5*pp ; 0.0 qq 0.0 ]

			cc1,cc2,cc3,cc4=GramSchmidtC(e1,e2,e3,e4,Cabm)
			cs1,cs2,cs3,cs4=GramSchmidtC(e5,e6,e7,e8,Cabm)

			CL[p,q,1,:,:] = cc1
			CL[p,q,2,:,:] = cc2
			CL[p,q,3,:,:] = cc3
			CL[p,q,4,:,:] = cc4

			CL[p,q,5,:,:] = cs1
			CL[p,q,6,:,:] = cs2
			CL[p,q,7,:,:] = cs3
			CL[p,q,8,:,:] = cs4

		end

	end

	return CL

end



function getnbasis(pmax,qmax,x,y,Nx,Ny,Lx,Ly)

	ncosx = zeros(pmax, Nx)
	nsinx = zeros(pmax, Nx)

	ncosy = zeros(qmax, Ny)
	nsiny = zeros(qmax, Ny)

	for p in 1:pmax
	    pp = 2.0*p*pi/Lx
	    for i in 1:Nx
	        ncosx[p,i] = cos(pp*x[i])*sqrt(2.0/Lx)
	        nsinx[p,i] = sin(pp*x[i])*sqrt(2.0/Lx)
	    end
	end

	for q in 1:qmax
	    qq = 2.0*pi*q/Ly
	    for j in 1:Ny
	        ncosy[q,j] = cos(qq*y[j])*sqrt(2.0/Ly)
	        nsiny[q,j] = sin(qq*y[j])*sqrt(2.0/Ly)
	    end
	end

	return ncosx, ncosy, nsinx, nsiny

end





function makeϵbf!(ϵb,p,q,cc00,cs,ncsxys,Nx,Ny,Lx,Ly)

    #ϵb1 = zeros(3,Nx,Ny)

    ncosx, ncosy, nsinx, nsiny = ncsxys
    
    if p == 0 && q == 0
        
        #println("p and q are zero")
        
        sqLxLy = sqrt(Lx*Ly)
        
    @inbounds    for a in 1:3, i in 1:Nx, j in 1:Ny
            for b in 1:3
                ϵb[b,a,i,j] = cc00[b][a]/sqLxLy
            end
            for b in 4:8
                ϵb[b,a,i,j] = 0.0
            end
        end
    elseif p == 0
        
        #println("only p is zero")
        
        sqLx = sqrt(Lx)
        
     @inbounds   for a in 1:3, i in 1:Nx, j in 1:Ny
            for b in 1:3
                ϵb[b,a,i,j] = cc00[b][a]*ncosy[q,j]/sqLx
            end
            for b in 4:6
                ϵb[b,a,i,j] = cc00[b-3][a]*nsiny[q,j]/sqLx
            end
            for b in 7:8
                ϵb[b,a,i,j] = 0.0
            end
        end
    elseif q == 0
        
        #println("only q is zero")
        
        sqLy = sqrt(Ly)
        
       @inbounds for a in 1:3, i in 1:Nx, j in 1:Ny
            for b in 1:3
                ϵb[b,a,i,j] = cc00[b][a]*ncosx[p,i]/sqLy
            end
            for b in 4:6
                ϵb[b,a,i,j] = cc00[b-3][a]*nsinx[p,i]/sqLy
            end
            for b in 7:8
                ϵb[b,a,i,j] = 0.0
            end
        end
    else

       # println("else!")

        @inbounds for a in 1:3, i in 1:Nx, j in 1:Ny

            @fastmath ϵb[1,a,i,j] = cs[p,q,1,1,a]*ncosx[p,i]*ncosy[q,j] + cs[p,q,1,2,a]*nsinx[p,i]*nsiny[q,j]
            @fastmath ϵb[2,a,i,j] = cs[p,q,2,1,a]*ncosx[p,i]*ncosy[q,j] + cs[p,q,2,2,a]*nsinx[p,i]*nsiny[q,j]
            @fastmath ϵb[3,a,i,j] = cs[p,q,3,1,a]*ncosx[p,i]*ncosy[q,j] + cs[p,q,3,2,a]*nsinx[p,i]*nsiny[q,j]
            @fastmath ϵb[4,a,i,j] = cs[p,q,4,1,a]*ncosx[p,i]*ncosy[q,j] + cs[p,q,4,2,a]*nsinx[p,i]*nsiny[q,j]

            @fastmath ϵb[5,a,i,j] = cs[p,q,5,1,a]*ncosx[p,i]*nsiny[q,j] + cs[p,q,5,2,a]*nsinx[p,i]*ncosy[q,j]
            @fastmath ϵb[6,a,i,j] = cs[p,q,6,1,a]*ncosx[p,i]*nsiny[q,j] + cs[p,q,6,2,a]*nsinx[p,i]*ncosy[q,j]
            @fastmath ϵb[7,a,i,j] = cs[p,q,7,1,a]*ncosx[p,i]*nsiny[q,j] + cs[p,q,7,2,a]*nsinx[p,i]*ncosy[q,j]
            @fastmath ϵb[8,a,i,j] = cs[p,q,8,1,a]*ncosx[p,i]*nsiny[q,j] + cs[p,q,8,2,a]*nsinx[p,i]*ncosy[q,j]

        end
        
   end

end



function get_C2m(C2)

	Cabm = zeros(3,3)

	for a in 1:3, b in 1:3
	    if a == 3 && b == 3
	        Cabm[a,b] = C2[a+3, b+3]
	    elseif a == 3
	        Cabm[a,b] = C2[a+3, b]
	    elseif b == 3
	        Cabm[a,b] = C2[a, b+3]
	    else
	        Cabm[a,b] = C2[a,b]
	    end
	end

	return Cabm

end


function getcs(pp,qq,Cabm)


	e1 = [ pp 0.0 0.0 ; 0.0 0.0 -0.5*qq ]
	e2 = [ 0.0 qq 0.0 ; 0.0 0.0 -0.5*pp ]
	e3 = [ 0.0 0.0 -0.5*qq ; pp 0.0 0.0 ]
	e4 = [ 0.0 0.0 -0.5*pp ; 0.0 qq 0.0 ]

	e5 = [ pp 0.0 0.0 ; 0.0 0.0 0.5*qq ]
	e6 = [ 0.0 qq 0.0 ; 0.0 0.0 0.5*pp ]
	e7 = [ 0.0 0.0 0.5*qq ; pp 0.0 0.0 ]
	e8 = [ 0.0 0.0 0.5*pp ; 0.0 qq 0.0 ]

	cc1,cc2,cc3,cc4=GramSchmidtC(e1,e2,e3,e4,Cabm)
	cs1,cs2,cs3,cs4=GramSchmidtC(e5,e6,e7,e8,Cabm)

	return cc1, cc2, cc3, cc4, cs1, cs2, cs3, cs4

end


function getc00(Cabm)

	cc00 = GramSchmidt00([1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0] ,Cabm)

	return cc00

end

function dotC(e1, e2, Cabm)

    eCe = 0.0
    
    for a in 1:3, b in 1:3
        eCe += e1[1,a]*Cabm[a,b]*e2[1,b] + e1[2,a]*Cabm[a,b]*e2[2,b]
    end
    
    return eCe
end

function dotC00(e1, e2, Cabm)
    
    eCe = 0.0
    
    for a in 1:3, b in 1:3
        eCe += e1[a]*Cabm[a,b]*e2[b]
    end
    
    return eCe
end

function normalizeC!(f1,Cabm)
    
    f1 ./= sqrt(dotC(f1,f1,Cabm))
        
end

function normalizeC00!(f1,Cabm)
    
    f1 ./= sqrt(dotC00(f1,f1,Cabm))
        
end

function GramSchmidtC(e1,e2,e3,e4,Cabm)
    
    f1 = e1
    normalizeC!(f1,Cabm)
    
    f2 = e2 - dotC(f1,e2,Cabm)*f1
    normalizeC!(f2,Cabm)
    
    f3 = e3 - dotC(f1,e3,Cabm)*f1 - dotC(f2,e3,Cabm)*f2
    normalizeC!(f3,Cabm)
    
    f4 = e4 - dotC(f1,e4,Cabm)*f1 - dotC(f2,e4,Cabm)*f2 - dotC(f3,e4,Cabm)*f3
    normalizeC!(f4,Cabm)
    
    return f1 , f2 , f3 , f4
    
end

function GramSchmidt00(e1,e2,e3,Cabm)
    
    f1 = e1
    normalizeC00!(f1,Cabm)
    
    f2 = e2 - dotC00(f1,e2,Cabm)*f1
    normalizeC00!(f2,Cabm)
    
    f3 = e3 - dotC00(f1,e3,Cabm)*f1 - dotC00(f2,e3,Cabm)*f2
    normalizeC00!(f3,Cabm)
    
    
    return f1, f2, f3
    
end








