function getH(P,N,dx,G2,A2,A4,a111,a112,a123,R2)


    Dxx_small = zeros(N,N)

    for i in 1:N
        if i-2 > 0
            Dxx_small[i-2,i] = -1.0/12.0
        end
        if i-1 > 0
            Dxx_small[i-1,i] = 4.0/3.0
        end
        Dxx_small[i,i] = -5.0/2.0
        if i+1 < N+1
            Dxx_small[i+1,i] = 4.0/3.0
        end
        if i+2 < N+1
            Dxx_small[i+2,i] = -1.0/12.0
        end
    end

    Dxx_small ./= dx^2;

    Dxx = zeros(3*N,3*N)

    Dxx = -Symmetric( [ G2[1,1]*Dxx_small G2[1,2]*Dxx_small G2[1,3]*Dxx_small ;
        G2[2,1]*Dxx_small G2[2,2]*Dxx_small G2[2,3]*Dxx_small ;
        G2[3,1]*Dxx_small G2[3,2]*Dxx_small G2[3,3]*Dxx_small ] );


    D2VP = zeros(3,3,N,N)

    D2V6o = zeros(3,3)
    D2V6 = zeros(3,3)

    A2ss = zeros(3,3)
    A4ss = zeros(3,3,3,3)

    for a in 1:3, b in 1:3
        A2ss[a,b] = 0.5*(A2[a,b] + A2[b,a])
        for c in 1:3, d in 1:3
            A4ss[a,b,c,d] = (A4[a,b,c,d] + A4[a,c,b,d] + A4[a,c,d,b]
                          + A4[b,a,c,d] + A4[c,a,b,d] + A4[c,a,d,b]
                          + A4[b,c,a,d] + A4[c,b,a,d] + A4[c,d,a,b]
                          + A4[b,c,d,a] + A4[c,b,d,a] + A4[c,d,b,a] )/12.0
        end
    end

    for i in 1:N
        
        for a in 1:3, b in 1:3    
        
            D2VP[a,b,i,i] += 2.0*A2ss[a,b]
        
            for c in 1:3, d in 1:3
            
                D2VP[a,b,i,i] += 12.0*A4ss[a,b,c,d]*P[c,i]*P[d,i]
            
            end
        
        end
        
        D2V6o[1,1] = d11V6(R2'*P[:,i], a111,a112,a123)
        D2V6o[1,2] = d12V6(R2'*P[:,i], a111,a112,a123)
        D2V6o[1,3] = d13V6(R2'*P[:,i], a111,a112,a123)
        
        D2V6o[2,1] = d21V6(R2'*P[:,i], a111,a112,a123)
        D2V6o[2,2] = d22V6(R2'*P[:,i], a111,a112,a123)
        D2V6o[2,3] = d23V6(R2'*P[:,i], a111,a112,a123)
        
        D2V6o[3,1] = d31V6(R2'*P[:,i], a111,a112,a123)
        D2V6o[3,2] = d32V6(R2'*P[:,i], a111,a112,a123)
        D2V6o[3,3] = d33V6(R2'*P[:,i], a111,a112,a123)
        
        for a in 1:3, b in 1:3
            
            D2V6[a,b] = 0.0
            
            for c in 1:3, d in 1:3
               
                D2V6[a,b] += R2[a,c]*R2[b,d]*D2V6o[c,d]
                
            end
            
            D2VP[a,b,i,i] += D2V6[a,b]
            
        end
        
    end

    Dxx +=  [  D2VP[1,1,:,:] D2VP[1,2,:,:] D2VP[1,3,:,:]  ;
        D2VP[2,1,:,:] D2VP[2,2,:,:] D2VP[2,3,:,:]  ;
        D2VP[3,1,:,:] D2VP[3,2,:,:] D2VP[3,3,:,:]   ]


    return Dxx

end


function getfreqs(H)

    return eigvals(Symmetric(H))

end

function getvecs(H)

    return eigvecs(Symmetric(H))

end