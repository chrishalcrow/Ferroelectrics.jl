function dxB(phi, a, i, ls)
     @inbounds 0.5*(phi[a,i+1] - phi[a,i-1])/ls
end

function d2xB(phi, a, i, ls)
    @inbounds (phi[a,i+1] - 2.0*phi[a,i] + phi[a,i-1])/ls^2
end

function dxD(phi, a, i,ls)
    @inbounds (-phi[a,i+2] + 8.0*phi[a,i+1] - 8.0*phi[a,i-1] + phi[a,i-2] )/(12.0*ls)
end

function d2xD(phi, a,i,ls)
    @inbounds (-phi[a,i+2] + 16.0*phi[a,i+1] - 30.0*phi[a,i] + 16.0*phi[a,i-1] - phi[a,i-2] )/(12.0*ls^2)
end


function dxDsing(phi,  i,ls)
    @inbounds (-phi[i+2] + 8.0*phi[i+1] - 8.0*phi[i-1] + phi[i-2] )/(12.0*ls)
end

function d2xDsing(phi, i,ls)
    @inbounds (-phi[i+2] + 16.0*phi[i+1] - 30.0*phi[i] + 16.0*phi[i-1] - phi[i-2] )/(12.0*ls^2)
end






function dxD2D(phi, a, i, j, lsx)
    @inbounds (-phi[a,i+2,j] + 8.0*phi[a,i+1,j] - 8.0*phi[a,i-1,j] + phi[a,i-2,j])/(12.0*lsx)
end
function dyD2D(phi, a, i, j, lsy)
    @inbounds (-phi[a,i,j+2] + 8.0*phi[a,i,j+1] - 8.0*phi[a,i,j-1] + phi[a,i,j-2])/(12.0*lsy)
end

function d2xD2D(phi, a, i, j, lsx)
    @inbounds (-phi[a,i+2,j] + 16.0*phi[a,i+1,j] - 30.0*phi[a,i,j] + 16.0*phi[a,i-1,j] - phi[a,i-2,j])/(12.0*lsx^2)
end
function d2yD2D(phi, a, i, j, lsx)
    @inbounds (-phi[a,i,j+2] + 16.0*phi[a,i,j+1] - 30.0*phi[a,i,j] + 16.0*phi[a,i,j-1,] - phi[a,i,j-2])/(12.0*lsx^2)
end

function dxdydiffD2D(phi, a, i, j, lsx, lsy)
    @inbounds (-phi[a,i+2,j+2] + 16.0*phi[a,i+1,j+1] - 30.0*phi[a,i,j] + 16.0*phi[a,i-1,j-1] - phi[a,i-2,j-2])/(12.0*lsx*lsy)
end


function dxD2Du(phi, a, b, i, j, lsx)
    @inbounds (-phi[a][b,i+2,j] + 8.0*phi[a][b,i+1,j] - 8.0*phi[a][b,i-1,j] + phi[a][b,i-2,j])/(12.0*lsx)
end
function dyD2Du(phi, a, b, i, j, lsy)
    @inbounds (-phi[a][b,i,j+2] + 8.0*phi[a][b,i,j+1] - 8.0*phi[a][b,i,j-1] + phi[a][b,i,j-2])/(12.0*lsy)
end

function d2xD2Du(phi, a, b, i, j, lsx)
    @inbounds (-phi[a][b,i+2,j] + 16.0*phi[a][b,i+1,j] - 30.0*phi[a][b,i,j] + 16.0*phi[a][b,i-1,j] - phi[a][b,i-2,j])/(12.0*lsx^2)
end
function d2yD2Du(phi, a,b, i, j, lsx)
    @inbounds (-phi[a][b,i,j+2] + 16.0*phi[a][b,i,j+1] - 30.0*phi[a][b,i,j] + 16.0*phi[a][b,i,j-1,] - phi[a][b,i,j-2])/(12.0*lsx^2)
end

function dxdydiffD2Du(phi, a,b, i, j, lsx, lsy)
    @inbounds (-phi[a][b,i+2,j+2] + 16.0*phi[a][b,i+1,j+1] - 30.0*phi[a][b,i,j] + 16.0*phi[a][b,i-1,j-1] - phi[a][b,i-2,j-2])/(12.0*lsx*lsy)
end



function getP2D(P,i,j)
    return [ P[1,i,j], P[2,i,j], P[3,i,j] ]
end


function getDP2D!(DPpt,P,i,j,dx,dy)
    
    for a in 1:3
        DPpt[a,1] = dxD2D(P,a,i,j,dx)
        DPpt[a,2] = dyD2D(P,a,i,j,dy)
    end
    
end

function getDDP!(DPpt,P,i,j,dx,dy)
    
    for a in 1:3

        DPpt[a,1,1] = d2xD2D(P,a,i,j,dx)
        DPpt[a,2,2] = d2yD2D(P,a,i,j,dy)
        DPpt[a,1,2] = (dxdydiffD2D(P, a, i, j, dx, dy) - DPpt[a,1,1] - DPpt[a,2,2])/2.0
        DPpt[a,2,1] = DPpt[a,1,2]

    end
    
end

function getDDPeng!(DPpt,P,i,j,dx,dy)
    
    for a in 1:3

        DPpt[a,1,1] = d2xD2D(P,1,a,i,j,dx)
        DPpt[a,2,2] = d2yD2D(P,1,a,i,j,dy)
        DPpt[a,1,2] = (dxdydiffD2D(P, 1, a, i, j, dx, dy) - DPpt[a,1,1] - DPpt[a,2,2])/2.0
        DPpt[a,2,1] = DPpt[a,1,2]

    end
    
end







function dxDu(phi, a, b, i, j, lsx)
    @inbounds (-phi[a,b,i+2,j] + 8.0*phi[a,b,i+1,j] - 8.0*phi[a,b,i-1,j] + phi[a,b,i-2,j])/(12.0*lsx)
end
function dyDu(phi, a, b, i, j, lsy)
    @inbounds (-phi[a,b,i,j+2] + 8.0*phi[a,b,i,j+1] - 8.0*phi[a,b,i,j-1] + phi[a,b,i,j-2])/(12.0*lsy)
end

function d2xDu(phi, a, b, i, j, lsx)
    @inbounds (-phi[a,b,i+2,j] + 16.0*phi[a,b,i+1,j] - 30.0*phi[a,b,i,j] + 16.0*phi[a,b,i-1,j] - phi[a,b,i-2,j])/(12.0*lsx^2)
end
function d2yDu(phi, a, b, i, j, lsx)
    @inbounds (-phi[a,b,i,j+2] + 16.0*phi[a,b,i,j+1] - 30.0*phi[a,b,i,j] + 16.0*phi[a,b,i,j-1,] - phi[a,b,i,j-2])/(12.0*lsx^2)
end

function dxdydiffDu(phi, a, b, i, j, lsx, lsy)
    @inbounds (-phi[a,b,i+2,j+2] + 16.0*phi[a,b,i+1,j+1] - 30.0*phi[a,b,i,j] + 16.0*phi[a,b,i-1,j-1] - phi[a,b,i-2,j-2])/(12.0*lsx*lsy)
end



function getP(P,i)
    return [ P[1,i], P[2,i], P[3,i] ]
end

function getDP!(DPpt,P,i,dx)
    
    for a in 1:3
        DPpt[a] = dxD(P,a,i,dx)
    end
    
end

function getDDP!(DDPpt,P,i,dx)
    
    for a in 1:3
        DDPpt[a] = d2xD(P,a,i,dx)
    end
    
end


function getQPP!(QPP,Q3,P0,Nx,Ny)

    for α in 1:6, i in 1:Nx, j in 1:Ny

        QPP[α,i,j] = 0.0

        for b in 1:3, c in 1:3
            QPP[α,i,j] += Q3[α,b,c]*P0[b,i,j]*P0[c,i,j]
        end

    end

end

function getQPP(Q3,P0,Nx,Ny)

    QPP = zeros(6,Nx,Ny)

    for α in 1:6, b in 1:3, c in 1:3, i in 1:Nx, j in 1:Ny
        #if α == 3
          #  QPP[α,i,j] = Q3[α+3,b,c]P0[b,i,j]*P0[c,i,j]
        #else
            QPP[α,i,j] += Q3[α,b,c]*P0[b,i,j]*P0[c,i,j]
        #end
    end

    return QPP

end

