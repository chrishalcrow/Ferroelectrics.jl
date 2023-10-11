

function dx(phi, a, i)
    @inbounds (-phi.field[a,i+2] + 8.0*phi.field[a,i+1] - 8.0*phi.field[a,i-1] + phi.field[a,i-2] )/(12.0*phi.grid.ls)
end

function d2x(phi, a,i)
    @inbounds (-phi.field[a,i+2] + 16.0*phi.field[a,i+1] - 30.0*phi.field[a,i] + 16.0*phi.field[a,i-1] - phi.field[a,i-2] )/(12.0*phi.grid.ls^2)
end




function getP(P,i)
    return SVector{3,Float64}(
        P.field[1,i],
        P.field[2,i],
        P.field[3,i]
         )
end

function getDP(P,i)
    
    return SVector{3,Float64}(
        dx(P,1,i),
        dx(P,2,i),
        dx(P,3,i)
    )
        
end

function getDDP(P,i)
    
    return SVector{3,Float64}(
        d2x(P,1,i),
        d2x(P,2,i),
        d2x(P,3,i)
    )
        
end

