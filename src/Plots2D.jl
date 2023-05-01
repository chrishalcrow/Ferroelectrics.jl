

function plot_strain_and_QPP(P0, ϵ0, Q3,Nx,Ny,x,y;scale=3.0)

    QPP = zeros(6,Nx,Ny)
    getQPP!(QPP,Q3,P0,Nx,Ny)

    fig = Figure(resolution = (6*300, 2*300).*scale)

    for a in 1:6


        heatmap( fig[1,a] , x, y, ϵ0[a,:,:] )
        heatmap( fig[2,a] , x, y, QPP[a,:,:] )

    end

  # end

    return fig

end



function plot_surfaces(G,n,m,x,y;scale=3.0)

    fig = Figure(resolution = (n*400, m*400).*scale)

    #for a in 1:n, b in 1:m

    ax = Axis3(fig[1,1])

        surface!( ax, x, y, G )

   # end

    return fig

end



function sortlengths(lengths)
    
    newlengths = fill(colorant"rgb(255,0,0)", (size(lengths)[1]));
    
    vacs = [ [1.0,1.0,1.0], [-1.0,1.0,1.0], [1.0,-1.0,1.0], [1.0,1.0,-1.0], [-1.0,-1.0,1.0], [-1.0,1.0,-1.0], [1.0,-1.0,-1.0], [-1.0,-1.0,-1.0] ]
    vaccols = [ colorant"rgb(255,255,255)",  colorant"rgb(0,255,255)",  colorant"rgb(255,0,255)" ,  colorant"rgb(255,255,0)" ,  colorant"rgb(0,0,255)" ,  colorant"rgb(0,255,0)" ,  colorant"rgb(255,0,0)",  colorant"rgb(0,0,0)"]

    vec = zeros(3)
    
    for i in 1:size(lengths)[1]  
        vec .= lengths[i]
        normalize!(vec);

        newlengths[i] = vaccols[argmin( [ norm(vec - vacs[i]) for i in 1:8 ] )]
    end
    
    return newlengths
    
end





function arrowplot(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq;scale=1.0)


    xlr = x[1:SFreq:Nx];
    ylr = y[1:SFreq:Ny];

    Plr0 = P0[:,1:SFreq:Nx,1:SFreq:Ny]
    Plr = 0.0.*Plr0

    for i in 1:size(Plr0)[2], j in 1:size(Plr0)[3]
        Plr[:,i,j] .= R2'*Plr0[:,i,j]
    end


    ps = [Point3f(xlr[i], ylr[j], 0.0) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];
    ns = [Point3f(Plr[1,i,j],Plr[2,i,j], Plr[3,i,j]) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];


    newlengths = sortlengths(ns)

    fig = Figure(resolution=(Ly/Lx*600,600).*scale, fontsize=12)


    axs = [Axis3(fig[1, i]; aspect=(Lx,Ly,min(Lx,Ly)), perspectiveness=0.0 , azimuth=0, elevation = pi/2,
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        zspinesvisible=false,
        zticksvisible=false,
        zticklabelsvisible=false,
        zlabelvisible=false
        
        ) for i=1:1]

    ascale = 1.5*scale

    arrows!(axs[1], ps, ns./3.0,  arrowsize=Vec3f(0.2*ascale, 0.2*ascale, 0.05*ascale), linewidth=0.02,  #arrowtail = Vec3f(1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0)),
        color= newlengths,
        align=:center)


    return fig

end


function arrow_and_contour_plot(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq;scale=1.0)


    xlr = x[1:SFreq:Nx];
    ylr = y[1:SFreq:Ny];

    Plr0 = P0[:,1:SFreq:Nx,1:SFreq:Ny]
    Plr = 0.0.*Plr0

    for i in 1:size(Plr0)[2], j in 1:size(Plr0)[3]
        Plr[:,i,j] .= R2'*Plr0[:,i,j]
    end


    ps = [Point3f(xlr[i], ylr[j], 0.0) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];
    ns = [Point3f(Plr[1,i,j],Plr[2,i,j], Plr[3,i,j]) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];


    newlengths = sortlengths(ns)

    fig = Figure(resolution=(2.0*Ly/Lx*600,600).*scale, fontsize=12)


    axs = Axis3(fig[1, 1]; aspect=(Lx,Ly,min(Lx,Ly)), perspectiveness=0.0 , azimuth=0, elevation = pi/2,
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        zspinesvisible=false,
        zticksvisible=false,
        zticklabelsvisible=false,
        zlabelvisible=false


        
        ) 

    ascale = 1.5*scale

    arrows!(axs, ps, ns./3.0,  arrowsize=Vec3f(0.2*ascale, 0.2*ascale, 0.05*ascale), linewidth=0.02,  #arrowtail = Vec3f(1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0)),
        color= newlengths,
        align=:center)


    ax2 = Axis(fig[1, 2];
ygridvisible = true,
        xgridvisible = false,
        yticklabelsvisible=false,
        xticksize=5

        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end

    contour!( ax2, y[3:Ny-2], -x[3:Nx-2], P0n[1,3:Nx-2,3:Ny-2]', levels=[0.0] )
    contour!( ax2, y[3:Ny-2], -x[3:Nx-2], P0n[2,3:Nx-2,3:Ny-2]', levels=[0.0] )
    contour!( ax2, y[3:Ny-2], -x[3:Nx-2], P0n[3,3:Nx-2,3:Ny-2]', levels=[0.0] )

    colgap!(fig.layout,1, -20)

    return fig

end






