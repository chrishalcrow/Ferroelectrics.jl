

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


function arrow_and_contour_plot(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq;scale=1.0,savefile=false,filename="temp.png")


    xlr = x[3:SFreq:Nx-2];
    ylr = y[3:SFreq:Ny-2];

    Plr0 = P0[:,3:SFreq:Nx-2,3:SFreq:Ny-2]
    Plr = 0.0.*Plr0

    for i in 1:size(Plr0)[2], j in 1:size(Plr0)[3]
        Plr[:,i,j] .= R2'*Plr0[:,i,j]
    end


    ps = [Point3f(xlr[i], ylr[j], 0.0) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];
    ns = [Point3f(Plr[1,i,j],Plr[2,i,j], Plr[3,i,j]) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];


    newlengths = sortlengths(ns)

    fig = Figure(resolution=(2.0*Lx/Ly*600,600).*scale, fontsize=12)


    axs = Axis3(fig[1, 1]; aspect=(Lx,Ly,1.0*min(Lx,Ly)), perspectiveness=0.0 , azimuth=-pi/2-0.00001, elevation = pi/2,
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        zspinesvisible=false,
        zticksvisible=false,
        zticklabelsvisible=false,
        zlabelvisible=false,
        yautolimitmargin = (0.02,0.02),
        xlabelsize = 12*scale,
        ylabelsize = 12*scale,
        ylabel=L"r",
        xlabel=L"s",
        tellwidth=true,tellheight=true,
        xticklabelsize = 10*scale,
        yticklabelsize = 10*scale





        
        ) 

    ascale = 1.0*scale

    arrows!(axs, ps, ns./10.0,  arrowsize=Vec3f(0.2*ascale, 0.2*ascale, 0.2), linewidth=0.1,  #arrowtail = Vec3f(1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0)),
        color= newlengths,
        align=:center)

    zlims!(-1.5,1.5)



    xcut=0

    ax2 = Axis3(fig[1, 2]; aspect=(Lx,Ly,1.0*min(Lx,Ly)), 
    perspectiveness=0.0 , azimuth=-pi/2-0.00001, elevation = pi/2,
    xgridvisible = false,
    #ygridvisible = false,
    zgridvisible = false,
    zspinesvisible=false,
    zticksvisible=false,
    zticklabelsvisible=false,
    zlabelvisible=false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    yticklabelsvisible=false,
    ylabelvisible=false,
    xlabelsize = 12*scale,
    xlabel=L"s",
    xticklabelsize = 10*scale


        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    
    zlims!(-1.5,1.5)



    ycut=0

    theLW = 5

   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue )



    colgap!(fig.layout,1, -61)

    if savefile == true
        imusingnotebook()
        save(filename, fig)
        #imusingterminal()
    end



    return fig

end





function arrow_and_contour_plot2(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq;scale=1.0,savefile=false,filename="temp.png")


    xlr = x[3:SFreq:Nx-2];
    ylr = y[3:SFreq:Ny-2];

    Plr0 = P0[:,3:SFreq:Nx-2,3:SFreq:Ny-2]
    Plr = 0.0.*Plr0

    for i in 1:size(Plr0)[2], j in 1:size(Plr0)[3]
        Plr[:,i,j] .= R2'*Plr0[:,i,j]
    end


    ps = [Point3f(xlr[i], ylr[j], 0.0) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];
    ns = [Point3f(Plr[1,i,j],Plr[2,i,j], Plr[3,i,j]) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];


    newlengths = sortlengths(ns)

    fig = Figure(resolution=(2*Lx/Ly*600,600).*scale, fontsize=12)

    

    axs = Axis3(fig[1, 1]; aspect=(Lx,Ly,1.0*min(Lx,Ly)),
    perspectiveness=0.0 , azimuth=-pi/2-0.00001, elevation = pi/2,
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        zspinesvisible=false,
        zticksvisible=false,
        zticklabelsvisible=false,
        zlabelvisible=false,
        yautolimitmargin = (0.02,0.02),
        xlabelsize = 24*scale,
        ylabelsize = 24*scale,
        ylabel=L"r",
        xlabel=L"s",
        ylabeloffset = 20,
        xlabeloffset = 20,
        tellwidth=true,tellheight=true,
        xticklabelsize = 20*scale,
        yticklabelsize = 20*scale





        
        ) 

    ascale = 1.0*scale

    arrows!(axs, ps, ns./10.0,  arrowsize=Vec3f(0.2*ascale, 0.2*ascale, 0.2), linewidth=0.1,  #arrowtail = Vec3f(1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0)),
        color= newlengths,
        align=:center)

    zlims!(-1.5,1.5)

    #Box(fig[1,1],   color = (:blue,0.2) )

    xcut=0

    ax2 = Axis3(fig[1, 2]; aspect=(Lx,Ly,1.0*min(Lx,Ly)), 
    perspectiveness=0.0 , azimuth=-pi/2-0.00001, elevation = pi/2,
    xgridvisible = false,
    #ygridvisible = false,
    zgridvisible = false,
    zspinesvisible=false,
    zticksvisible=false,
    zticklabelsvisible=false,
    zlabelvisible=false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    xlabeloffset = 20,
    yticklabelsvisible=false,
    ylabelvisible=false,
    xlabelsize = 24*scale,
    xlabel=L"s",
    xticklabelsize = 20*scale


        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    
    zlims!(-1.5,1.5)



    ycut=0

    theLW = 5

   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red, label="P1")
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green, label="P2")
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue , label="P3")

   #Box(fig[1,2],   color = (:red,0.2) )



   colsize!(fig.layout, 1, Auto(1.0))
    colgap!(fig.layout,1, -62)

    if savefile == true
        imusingnotebook()
        save(filename, fig)
        #imusingterminal()
    end

    resize_to_layout!(fig)

    return fig

end







function arrow_and_contour_plot2D(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq;scale=1.0,savefile=false,filename="temp.png")


    xlr = x[3:SFreq:Nx-2];
    ylr = y[3:SFreq:Ny-2];

    Plr0 = P0[:,3:SFreq:Nx-2,3:SFreq:Ny-2]
    Plr = 0.0.*Plr0

    for i in 1:size(Plr0)[2], j in 1:size(Plr0)[3]
        Plr[:,i,j] .= R2'*Plr0[:,i,j]
    end


    ps = [Point3f(xlr[i], ylr[j], 0.0) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];
    ns = [Point3f(Plr[1,i,j],Plr[2,i,j], Plr[3,i,j]) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];


    newlengths = sortlengths(ns)

    fig = Figure(resolution=(2*Lx/Ly*600,600).*scale, fontsize=12)

    

    axs = Axis3(fig[1, 1]; aspect=(Lx,Ly,1.0*min(Lx,Ly)),
    perspectiveness=0.0 , azimuth=-pi/2-0.00001, elevation = pi/2,
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        zspinesvisible=false,
        zticksvisible=false,
        zticklabelsvisible=false,
        zlabelvisible=false,
        yautolimitmargin = (0.02,0.02),
        xlabelsize = 12*scale,
        ylabelsize = 12*scale,
        ylabel="r",
        xlabel="s",
        tellwidth=true,tellheight=true,
        xticklabelsize = 10*scale,
        yticklabelsize = 10*scale





        
        ) 

    ascale = 1.0*scale

    arrows!(axs, ps, ns./10.0,  arrowsize=Vec3f(0.2*ascale, 0.2*ascale, 0.2), linewidth=0.1,  #arrowtail = Vec3f(1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0)),
        color= newlengths,
        align=:center)

    zlims!(-1.5,1.5)

    Box(fig[1,1],   color = (:blue,0.2) )

    xcut=0

    ax2 = Axis(fig[1, 2]; aspect=Lx/Ly,
    xgridvisible = false,
    #ygridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    yticklabelsvisible=false,
    ylabelvisible=false,
    xlabelsize = 12*scale,
    xlabel=L"s",
    xticklabelsize = 10*scale


        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    

    ycut=0

    theLW = 5

   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue )

   Box(fig[1,2],   color = (:red,0.2) )

   #for (i, label) in enumerate(["sleep", "awake", "test"])
   #     Box(fig[2, 2], color = :gray90)
   #     Label(fig[2, 2], label, rotation = 0.0, tellheight = false)
   # end

   colsize!(fig.layout, 1, Auto(1.0))
   # colgap!(fig.layout,1, -62)

    if savefile == true
        imusingnotebook()
        save(filename, fig)
        #imusingterminal()
    end

    resize_to_layout!(fig)

    return fig

end











function testfn()

    println("hello")

end


function SLf1()


    #fig=Figure(resolution=(2*Lx/Ly*600,600).*scale, fontsize=12)

    fig=Figure

    img = load("line.png");

    ImX, ImY = size(img)

    image(fig[1, 2], img',

        axis = (aspect = DataAspect(),
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        xticksvisible=false,
        yticksvisible=false,
        yreversed=true,
    ))

    innerradius = ImX*0.098
    outerradius = ImX*0.433

    howthick = 6

    vlines!([34], ymin = 0.05, ymax = 0.95, linewidth=howthick)
    vlines!([474], ymin = 0.05, ymax = 0.95,  linewidth=howthick)

    hlines!([22], xmin = 0.06, xmax = 0.94,  color = :red, linewidth=howthick)
    hlines!([483], xmin = 0.06, xmax = 0.94, color = :green, linewidth=howthick)

    themarkersize=25

    scatter!(34,475, marker = :utriangle, markersize=themarkersize )
    scatter!(40,22, marker = :ltriangle, markersize=themarkersize, color=:red )
    scatter!(474,32, marker = :dtriangle, markersize=themarkersize, color=:orange )
    scatter!(465,484, marker = :rtriangle, markersize=themarkersize, color=:green )


    img = load("sky.png");
    size(img)

    ImX, ImY = size(img)

    image(fig[1,1], img',

        axis = (aspect = DataAspect(),
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        xticksvisible=false,
        yticksvisible=false,
        yreversed=true
    ))

    innerradius = ImX*0.098
    outerradius = ImX*0.433

    howthick = 6

    arc!(Point2f(ImX*0.5), innerradius, -3π/2+0.2, π/2-0.2,linewidth=howthick)
    arc!(Point2f(ImX*0.5), outerradius, -3π/2+0.05, π/2-0.05,linewidth=howthick)

    vlines!([243], ymin = 0.595, ymax = 0.93, color = :red, linewidth=howthick)
    vlines!([265], ymin = 0.595, ymax = 0.93, color = :green, linewidth=howthick)

    scatter!(270,302, marker = :dtriangle, markersize=themarkersize )
    scatter!(243,313, marker = :dtriangle,  markersize=themarkersize, color=:red )
    scatter!(234,475, marker = :rtriangle,markersize=themarkersize, color=:orange )
    scatter!(265,465, marker = :utriangle, markersize=themarkersize, color=:green )


    return fig

end



function SLf2(PL,Nx,Ny,dx,dy,x,y,Lx,Ly,R2L,SFreq;scale=1.0,makepngs=true)

    fig=Figure(resolution=(0.5*6.161, 6.2) .* (72*scale), fontsize=12)


    G1 = fig[1,1] = GridLayout()
    G2 = fig[1,2] = GridLayout()
    G3 = fig[2,1] = GridLayout()
    G4 = fig[2,2] = GridLayout()
    G5 = fig[3,1] = GridLayout()
    G6 = fig[3,2] = GridLayout()
    G7 = fig[4,1] = GridLayout()
    G8 = fig[4,2] = GridLayout()

    theticklabelsize = 8*scale
    xcut=0
    ycut=0
    theLW=1*scale

    P0  = zeros(3,Nx,Ny)
    R2 = zeros(3,3)

    
    P0 = PL[1,:,:,:]
    R2 = R2L[1,:,:]

    if makepngs == true
        export_arrow(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq,scale=1.0,filename="0_6.png")
    end
    img = load("0_6.png");
    size(img)

    ImX, ImY = size(img)


    ytks = [0,500/3,2*500/3,500]
    pytks = [0,1,2,3]

    image(G1[1, 1], img',

        axis = (aspect = DataAspect(),
        xticklabelsvisible=false,
        xticksvisible=false,
        ylabelpadding=-8,
        yticksize=3,
        yticksvisible=true,
        yreversed=true,
        ylabel=L"r",
        yticklabelsvisible=true,
        ylabelvisible=true,
        xlabelvisible=false,
        yticks = (ytks,[latexstring(pytks[a]) for a in 1:length(ytks)])
        
    ))

    ax2 = Axis(G1[1, 2]; aspect=Lx/Ly,
    xgridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    yticklabelsvisible=false,
    xlabelsize = 12*scale,
    xticklabelsize = 10*scale,
    xticksvisible=false,
    xticklabelsvisible=false
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue )


   colgap!(G1,1, 0.0*scale)
   Label(G1[1,1:2,Top()], L"N=0")




   P0 = PL[2,:,:,:]
   if makepngs == true
        export_arrow(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq,scale=1.0,filename="6_6.png")
   end
   img = load("6_6.png");



   image(G2[1, 1], img',

        axis = (aspect = DataAspect(),
        xticklabelsvisible=false,
        xticksvisible=false,
        ylabelpadding=-8,
        yticksize=3,
        yticksvisible=false,
        yreversed=true,
        yticklabelsvisible=false,
        ylabelvisible=false,
        xlabelvisible=false,
        yticks = (ytks,[latexstring(pytks[a]) for a in 1:length(ytks)])
        
    ))

    ax2 = Axis(G2[1, 2]; aspect=Lx/Ly,
    xgridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    yticklabelsvisible=false,
    xlabelsize = 12*scale,
    xticklabelsize = 10*scale,
    xticksvisible=false,
    xticklabelsvisible=false,
    ylabelvisible=true,
    ylabel=latexstring("\\mathbf{s}=(1,1,1),\\;\\; \\mathbf{r}=(1,-1,0)"),
    ylabelrotation=3pi/2,
    yaxisposition=:right
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue )

   colgap!(G2,1, 0.0*scale)
   Label(G2[1,1:2,Top()], L"N=1")


      P0 = PL[3,:,:,:]
    R2 = R2L[2,:,:,]

    if makepngs == true
   export_arrow(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq,scale=1.0,filename="1_6.png")
    end

   img = load("1_6.png");
    size(img)

    ImX, ImY = size(img)


    ytks = [0,500/3,2*500/3,500]
    pytks = [0,1,2,3]

    image(G3[1, 1], img',

        axis = (aspect = DataAspect(),
        xticklabelsvisible=false,
        xticksvisible=false,
        ylabelpadding=-8,
        yticksize=3,
        yticksvisible=true,
        yreversed=true,
        ylabel=L"r",
        yticklabelsvisible=true,
        ylabelvisible=true,
        xlabelvisible=false,
        yticks = (ytks,[latexstring(pytks[a]) for a in 1:length(ytks)])
        
    ))

    ax2 = Axis(G3[1, 2]; aspect=Lx/Ly,
    xgridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    yticklabelsvisible=false,
    xlabelsize = 12*scale,
    xticklabelsize = 10*scale,
    xticksvisible=false,
    xticklabelsvisible=false
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue )


   colgap!(G3,1, 0.0*scale)
   Label(G3[1,1:2,Top()], L"N=1/6")






   P0 = PL[4,:,:,:]

   if makepngs == true
   export_arrow(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq,scale=1.0,filename="5_6.png")
   end
   
   img = load("5_6.png");


   image(G4[1, 1], img',

        axis = (aspect = DataAspect(),
        xticklabelsvisible=false,
        xticksvisible=false,
        ylabelpadding=-8,
        yticksize=3,
        yticksvisible=false,
        yreversed=true,
        yticklabelsvisible=false,
        ylabelvisible=false,
        xlabelvisible=false,
        yticks = (ytks,[latexstring(pytks[a]) for a in 1:length(ytks)])
        
    ))

    ax2 = Axis(G4[1, 2]; aspect=Lx/Ly,
    xgridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    yticklabelsvisible=false,
    xlabelsize = 12*scale,
    xticklabelsize = 10*scale,
    xticksvisible=false,
    xticklabelsvisible=false
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue )

   colgap!(G4,1, 0.0*scale)
   Label(G4[1,1:2,Top()], L"N=5/6")








    R2 = R2L[3,:,:]
   P0 = PL[5,:,:,:]
   if makepngs == true
    export_arrow(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq,scale=1.0,filename="2_6.png")
   end 

    img = load("2_6.png");

    ImX, ImY = size(img)


    ytks = [0,500/3,2*500/3,500]
    pytks = [0,1,2,3]

    image(G5[1, 1], img',

        axis = (aspect = DataAspect(),
        xticklabelsvisible=false,
        xticksvisible=false,
        ylabelpadding=-8,
        yticksize=3,
        yticksvisible=true,
        yreversed=true,
        ylabel=L"r",
        yticklabelsvisible=true,
        ylabelvisible=true,
        xlabelvisible=false,
        yticks = (ytks,[latexstring(pytks[a]) for a in 1:length(ytks)])
        
    ))

    ax2 = Axis(G5[1, 2]; aspect=Lx/Ly,
    xgridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    yticklabelsvisible=false,
    xlabelsize = 12*scale,
    xticklabelsize = 10*scale,
    xticksvisible=false,
    xticklabelsvisible=false
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue )


   colgap!(G5,1, 0.0*scale)
   Label(G5[1,1:2,Top()], L"N=2/6")

   


   P0 = PL[6,:,:,:]

   if makepngs == true
        export_arrow(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq,scale=1.0,filename="4_6.png")
   end
    img = load("4_6.png");



   image(G6[1, 1], img',

        axis = (aspect = DataAspect(),
        xticklabelsvisible=false,
        xticksvisible=false,
        ylabelpadding=-8,
        yticksize=3,
        yticksvisible=false,
        yreversed=true,
        yticklabelsvisible=false,
        ylabelvisible=false,
        xlabelvisible=false,
        yticks = (ytks,[latexstring(pytks[a]) for a in 1:length(ytks)])
        
    ))

    ax2 = Axis(G6[1, 2]; aspect=Lx/Ly,
    xgridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    yticklabelsvisible=false,
    xlabelsize = 12*scale,
    xticklabelsize = 10*scale,
    xticksvisible=false,
    xticklabelsvisible=false
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue )

   xtks = [0,80,160,240]
   pxtks = [0,1,2,3]

   colgap!(G6,1, 0.0*scale)
   Label(G6[1,1:2,Top()], L"N=4/6")

    R2 = R2L[4,:,:]
   P0 = PL[7,:,:,:]
   if makepngs == true
    export_arrow(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq,scale=1.0,filename="3_6.png")
    end
img = load("3_6.png");

    ImX, ImY = size(img)




    image(G7[1, 1], img',

        axis = (aspect = DataAspect(),
        ylabelpadding=-8,
        yticksize=3,
        yticksvisible=true,
        yreversed=true,
        ylabel=L"r",
        yticklabelsvisible=true,
        ylabelvisible=true,
        xlabelpadding=-10,
        xlabel = L"s",
        yticks = (ytks,[latexstring(pytks[a]) for a in 1:length(ytks)]),
        xticks = (xtks,[latexstring(pxtks[a]) for a in 1:length(xtks)])
        
    ))

    ax2 = Axis(G7[1, 2]; aspect=Lx/Ly,
    xgridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    yticklabelsvisible=false,
    xlabelsize = 12*scale,
    xticklabelsize = 10*scale,
    xticksvisible=false,
    xticklabelsvisible=false
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue )


   colgap!(G7,1, 0.0*scale)
   Label(G7[1,1:2,Top()], L"N=1/2")




   P0 = PL[8,:,:,:]
   if makepngs == true
   export_arrow(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq,scale=1.0,filename="3_6_2.png")
   end
   img = load("3_6_2.png");



   image(G8[1, 1], img',

        axis = (aspect = DataAspect(),
        ylabelpadding=-8,
        yticksize=3,
        yticksvisible=false,
        yreversed=true,
        yticklabelsvisible=false,
        ylabelvisible=false,      
        xlabel = L"s",
        xlabelpadding=-10,  
        yticks = (ytks,[latexstring(pytks[a]) for a in 1:length(ytks)]),
        xticks = (xtks,[latexstring(pxtks[a]) for a in 1:length(xtks)])
        
    ))

    ax2 = Axis(G8[1, 2]; aspect=Lx/Ly,
    xgridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    yticklabelsvisible=false,
    xlabelsize = 12*scale,
    xticklabelsize = 10*scale,
    xticksvisible=false,
    xticklabelsvisible=false
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue )


   colgap!(G8,1, 0.0*scale)
   Label(G8[1,1:2,Top()], L"N=3/6")




   resize_to_layout!(fig)

    return fig

end








function SLEXT(;scale=1.0)

    fig=Figure(resolution=(0.5*6.161, 1.5) .* (72*scale), fontsize=12)



    img = load("3_6.png");
    size(img)

    ImX, ImY = size(img)

   # ytks = [0,500/3,2*500/3,500]
   # pytks = [0,1,2,3]

    image(fig[1,1], img',

        axis = (aspect = DataAspect(),
        xticklabelsvisible=false,
        xticksvisible=false,
        ylabelpadding=-8,
        yticksize=3,
        yticksvisible=true,
        yreversed=true,
        ylabel=L"r",
        yticklabelsvisible=true,
        ylabelvisible=true,
        xlabelvisible=false,
      #  yticks = (ytks,[latexstring(pytks[a]) for a in 1:length(ytks)])
        
    ))

    img = load("EXT.png");
    size(img)

    ImX, ImY = size(img)

   # ytks = [0,500/3,2*500/3,500]
   # pytks = [0,1,2,3]

    image(fig[1,2], img',

        axis = (aspect = DataAspect(),
        xticklabelsvisible=false,
        xticksvisible=false,
        ylabelpadding=-8,
        yticksize=3,
        yticksvisible=false,
        yreversed=true,
        ylabel=L"r",
        yticklabelsvisible=false,
        ylabelvisible=false,
        xlabelvisible=false,
      #  yticks = (ytks,[latexstring(pytks[a]) for a in 1:length(ytks)])
        
    ))

    poly!(Point2f[(50,190),(50,320),(190,320),(190,190)], color=(colorant"rgb(255,0,255)",0.5)
    , strokecolor=colorant"rgb(255,0,255)", strokewidth=3.0
    
    )


    img = load("3_6_2.png");
    size(img)

    ImX, ImY = size(img)

   # ytks = [0,500/3,2*500/3,500]
   # pytks = [0,1,2,3]

    image(fig[1,3], img',

        axis = (aspect = DataAspect(),
        xticklabelsvisible=false,
        xticksvisible=false,
        ylabelpadding=-8,
        yticksize=3,
        yticksvisible=false,
        yreversed=true,
        ylabel=L"r",
        yticklabelsvisible=false,
        ylabelvisible=false,
        xlabelvisible=false,
      #  yticks = (ytks,[latexstring(pytks[a]) for a in 1:length(ytks)])
        
    ))






    colgap!(fig.layout,1, -20*scale)
    colgap!(fig.layout,2, -20*scale)

    return fig




end

















function SLden(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq,LS;scale=1.0)

    fig=Figure(resolution=(0.5*6.161, 2) .* (72*scale), fontsize=12)
    BD = zeros(Nx,Ny)
    baryon(P0, x,y,BD=BD)

    G1 = fig[1,1] = GridLayout()


    theticklabelsize = 8*scale
    xcut=0
    ycut=0
    theLW=1*scale


    export_arrow(P0[:,31:131,:],101,Ny,dx,dy,x,y,Lx*101/161,Ly,R2,SFreq,scale=1.0,filename="f2_6.png")

    img = load("f2_6.png");
    size(img)

    ImX, ImY = size(img)
    println(ImX, ", ", ImY)

    xtks = [0,340/3,2*340/3,340]
    pxtks = round.(xtks.*(LS*dx*101/ImY), digits=2)


    ytks = [0,500/3,2*500/3,500]
    pytks = reverse(round.(ytks.*(LS*dy*161/ImX), digits=2))

    image(G1[1, 1], img',

        axis = (aspect = DataAspect(),

        ylabelpadding=-8,
        yticksize=3,
        yticksvisible=true,
        yreversed=true,
        ylabel=L"r",
        ylabelsize=8*scale,
        xlabelsize=8*scale,
        yticklabelsvisible=true,
        ylabelvisible=true,
        xautolimitmargin = (0.0,0.0),
        yticks = (ytks,[latexstring(pytks[a]) for a in 1:length(ytks)]),
        xlabel = L"s",
        xticks = (xtks,[latexstring(pxtks[a]) for a in 1:length(xtks)])
        
    ))
    ytks = ([0,500/3,2*500/3,500] .- 250)./500*Ly 

    ax2 = Axis(G1[1, 2]; aspect=101/161*Lx/Ly,
    xgridvisible = false,
    yautolimitmargin = (0.01,0.01),
    xlabelsize = 12*scale,
    xticklabelsize = 10*scale,
    yticksvisible=true,
    xticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
    yticks = ytks
        );

    xcut = 28

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    

    hm=heatmap!(ax2,x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut],BD[-2+3+xcut:-2+Nx-2-xcut,3:+Ny-2]
    ,colormap=Reverse(:CMRmap))


   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue )

   Colorbar(G1[1,3],hm, tellheight=true, label=L"\text{Skyrme charge density}", labelrotation=3*pi/2)
   #Label(G1[1,1:2,Top()], L"N=0")


   colgap!(G1,1, -30.0*scale)
   colgap!(G1,2, -15.0*scale)

   resize_to_layout!(fig)

    return fig

end





function stringfig(EL,PL,R2,x,y,stL;scale=3.0)

    Nx = size(x)[1]
    Ny = size(y)[1]

    ytks = [11.8,12.0,12.2,12.4]
    pytks = [ latexstring(ytks[a]) for a in 1:4 ]

    xtks = [2,4,6,8]
    pxtks = [ latexstring(xtks[a]) for a in 1:4 ]


    fig=Figure(resolution=(0.5*6.161, 2) .* (72*scale), fontsize=12,
    
    #Axis = ( 
    #yticks = (ytks,pytks),
    #xticks = (xtks,pxtks)
    #)
    
    
    )

    ax1 = Axis(fig[1,1:5];
     ylabel=L"\text{Energy}(n)",
     yticks = (ytks,pytks),
    xticks = (xtks,pxtks)
     )

    lines!(ax1, EL )



    theLW = 1

    c1 = colorant"rgb(125,255,125)"
    c2 = colorant"rgb(255,0,255)"
    c3 = colorant"rgb(0,255,255)"
    c4 = colorant"rgb(255,125,125)"
    c5 = colorant"rgb(125,125,255)"

    cL = [c1,c2,c3,c4,c5]

    for a in 1:5
        scatter!(ax1, stL[a],EL[stL[a]],color=cL[a])
    end

    arrows!([1.2],[11.9],[-0.14],[EL[1]-11.9+0.06])
    text!(1.2,11.905, text=L"N = 0")

    arrows!([8.8],[12.1],[9.0-8.8-0.06],[EL[9]-12.1-0.06])
    text!(8.8-0.3,12.1-0.1, text=L"N = 1")

    arrows!([5.5],[12.2],[6-5.5-0.04],[EL[6]-12.2-0.04])
    text!(5.5-0.6,12.2-0.1, text=L"N \text{ undefined}")



    xcut=34
    ycut=0

    ax2 = Axis(fig[2,1]; aspect=((Nx-2*xcut)/Ny),

    rightspinecolor = c1,
    leftspinecolor = c1,
    topspinecolor = c1,
    bottomspinecolor = c1,

    xgridvisible = false,
    ygridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    xticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'PL[1,:,i,j]
    end    
    contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
    , linewidth = theLW, color = :red)
    contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :green)
    contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :blue )


    ax3 = Axis(fig[2,2]; aspect=((Nx-2*xcut)/Ny),
    rightspinecolor = c2,
    leftspinecolor = c2,
    topspinecolor = c2,
    bottomspinecolor = c2,
    xgridvisible = false,
    ygridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    xticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'PL[2,:,i,j]
    end    
    contour!( ax3, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
    , linewidth = theLW, color = :red)
    contour!( ax3, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :green)
    contour!( ax3, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :blue )


    ax4 = Axis(fig[2,3]; aspect=((Nx-2*xcut)/Ny),
    rightspinecolor = c3,
    leftspinecolor = c3,
    topspinecolor = c3,
    bottomspinecolor = c3,
    xgridvisible = false,
    ygridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    xticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'PL[3,:,i,j]
    end    
    contour!( ax4, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
    , linewidth = theLW, color = :red)
    contour!( ax4, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :green)
    contour!( ax4, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :blue )



    ax5 = Axis(fig[2,4]; aspect=((Nx-2*xcut)/Ny),
    rightspinecolor = c4,
    leftspinecolor = c4,
    topspinecolor = c4,
    bottomspinecolor = c4,
    xgridvisible = false,
    ygridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    xticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'PL[4,:,i,j]
    end    
    contour!( ax5, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
    , linewidth = theLW, color = :red)
    contour!( ax5, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :green)
    contour!( ax5, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :blue )
        



    ax6 = Axis(fig[2,5]; aspect=((Nx-2*xcut)/Ny),
    xgridvisible = false,
    rightspinecolor = c5,
    leftspinecolor = c5,
    topspinecolor = c5,
    bottomspinecolor = c5,
    ygridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    xticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'PL[5,:,i,j]
    end    
    contour!( ax6, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
    , linewidth = theLW, color = :red)
    contour!( ax6, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :green)
    contour!( ax6, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :blue )
        

    rowgap!(fig.layout,1, -0.0*scale)


    return fig




end







function stringfig2(EL,PL,R2,x,y,stL;scale=3.0)

    Nx = size(x)[1]
    Ny = size(y)[1]

    ytks = [9.95,10.05,10.15]
    pytks = [ latexstring(ytks[a]) for a in 1:3 ]

    xtks = [2,4,6]
    pxtks = [ latexstring(xtks[a]) for a in 1:3 ]


    fig=Figure(resolution=(0.5*6.161, 2) .* (72*scale), fontsize=12 )

    ax1 = Axis(fig[1,1:5];
     ylabel=L"\text{Energy}(n)",
     yticks = (ytks,pytks),
    xticks = (xtks,pxtks)
     )

    lines!(ax1, EL )



    theLW = 1

    c1 = colorant"rgb(125,255,125)"
    c2 = colorant"rgb(255,0,255)"
    c3 = colorant"rgb(0,255,255)"
    c4 = colorant"rgb(255,125,125)"
    c5 = colorant"rgb(125,125,255)"

    cL = [c1,c2,c3,c4,c5]

    for a in 1:5
        scatter!(ax1, stL[a],EL[stL[a]],color=cL[a])
    end


    xcut=34
    ycut=0

    ax2 = Axis(fig[2,1]; aspect=((Nx-2*xcut)/Ny),

    rightspinecolor = c1,
    leftspinecolor = c1,
    topspinecolor = c1,
    bottomspinecolor = c1,

    xgridvisible = false,
    ygridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    xticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'PL[1,:,i,j]
    end    
    contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
    , linewidth = theLW, color = :red)
    contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :green)
    contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :blue )


    ax3 = Axis(fig[2,2]; aspect=((Nx-2*xcut)/Ny),
    rightspinecolor = c2,
    leftspinecolor = c2,
    topspinecolor = c2,
    bottomspinecolor = c2,
    xgridvisible = false,
    ygridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    xticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'PL[2,:,i,j]
    end    
    contour!( ax3, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
    , linewidth = theLW, color = :red)
    contour!( ax3, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :green)
    contour!( ax3, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :blue )


    ax4 = Axis(fig[2,3]; aspect=((Nx-2*xcut)/Ny),
    rightspinecolor = c3,
    leftspinecolor = c3,
    topspinecolor = c3,
    bottomspinecolor = c3,
    xgridvisible = false,
    ygridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    xticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'PL[3,:,i,j]
    end    
    contour!( ax4, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
    , linewidth = theLW, color = :red)
    contour!( ax4, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :green)
    contour!( ax4, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :blue )



    ax5 = Axis(fig[2,4]; aspect=((Nx-2*xcut)/Ny),
    rightspinecolor = c4,
    leftspinecolor = c4,
    topspinecolor = c4,
    bottomspinecolor = c4,
    xgridvisible = false,
    ygridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    xticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'PL[4,:,i,j]
    end    
    contour!( ax5, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
    , linewidth = theLW, color = :red)
    contour!( ax5, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :green)
    contour!( ax5, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :blue )
        



    ax6 = Axis(fig[2,5]; aspect=((Nx-2*xcut)/Ny),
    xgridvisible = false,
    rightspinecolor = c5,
    leftspinecolor = c5,
    topspinecolor = c5,
    bottomspinecolor = c5,
    ygridvisible = false,
    yautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    xticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'PL[5,:,i,j]
    end    
    contour!( ax6, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
    , linewidth = theLW, color = :red)
    contour!( ax6, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :green)
    contour!( ax6, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
    , linewidth = theLW, color = :blue )
        

    rowgap!(fig.layout,1, -0.0*scale)


    return fig




end



















function export_arrow(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq;scale=1.0, filename="temp.png")



    xlr = x[3:SFreq:Nx-2];
    ylr = y[3:SFreq:Ny-2];

    Plr0 = P0[:,3:SFreq:Nx-2,3:SFreq:Ny-2]
    Plr = 0.0.*Plr0

    for i in 1:size(Plr0)[2], j in 1:size(Plr0)[3]
        Plr[:,i,j] .= R2'*Plr0[:,i,j]
    end


    ps = [Point3f(xlr[i], ylr[j], 0.0) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];
    ns = [Point3f(Plr[1,i,j],Plr[2,i,j], Plr[3,i,j]) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];


    newlengths = sortlengths(ns)

    fig = Figure(fontsize=12, resolution=(1.15*Lx/Ly*600,600).*scale, )


# FALSE PARADISE
                    
        axs = Axis3(fig[1, 1]; aspect=(Lx,Ly,1.0*min(Lx,Ly)), 
            perspectiveness=0.0 , azimuth=-pi/2, elevation = pi/2,

        zgridvisible = false,
        zspinesvisible=false,
        zticksvisible=false,
        zticklabelsvisible=false,
        zlabelvisible=false,


        tellwidth=true,tellheight=true,
        
        xspinesvisible = false,
        xticksvisible = false,
        xlabelvisible = false,
        xgridvisible = false,
        xticklabelsize = 10*scale,
        xticklabelsvisible = false,

        ylabelvisible = false,
        yspinesvisible = false,
        yticksvisible = false,
        ygridvisible = false,
        yautolimitmargin = (0.02,0.02),
        yticklabelsvisible = false

        ) 




    ascale = 1.0*scale

    arrows!(axs, ps, ns./10.0,  arrowsize=Vec3f(0.2*ascale, 0.2*ascale, 0.2), linewidth=0.1,  #arrowtail = Vec3f(1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0)),
        color= newlengths,
        align=:center)



    #ablines!([1, 2, 3], [1, 1.5, 2], color = [:red, :orange, :pink], linestyle=:dash, linewidth=2)


#cam = cameracontrols(axs.scene)
#update_cam!(axs.scene, [5,5,-10])

resize_to_layout!(fig)

   # zlims!(-1.5,1.5)


    save(filename,axs.scene)
    #save(filename,fig)


    return axs.scene


end


















function arrow_and_contour_and_surface_plot(P0,Nx,Ny,dx,dy,x,y,Lx,Ly,R2,SFreq;scale=1.0,savefile=false,filename="temp.png")

    BD = zeros(Nx,Ny)
    baryon(P0, x,y,BD=BD)

    xlr = x[3:SFreq:Nx-2];
    ylr = y[3:SFreq:Ny-2];

    Plr0 = P0[:,3:SFreq:Nx-2,3:SFreq:Ny-2]
    Plr = 0.0.*Plr0

    for i in 1:size(Plr0)[2], j in 1:size(Plr0)[3]
        Plr[:,i,j] .= R2'*Plr0[:,i,j]
    end


    ps = [Point3f(xlr[i], ylr[j], 0.0) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];
    ns = [Point3f(Plr[1,i,j],Plr[2,i,j], Plr[3,i,j]) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];


    newlengths = sortlengths(ns)

    fig = Figure(resolution=(2.0*Lx/Ly*600,600).*scale, fontsize=12)


    axs = Axis3(fig[1, 1]; aspect=(Lx,Ly,1.0*min(Lx,Ly)), 
        perspectiveness=0.0 , azimuth=-pi/2-0.00001, elevation = pi/2,
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        zspinesvisible=false,
        zticksvisible=false,
        zticklabelsvisible=false,
        zlabelvisible=false,
        yautolimitmargin = (0.02,0.02),
        xlabelsize = 12*scale,
        ylabelsize = 12*scale,
        ylabel=L"r",
        xlabel=L"s",
        tellwidth=true,tellheight=true,
        xticklabelsize = 10*scale,
        yticklabelsize = 10*scale
        
        ) 

    ascale = 1.0*scale

    arrows!(axs, ps, ns./10.0,  arrowsize=Vec3f(0.2*ascale, 0.2*ascale, 0.2), linewidth=0.1,  #arrowtail = Vec3f(1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0)),
        color= newlengths,
        align=:center)

    zlims!(-1.5,1.5)



    xcut=0

    ax2 = Axis3(fig[1, 2]; aspect=(Lx,Ly,1.0*min(Lx,Ly)), 
    perspectiveness=0.0 , azimuth=-pi/2-0.00001, elevation = pi/2,
    xgridvisible = false,
    #ygridvisible = false,
    zgridvisible = false,
    zspinesvisible=false,
    zticksvisible=false,
    zticklabelsvisible=false,
    zlabelvisible=false,
    yautolimitmargin = (0.01,0.01),
    xautolimitmargin = (0.01,0.01),
    yticksvisible=false,
    yticklabelsvisible=false,
    ylabelvisible=false,
    xlabelsize = 12*scale,
    xlabel=L"s",
    xticklabelsize = 10*scale,



        );

    P0n = zeros(3,Nx,Ny)
    for i in 1:Nx, j in 1:Ny
        P0n[:,i,j] = R2'P0[:,i,j]
    end    
    
    zlims!(-1.5,1.5)
    ycut=0

    hm=heatmap!(ax2,x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut],BD[3:Nx-2,3:Ny-2],colormap=:Blues_9)

 

    theLW = 5

   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[1,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0] 
   , linewidth = theLW, color = :red)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[2,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :green)
   contour!( ax2, x[3+xcut:Nx-2-xcut], y[3+ycut:Ny-2-ycut], P0n[3,3+xcut:Nx-2-xcut,3+ycut:Ny-2-ycut], levels=[0.0]
   , linewidth = theLW, color = :blue )


   Colorbar(fig[1,3],hm, tellheight=true)


   colgap!(fig.layout,1, -0*scale)
   colgap!(fig.layout,2, -0*scale)

    if savefile == true
        imusingnotebook()
        save(filename, fig)
        #imusingterminal()
    end



    return fig

end






