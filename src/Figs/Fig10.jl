function getfig10data(N,dx)

    x = setgrid(N,dx)

   # FIGURE 10 stuff

    ss = [ 1.0,1.0,1.0]./sqrt(3.0)
    rr = [ 0.0, 1.0, -1.0 ]./sqrt(2.0)

    R2 = getR2(ss,rr)

    A2 = zeros(3,3)
    G2 = zeros(3,3)

    A4 = zeros(3,3,3,3)

    (PV, A0, a111, a112, a123, V0) = set_parameters!(G2,A2,A4, 250, R2, "oct");


    A2s = zeros(3,3)
    A4s = zeros(3,3,3,3)

    for a in 1:3, b in 1:3
        A2s[a,b] = (A2[a,b] + A2[b,a])/2.0
        for c in 1:3, d in 1:3
            A4s[a,b,c,d] = (A4[a,b,c,d] + A4[b,a,c,d] + A4[b,c,a,d] + A4[b,c,d,a])/4.0
        end
    end


    tf=100
    EList=zeros(tf)

    scaler=0.6

    x1L = -2.0*scaler
    x1R = 1.0*scaler

    x2L = -1.0*scaler
    x2R = 2.0*scaler

    Pinxyz = zeros(3,N)

    Pinxyz[1,:] = (1.0 .+ tanh.(x .+ x1L) - tanh.(x .+ x1R))./sqrt(2.0)
    Pinxyz[2,:] = (1.0 .+ tanh.(x .+ x2L) - tanh.(x .+ x2R))./sqrt(2.0)
    Pinxyz[3,:] = 0.0.*x 


    dedp = zeros(3,N)

    Pinsrt = srtfxyz(Pinxyz,R2)

    sol = flowRAKf!(Pinsrt, (N,dx,G2,A2s,A4s,R2, a111, a112, a123,tf,EList,V0,dedp) )

    return sol

end

function makefig10(Pplot1, x, N,R2, G2, A2, A4, a111, a112, a123, V0; scale=1.0)

    spcolor = RGBf(0.7,0.7,0.7)



    f = Figure(backgroundcolor = :white, resolution = (1.0*0.5*6.161, 1.0) .* (72*scale), fontsize=8*scale)

    ax1 = Axis(f[1,1],
        #title = L"6+6=12",
        bottomspinecolor = spcolor,
        topspinecolor = spcolor,
        leftspinecolor = spcolor,
        rightspinecolor = spcolor,
        xtickcolor = spcolor,
        ytickcolor = spcolor,
        titlesize = 12*scale,
      #  xlabel = L"x",
        xlabelpadding = 0,
        xgridwidth = 0.5*scale,
        ygridwidth = 0.5*scale,
        ygridvisible = false,
        xgridvisible = false

    )
        
   thelabels = [L"P_x", L"P_y", L"P_z" ]
    
    xp = x

    thecolors = [ colorant"#1f77b4", colorant"#ff7f0e", colorant"#2ca02c" ]
    thelinestyles = [:solid,:solid,:solid]

    ys = [ Pplot1[a,:] for a in 1:3 ]
    ps = [ lines!(ax1, xp, ys[a], label=thelabels[a], linewidth=1.0*scale, color=thecolors[a]) for a in 1:3 ]
    
    axislegend(position = :lb, padding=[10,10,10,10],rowgap=0)

    xtks = [-15,-10,-5, 0, 5,10,15]
    ytks = [-1.0,-0.8,-0.4,0.0,0.4,0.8]
        
    if length(xtks) != 0
        ax1.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax1.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end



    ylims!(ax1, -0.9, 0.9)
    xlims!(ax1, -17, 17)






    #colgap!(f.layout,1, 10)

    return f
end