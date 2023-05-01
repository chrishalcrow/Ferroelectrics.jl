

function makefig11(Pplot1, x, N,R2, G2, A2, A4, a111, a112, a123, V0; scale=1.0)

    spcolor = RGBf(0.7,0.7,0.7)



    f = Figure(backgroundcolor = :white, resolution = (0.5*6.161, 0.9) .* (72*scale), fontsize=8*scale)

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

    xtks = [-8,-4,0,4,8]
    ytks = [-1.0,-0.8,-0.5,0.0,0.5,0.8]
        
    if length(xtks) != 0
        ax1.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax1.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end



    ylims!(ax1, -0.75, 0.75)
    xlims!(ax1, -10, 10)






    #colgap!(f.layout,1, 10)

    return f
end