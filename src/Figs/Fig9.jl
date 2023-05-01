
function flowRAKfig9!(phi,p)
    
    tf = p[10]
    EList = zeros(tf)

    # p = (N,dx,G2,A2,A4,R2, a111, a112, a123,tf,Em0p,V0)

    tspan = (0.0, tf)
    engtimes = [ 1.0*i for i=1:tf ]

    cb = PresetTimeCallback(engtimes, affect_checkeng!, save_positions=(false,true))
    #cb = DiscreteCallback(condition,affect!, save_positions=(false,true))

    prob = ODEProblem(dEduf!,phi,tspan,p)
    sol = solve(prob, Tsit5(),reltol=1e-8, abstol=1e-8, save_everystep=false, callback=cb)
    #sol = solve(prob, Tsit5(),reltol=1e-8, abstol=1e-8, save_everystep=true, callback=cb, tstops = engtimes )

    #if sol.t[end] == tf
    #   println("AGGGG");
    #else
    #    println(sol.t[end])
    #end

    return sol
    
end




function makefig9(Pplot1, Pplot2, Pplot3, Pplot4,Pplot5, Pplot6, Pplot7, Pplot8, x, N,R2, G2, A2, A4, a111, a112, a123, V0; scale=1.0)

    spcolor = RGBf(0.7,0.7,0.7)



    f = Figure(backgroundcolor = :white, resolution = (0.5*6.161, 3.0) .* (72*scale), fontsize=8*scale)


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
    
    
    xtks = [-3, 0, 3]
    ytks = [-0.7,0,0.7]
        
    if length(xtks) != 0
        ax1.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax1.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end

    ylims!(ax1, -0.9, 0.9)
    xlims!(ax1, -4, 4)


    ax2 = Axis(f[2,1],
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

    ys = [ Pplot2[a,:] for a in 1:3 ]
    ps = [ lines!(ax2, xp, ys[a], label=thelabels[a], linewidth=1.0*scale, color=thecolors[a]) for a in 1:3 ]
    
    
    xtks = [-3, 0, 3]

        
    if length(xtks) != 0
        ax2.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax2.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end



    ylims!(ax2, -0.9, 0.9)
    xlims!(ax2, -4, 4)



    ax3 = Axis(f[3,1],
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

    ys = [ Pplot3[a,:] for a in 1:3 ]
    ps = [ lines!(ax3, xp, ys[a], label=thelabels[a], linewidth=1.0*scale, color=thecolors[a]) for a in 1:3 ]
    
    
    xtks = [-3, 0, 3]

        
    if length(xtks) != 0
        ax3.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax3.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end



    ylims!(ax3, -0.9, 0.9)
    xlims!(ax3, -4, 4)



    ax4 = Axis(f[4,1],
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

    ys = [ Pplot4[a,:] for a in 1:3 ]
    ps = [ lines!(ax4, xp, ys[a], label=thelabels[a], linewidth=1.0*scale, color=thecolors[a]) for a in 1:3 ]
    
    
    xtks = [-3, 0, 3]

        
    if length(xtks) != 0
        ax4.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax4.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end



    ylims!(ax4, -0.9, 0.9)
    xlims!(ax4, -4, 4)








    ax5 = Axis(f[1,2],
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
        xgridvisible = false,
        yticklabelsvisible=false

    )
        
   thelabels = [L"P_x", L"P_y", L"P_z" ]
    
    xp = x

    thecolors = [ colorant"#1f77b4", colorant"#ff7f0e", colorant"#2ca02c" ]
    thelinestyles = [:solid,:solid,:solid]

    ys = [ Pplot5[a,:] for a in 1:3 ]
    ps = [ lines!(ax5, xp, ys[a], label=thelabels[a], linewidth=1.0*scale, color=thecolors[a]) for a in 1:3 ]
    
    
    xtks = [-3, 0, 3]

        
    if length(xtks) != 0
        ax5.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax5.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end

    ylims!(ax5, -0.9, 0.9)
    xlims!(ax5, -4, 4)


    ax6 = Axis(f[2,2],
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
        xgridvisible = false,
        yticklabelsvisible=false

    )
        
   thelabels = [L"P_x", L"P_y", L"P_z" ]
    
    xp = x

    thecolors = [ colorant"#1f77b4", colorant"#ff7f0e", colorant"#2ca02c" ]
    thelinestyles = [:solid,:solid,:solid]

    ys = [ Pplot6[a,:] for a in 1:3 ]
    ps = [ lines!(ax6, xp, ys[a], label=thelabels[a], linewidth=1.0*scale, color=thecolors[a]) for a in 1:3 ]
    
    
    xtks = [-3, 0, 3]

        
    if length(xtks) != 0
        ax6.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax6.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end



    ylims!(ax6, -0.9, 0.9)
    xlims!(ax6, -4, 4)



    ax7 = Axis(f[3,2],
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
        xgridvisible = false,
        yticklabelsvisible=false

    )
        
   thelabels = [L"P_x", L"P_y", L"P_z" ]
    
    xp = x

    thecolors = [ colorant"#1f77b4", colorant"#ff7f0e", colorant"#2ca02c" ]
    thelinestyles = [:solid,:solid,:solid]

    ys = [ Pplot7[a,:] for a in 1:3 ]
    ps = [ lines!(ax7, xp, ys[a], label=thelabels[a], linewidth=1.0*scale, color=thecolors[a]) for a in 1:3 ]
    
    
    xtks = [-3, 0, 3]

        
    if length(xtks) != 0
        ax7.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax7.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end



    ylims!(ax7, -0.9, 0.9)
    xlims!(ax7, -4, 4)



    ax8 = Axis(f[4,2],
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
        xgridvisible = false,
        yticklabelsvisible=false

    )
        
   thelabels = [L"P_x", L"P_y", L"P_z" ]
    
    xp = x

    thecolors = [ colorant"#1f77b4", colorant"#ff7f0e", colorant"#2ca02c" ]
    thelinestyles = [:solid,:solid,:solid]

    ys = [ Pplot8[a,:] for a in 1:3 ]
    ps = [ lines!(ax8, xp, ys[a], label=thelabels[a], linewidth=1.0*scale, color=thecolors[a]) for a in 1:3 ]
    
    
    xtks = [-3, 0, 3]

        
    if length(xtks) != 0
        ax8.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax8.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end



    ylims!(ax8, -0.9, 0.9)
    xlims!(ax8, -4, 4)






    #colgap!(f.layout,1, 10)

    return f
end