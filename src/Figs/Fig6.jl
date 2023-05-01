function getfig6data(PV, N, dx, G2, A2, A4, R2, a111, a112, a123, V0, dt, n)

    x = setgrid(N,dx)

    DP = zeros(3,N)
    P = ConnectVacThroughAnotherVac(x, PV, -[1.5,1.0,1.0])

    flow!(P, DP, N, dx, G2, A2, A4, R2, a111, a112, a123, V0, dt, n) 

    return P

end

function makefig6(Pplot1, x, N,R2, G2, A2, A4, a111, a112, a123, V0,PN,PB; scale=1.0)

    spcolor = RGBf(0.7,0.7,0.7)

    
    xc = -1.7:0.02:1.7
    yc =  -1.7:0.02:1.7
    VN = zeros(size(xc)[1],size(yc)[1])

    for a in 1:size(xc)[1], b in 1:size(yc)[1]
        VN[a,b] = min(3,PotEng( xc[a]*PN + yc[b]*PB, R2, G2, A2, A4, a111, a112, a123, V0 ))
    end




    f = Figure(backgroundcolor = :white, resolution = (0.45*6.161, 1.3) .* (72*scale), fontsize=8*scale)

    ax1 = Axis(f[1,1],
        #title = L"6+6=12",
        bottomspinecolor = spcolor,
        topspinecolor = spcolor,
        leftspinecolor = spcolor,
        rightspinecolor = spcolor,
        xtickcolor = spcolor,
        ytickcolor = spcolor,
        titlesize = 12*scale,
        xlabel = L"x",
        xlabelpadding = 0,
        ytickalign = 1.0,
        xgridwidth = 0.5*scale,
        ygridwidth = 0.5*scale

    )
        
   thelabels = [L"P_I", L"P_N", L"P_B" ]
    
    xp = x

    thelinestyles = [:solid,:solid,:solid]

    ys = [ Pplot1[a,:] for a in 1:3 ]
    ps = [ lines!(ax1, xp, ys[a], label=thelabels[a], linewidth=1.0*scale, linestyle=thelinestyles[a]) for a in 1:3 ]
    
    axislegend(position = :rb, padding=[10,10,3,3],rowgap=0)

    xtks = [-8,-4,0,4,8]
    ytks = [-1.0,-0.5,0.0,0.5,1.0]
        
    if length(xtks) != 0
        ax1.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax1.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end

    scatter!(xp[Int((N-1)/2)],0.0, marker = :star5, color=:red, markersize=8*scale)

    ylims!(ax1, -1.1, 1.1)
    xlims!(ax1, -3, 3)

    ax2 = Axis(f[1,2],
        ylabelpadding = -scale*5.0,
        xlabelpadding = 0.0,
ytickalign = 1.0,
# yticklabelsvisible = false,
        bottomspinecolor = spcolor,
        topspinecolor = spcolor,
        leftspinecolor = spcolor,
        rightspinecolor = spcolor,
        xtickcolor = spcolor,
        ytickcolor = spcolor,
xgridwidth = 0.5*scale,
        ygridwidth = 0.5*scale,
    xlabel = L"P_N",
    ylabel = L"P_B"
    )

   # minL = 0.1
   # maxL = 7.0
   # nlev = 16

    ylims!(ax2, -1.55, 1.55)
   xlims!(ax2, -1.1,1.1)


    xtks = [-1,0,1]
    ytks = [-1,0,1]

    if length(xtks) != 0
        ax2.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax2.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end
       

    hm = contour!(ax2,xc,yc,VN,
        colormap = :rainbow1,
    levels = [0.2,0.6,1.0,1.5,2.0,2.9],
      # levels= minL:(maxL-minL)/nlev:maxL,
        linewidth = 0.5*scale
    )

       # scatter!(0.0,0.0, marker = :star5, color=:blue, markersize=4*scale)
    scatter!(0.0,1.0, marker = :star5, color=:red, markersize=4*scale)

    #colbtks = [1.12 + (i-1)*0.2 for i in 1:10]
    #=
    colb = Colorbar(f[1,3],hm
     #   ,size = 5*scale
    #  ,vertical=false
    #€  ,ticks = (colbtks,[latexstring(round(colbtks[a],digits=1)) for a in 1:length(colbtks)])
        ,tellheight = true
        ,tellwidth = true
     #   ,flipaxis = false
    )
=#



    colgap!(f.layout,1, 15)

    return f
end