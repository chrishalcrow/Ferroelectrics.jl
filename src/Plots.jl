


function plot_field(P; 
	savefile = false,
	filename = "temp.pdf",
	scale::Float64=1.0,
	xlimL = P.x[1],
	xlimR = P.x[end],
	ylimT = maximum(P.field)*1.2,
	ylimB = minimum(P.field)*1.2,
	titlename = "",
	xtks = [],
	ytks = []
	)
	

	
	println(xtks)

    spcolor = RGBf(0.7,0.7,0.7)

    f1 = Figure(backgroundcolor = :white, resolution = (6.161, 2.4) .* (72*scale), fontsize=10*scale)

    ax = Axis(f1[1,1],
		title = titlename,
		#title = titlename,
        bottomspinecolor = spcolor,
        topspinecolor = spcolor,
        leftspinecolor = spcolor,
        rightspinecolor = spcolor,
        xtickcolor = spcolor,
        ytickcolor = spcolor,
        titlesize = 12*scale
    )
	
	if length(xtks) != 0
		ax.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
	end
	if length(ytks) != 0
		ax.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
	end

	
    xp = P.x
	ys = [ P.field[a,:] for a in 1:3 ]

	ps = [ lines!(ax, xp, ys[a], linewidth=1.5*scale) for a in 1:3 ] 

    ylims!(ax, ylimB, ylimT)
    xlims!(ax, xlimL, xlimR)

	Legend(f1[1,2],
	    framevisible = false,
	    ps,
	    [L"P_1(s)",L"P_2(s)",L"P_3(s)"],
        padding = (0.0f0, 0.0f0, 0.0f0, 0.0f0),
		linepoints = [Point2f(1.0-scale, 0.5), Point2f(1.0, 0.5)]
    )
	
   
   	colgap!(f1.layout,20*scale) 

   	if savefile == true
   		save(filename,f1)
   		return
   	end

    return f1
	
end

function plot_field(phi, x; 
	savefile = false,
	filename = "temp.pdf",
	scale::Float64=3.0,
	xlimL = x[1],
	xlimR = x[end],
	ylimT = maximum(phi)*1.2,
	ylimB = minimum(phi)*1.2,
	titlename = "",
	xtks = [],
	ytks = []
	)
	
	
	println(xtks)

    spcolor = RGBf(0.7,0.7,0.7)

    f1 = Figure(backgroundcolor = :white, resolution = (6.161, 2.4) .* (72*scale), fontsize=10*scale)

    ax = Axis(f1[1,1],
		title = titlename,
		#title = titlename,
        bottomspinecolor = spcolor,
        topspinecolor = spcolor,
        leftspinecolor = spcolor,
        rightspinecolor = spcolor,
        xtickcolor = spcolor,
        ytickcolor = spcolor,
        titlesize = 12*scale
    )
	
	if length(xtks) != 0
		ax.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
	end
	if length(ytks) != 0
		ax.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
	end

	
    xp = x
	ys = [ phi[a,:] for a in 1:3 ]

	ps = [ lines!(ax, xp, ys[a], linewidth=1.5*scale) for a in 1:3 ] 

    ylims!(ax, ylimB, ylimT)
    xlims!(ax, xlimL, xlimR)

	Legend(f1[1,2],
	    framevisible = false,
	    ps,
	    [L"P_1(s)",L"P_2(s)",L"P_3(s)"],
        padding = (0.0f0, 0.0f0, 0.0f0, 0.0f0),
		linepoints = [Point2f(1.0-scale, 0.5), Point2f(1.0, 0.5)]
    )
	
   
   	colgap!(f1.layout,20*scale) 

   	if savefile == true
   		save(filename,f1)
   		return
   	end

    return f1
	
end