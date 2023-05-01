function getfig8data(N,dx)

    x = setgrid(N,dx)

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


    tf=1000
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


    Pinsrt = srtfxyz(Pinxyz,R2)

    dedp = zeros(3,N)

    sol = xyzfsrt(flowRAKf!(Pinsrt, (N,dx,G2,A2s,A4s,R2, a111, a112, a123,tf,EList,V0,dedp) ),R2)


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

D2V = zeros(3,3)
Ppt = srtfxyz([-1.0/sqrt(2.0),-1.0/sqrt(2.0),0.0],R2)

D2V6o = zeros(3,3)
D2V6 = zeros(3,3)

for a in 1:3, b in 1:3    
    
    D2V[a,b] += 2.0*A2ss[a,b]

    for c in 1:3, d in 1:3

        D2V[a,b] += 12.0*A4ss[a,b,c,d]*Ppt[c]*Ppt[d]

    end
    
    D2V6o[1,1] = d11V6(R2'*Ppt, a111,a112,a123)
    D2V6o[1,2] = d12V6(R2'*Ppt, a111,a112,a123)
    D2V6o[1,3] = d13V6(R2'*Ppt, a111,a112,a123)
    
    D2V6o[2,1] = d21V6(R2'*Ppt, a111,a112,a123)
    D2V6o[2,2] = d22V6(R2'*Ppt, a111,a112,a123)
    D2V6o[2,3] = d23V6(R2'*Ppt, a111,a112,a123)
    
    D2V6o[3,1] = d31V6(R2'*Ppt, a111,a112,a123)
    D2V6o[3,2] = d32V6(R2'*Ppt, a111,a112,a123)
    D2V6o[3,3] = d33V6(R2'*Ppt, a111,a112,a123)
    
    for a in 1:3, b in 1:3
        
        D2V6[a,b] = 0.0
        
        for c in 1:3, d in 1:3
           
            D2V6[a,b] += R2[a,c]*R2[b,d]*D2V6o[c,d]
            
        end
        
        D2V[a,b] += D2V6[a,b]
        
    end
    
end

eigvecs(D2V)

eigbasisP=xyzfsrt(srtfxyz(sol,R2),eigvecs(D2V))


    return sol, eigbasisP

end

function makefig8(Pplot1,Pplot2, x, N,R2, G2, A2, A4, a111, a112, a123, V0; scale=1.0)

    spcolor = RGBf(0.7,0.7,0.7)



    f = Figure(backgroundcolor = :white, resolution = (1.1*0.5*6.161, 1.5) .* (72*scale), fontsize=8*scale)

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
    
    axislegend(position = :lb, padding=[10,10,5,5],rowgap=0)

    xtks = [-6.0,-3.0,0.0,3.0,6.0]
    ytks = [-1.0,-0.5,0.0,0.5,1.0]
        
    if length(xtks) != 0
        ax1.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax1.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end



    ylims!(ax1, -0.9, 0.9)
    xlims!(ax1, -7, 7)



    ax2 = Axis(f[1,2],
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
        xgridwidth = 0.5*scale,
        ygridwidth = 0.5*scale,
        ygridvisible = false,
        xgridvisible = false

    )
        
   thelabels = [L"\mu_1", L"\mu_2", L"\mu_3" ]
    
    xp = x

    thelinestyles = [:solid,:solid,:solid]
    #thecolors = [:black,:red,:violet]
    thecolors = [ colorant"#d62728", colorant"#9467bd", colorant"#8c564b" ]

    ys = [ Pplot2[a,:] for a in 1:3 ]
    ps = [ lines!(ax2, xp, ys[a], label=thelabels[a], linewidth=1.0*scale, color=thecolors[a]) for a in 1:3 ]
    
    axislegend(position = :lb, padding=[10,10,5,5],rowgap=0)

    xtks = [-6.0,-3.0,0.0,3.0,6.0]
    ytks = [-1.0,-0.5,0.0,0.5,1.0]
        
    if length(xtks) != 0
        ax2.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
    end
    if length(ytks) != 0
        ax2.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
    end



    ylims!(ax2, -1.1, 1.1)
    xlims!(ax2, -7, 7)



    #colgap!(f.layout,1, 10)

    return f
end