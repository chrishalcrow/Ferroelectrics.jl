
function gradient_flow!(P, sums; type=1, shift=2, steps = 1, dt = P.grid.ls/500.0, tolerance = 0.0, checks = max(100,steps), print_stuff = true, dEdp = zeros(3,P.grid.lp) )


    if tolerance == 0 && checks > steps
        checks = steps
    end
    
    if print_stuff == true
        println("initial: energy: ", energy(P,sums,type=type,shift=shift) )
    end

    counter = 0
    prev_error = 1.0e9
    
    while counter < steps
        
        gradient_flow_for_n_steps!(P,dEdp,checks,dt,sums,type)
        
        err = max_abs_err(dEdp)
        if err > 3*prev_error
            error("Suspected numerical blowup. Please use a smaller dt. Currently, dt = ", dt)
        end
        prev_error = err

        counter += checks
        
        if print_stuff == true
            println("after ", counter, " steps, error = ", round(err, sigdigits=4))
        end

        if tolerance != 0.0    # => we are in tol mode    
            if err < tolerance
                counter = steps + 1    # => end the while loop
            else
                steps += checks    # => continue the while loop
            end
        end

		if isnan(err)
			println("error = NaN, choose smaller dt")
			break
		end

    end

    if print_stuff == true
        println("final energy: ", energy(P,sums,type=type,shift=shift) )
    end

    return

end


function gradient_flow_for_n_steps!(ϕ,dEdp,n,dt,sums,type)
    for _ in 1:n
        gradient_flow_1_step!(ϕ,dEdp,dt,sums,type)
    end
end

function gradient_flow_1_step!(P, dEdp, dt,sums,type)

	getdEdP!(P,dEdp,sums,type)
    P.field .-= dt.*dEdp;
   
end


function FastLinearConvolutionGradient(f,g,n)
    if n == 1 && size(f)[1] == length(f) && size(g)[1] == length(g)
        N = length(f)

        f_pad = [f; zeros(N)]   
        g_pad = [g; zeros(1); reverse(g[2:N])]
    else
        N = size(f)[2]
        
        f_pad = [f[1:n,:] zeros(n, N)] 
        g_pad = [g[1:n,:] zeros(n, 1) g[1:n,end:-1:2]]
    end
   
    return FastCyclicConv(f_pad, g_pad, n)
end


function getdEdP!(P,dedp,sums,type)

	G2 = P.parameters.G2
	A2s = P.parameters.A2s
	A4s = P.parameters.A4s
	norm = P.parameters.rescaling

	Threads.@threads for i in 3:P.grid.lp-2

		@inbounds for a in 1:3
			dedp[a,i] = 0.0
		end

		Ppt = getP(P,i)
		DDPpt = getDDP(P,i)
					
		@inbounds for a in 1:3, b in 1:3
		
			dedp[a,i] += -G2[a,b]*DDPpt[b]
			dedp[a,i] += A2s[a,b]*Ppt[b]
			
			for c in 1:3, d in 1:3
				dedp[a,i] += (A4s[a,b,c,d])*Ppt[b]*Ppt[c]*Ppt[d]
			end
			
		end
		
	end
	if P.parameters.is_electrostatic
		N = P.grid.lp
		# Approach #1: capacitor model. P.pars[10] should be abs(g1), energy calculation for the whole sample
		if type == 1
			Threads.@threads for i in 3:N-2
				dedp[1,i] += (Ppt[1])/(8.8541878128e-12*norm)
			end
		elseif type == 2
			#Approach #2: interacting planes infinite in y direction and width b (P.b in code) in z direction. P.pars[10] should be abs(g1)
			for i in 3:N-2
				sum = zeros(3)
				for j in 3:N-2
					sum += [2*P.field[1,j]*sums[1,abs(i-j)+1],      0, 		2*P.field[3,j]*sums[2,abs(i-j)+1]]
				end
				dedp[:,i] += sum*P.b/(4*pi*8.8541878128e-12*norm*P.grid.ls) #in denominator P.ls, because gradient is computed for energy density (J/m^3)
			end
		elseif type == 3
			dedp[1:2:3, 3:N-2] += 2*real(FastLinearConvolutionGradient(P.field[1:2:3,3:N-2],sums[:,1:N-4],2)[:,1:N-4])*P.b/(4*pi*8.8541878128e-12*norm*P.grid.ls)
	
		end
	end
end



function max_abs_err(A)
    return maximum(abs.(A)) 
end






"""
    arrested_newton_flow!(skyrmion; skyrmion_dot, steps = n, tolerance = tol, dt=ls^2/80.0, frequency_of_checking_tolerance = freq, print_stuff = true)
    
Applies an arrested Newton flow to `skyrmion` whose initial time derivative field is skyrmion_dot with timestep `dt`, either for `n` steps or until the error falls below `tol`. The error is checked every `checks` steps.

See also [`gradient_flow!`, `newton_flow!`]
"""
function arrested_newton_flow!(ϕ; ϕd=zeros(3, ϕ.grid.lp), dt=ϕ.grid.ls/50.0, steps=1, tolerance = 0.0, checks = max(100,steps), print_stuff = true)

    if tolerance == 0 && checks > steps
        checks = steps
    end

    energy_density = zeros(ϕ.grid.lp)
    old_field = deepcopy(ϕ.field);

    dEdp = zeros(3, ϕ.grid.lp)
    dEdp2 = zeros(3, ϕ.grid.lp)
    dEdp3 = zeros(3, ϕ.grid.lp)
    dEdp4 = zeros(3, ϕ.grid.lp)
    ϕ2 = deepcopy(ϕ)


    counter = 0
    while counter < steps

        arrested_newton_flow_for_n_steps!(ϕ,ϕ2,ϕd,old_field,dEdp,dEdp2,dEdp3,dEdp4,dt,energy_density,checks, EnergyANF(ϕ,energy_density))
        error = max_abs_err(dEdp)
        counter += checks

        if print_stuff == true 
            println("after ", counter, " steps, error = ", round(error, sigdigits=4), " energy = ", round(sum(energy_density)*ϕ.grid.ls, sigdigits=8) )
            #println( round(error, sigdigits=8), "," )
        end

        if tolerance != 0.0    # => we are in tol mode
            if error < tolerance 
                counter = steps + 1    # => end the while loop
            else
                steps += checks
            end
        end

    end

    return

end
 
function arrested_newton_flow_for_n_steps!(ϕ,ϕ2,ϕd,old_field,dEdp1,dEdp2,dEdp3,dEdp4,dt,energy_density,n, initial_energy)

    new_energy = initial_energy
    
    for _ in 1:n

        old_energy = new_energy
        old_field .= ϕ.field

        newton_flow_for_1_step!(ϕ,ϕ2,ϕd,dEdp1,dEdp2,dEdp3,dEdp4,dt)
        new_energy = EnergyANF(ϕ,energy_density)

		#println("Old, ", old_energy, " and new, ", new_energy)

        if new_energy > old_energy

			println("ARREST")

            fill!(ϕd, 0.0);
            ϕ.field .= old_field;

            if new_energy > 1.2*old_energy
                error("Suspected numerical blow-up. Please use smaller dt. Currently, dt = ", dt)
            end
    
        end

    end

end

function newton_flow_for_1_step!(sk, sk2, skd ,dEdp1, dEdp2, dEdp3, dEdp4, dt)

    getdEdP!(sk, dEdp1)
    sk2.field .= sk.field .+ (0.5*dt).*skd

    #sk.pion_field .+= (0.5*dt).*skd
    #getdEdp!(sk2, dEdp2, sk, dEdp1, dt)
    getdEdP!(sk2, dEdp2)
	sk.field .= sk2.field - (0.5*dt)^2 .*dEdp1

    ##sk.pion_field .+= (0.5*dt)^2 .*dEdp1
    #getdEdp3!(sk, dEdp3, sk2, skd, dEdp1, dEdp2, dt)
    getdEdP!(sk, dEdp3)
	sk2.field .= sk.field .+ (0.5*dt).*skd .- (0.5*dt)^2 .*(4.0 .* dEdp2 .- dEdp1)

    ##sk.pion_field .+= (0.5*dt).*skd .+ (0.5*dt)^2 .*(4.0 .* dEdp2 .- dEdp1)
    #getdEdp4!(sk2, dEdp4, sk, dEdp1, dEdp2, dEdp3, skd, dt)
    getdEdP!(sk2, dEdp4)
	sk.field .= sk2.field + dt.*(  (5/6*dt).*dEdp2  .- (dt/6).*( dEdp1 .+ dEdp3 ) )  

	skd .-= (dt/6.0).*(dEdp1 .+ 2.0.*dEdp2 .+ 2.0.*dEdp3 .+ dEdp4)

    #RESET field to original value: sk.pion_field .-= (0.5*dt).*skd + (0.5*dt)^2 *(4.0 .* dEdp2 - dEdp1) + (0.5*dt)^2 .*dEdp1 + (0.5*dt).*skd
    #Then update: sk.pion_field .+= dt.*(skd + dt/6.0 .*( dEdp1 .+ dEdp2 .+ dEdp3 ) ), combined into:
    ###sk.pion_field .-=  dt.*(  (5/6*dt).*dEdp2  .- (dt/6).*( dEdp1 .+ dEdp3 ) )
    
    #
   
    #orthog_skd_and_sk_and_normer!(skd,sk)
    #normer!(sk)

    ##orthog_skd_and_sk_and_normer!(skd,sk)
   
end
