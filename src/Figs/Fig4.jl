
function getf4data(N)
	
	counter = 1
	
	ls = 2.0/N
	
	
	
	#xq = -1.0 + ls/2.0 : ls : 1.0 - ls/2.0
	#yq = -1.0 + ls/2.0 : ls : 1.0 - ls/2.0
	#zq = -1.0
	
	Pts = zeros(3, N*N*3)
	
	for i in 1:N, j in 1:N
		
		xq = -1.0 - ls/2.0 + ls*i
		yq = -1.0 - ls/2.0 + ls*j
		zq = -1.0
		
		rq = sqrt( xq^2 + yq^2 + zq^2 )
		
		Pts[:,counter] .= [ xq, yq, zq ]./rq
		counter += 1
			
	end
	
	for i in 1:N, k in 1:N/2
		
		xq = -1.0 - ls/2.0 + ls*i
		zq = -1.0 - ls/2.0 + ls*k
		yq = -1.0
		
		rq = sqrt( xq^2 + yq^2 + zq^2 )
		
		Pts[:,counter] .= [ xq, yq, zq ]./rq
		counter += 1
		
		xq = -1.0 - ls/2.0 + ls*i
		zq = -1.0 - ls/2.0 + ls*k
		yq = 1.0
		
		rq = sqrt( xq^2 + yq^2 + zq^2 )
		
		Pts[:,counter] .= [ xq, yq, zq ]./rq
		counter += 1
		
	end
	
	for j in 1:N, k in 1:N/2
		
		yq = -1.0 - ls/2.0 + ls*j
		zq = -1.0 - ls/2.0 + ls*k
		xq = -1.0
		
		rq = sqrt( xq^2 + yq^2 + zq^2 )
		
		Pts[:,counter] .= [ xq, yq, zq ]./rq
		counter += 1
		
		yq = -1.0 - ls/2.0 + ls*j
		zq = -1.0 - ls/2.0 + ls*k
		xq = 1.0
		
		rq = sqrt( xq^2 + yq^2 + zq^2 )
		
		Pts[:,counter] .= [ xq, yq, zq ]./rq
		counter += 1
		
	end
	
	
	#n=10
	#aspect=(1,1,1)
	#perspectiveness=0.5
	#x, y, z = 
	fig = Figure()#; resolution=(400, 400))
    ax1 = Axis3(fig[1, 1]  )
    scatter!(ax1, Pts[1,:], Pts[2,:], Pts[3,:]; markersize=4)
	display(fig)

	return Pts[:,1:counter-1]
	


end
	
	
	
	
function getPts(N)
	
	counter = 1
	
	ls = 2.0/N
	
	Pts = zeros(3, N*N*3+2*N)
	
	for i in 1:N, j in 1:N
		
		xq = -1.0 - ls/2.0 + ls*i
		yq = -1.0 - ls/2.0 + ls*j
		zq = -1.0
		
		rq = sqrt( xq^2 + yq^2 + zq^2 )
		
		Pts[:,counter] .= [ xq, yq, zq ]./rq
		counter += 1
			
	end
	
	for i in 1:N, k in 1:(N+1)/2
		
		xq = -1.0 - ls/2.0 + ls*i
		zq = -1.0 - ls/2.0 + ls*k
		yq = -1.0
		
		rq = sqrt( xq^2 + yq^2 + zq^2 )
		
		Pts[:,counter] .= [ xq, yq, zq ]./rq
		counter += 1
		
		xq = -1.0 - ls/2.0 + ls*i
		zq = -1.0 - ls/2.0 + ls*k
		yq = 1.0
		
		rq = sqrt( xq^2 + yq^2 + zq^2 )
		
		Pts[:,counter] .= [ xq, yq, zq ]./rq
		counter += 1
		
	end
	
	for j in 1:N, k in 1:(N+1)/2
		
		yq = -1.0 - ls/2.0 + ls*j
		zq = -1.0 - ls/2.0 + ls*k
		xq = -1.0
		
		rq = sqrt( xq^2 + yq^2 + zq^2 )
		
		Pts[:,counter] .= [ xq, yq, zq ]./rq
		counter += 1
		
		yq = -1.0 - ls/2.0 + ls*j
		zq = -1.0 - ls/2.0 + ls*k
		xq = 1.0
		
		rq = sqrt( xq^2 + yq^2 + zq^2 )
		
		Pts[:,counter] .= [ xq, yq, zq ]./rq
		counter += 1
		
	end
	
	
	#n=10
	#aspect=(1,1,1)
	#perspectiveness=0.5
	#x, y, z = 
	fig = Figure()#; resolution=(400, 400))
    ax1 = Axis3(fig[1, 1]  )
    scatter!(ax1, Pts[1,:], Pts[2,:], Pts[3,:]; markersize=4)
	display(fig)

	return Pts[:,1:counter-1]
	


end
	
	
	
	
