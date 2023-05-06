function makeString(P1,P2,N,x,y)

    Nx = size(x)[1]
    Ny = size(y)[1]

    Pstring = zeros(N,3,Nx,Ny)

    for a in 1:N, b in 1:3, i in 1:Nx, j in 1:Ny 
        L = (a-1)/(N-1)
        Pstring[a,b,i,j] = (1.0-L)*P1[b,i,j] + L*P2[b,i,j]
    end

    return Pstring

end


function dist(P1,P2)

    return sqrt(sum((P1 - P2).^2))

end

function string_length(Pstring)

    N = size(Pstring)[1]

    length = 0.0
    for a in 1:N-1

        length += dist(Pstring[a+1,:,:,:], Pstring[a,:,:,:])

    end

    return length

end

function string_positions(Pstring)

    N = size(Pstring)[1]

    pos = zeros(N)

    for a in 2:N
        pos[a] = pos[a-1] + dist(Pstring[a,:,:,:], Pstring[a-1,:,:,:])
    end

    return pos

end

function resample_string(Pstring)

    N = size(Pstring)[1]

    current_pos = string_positions(Pstring)
    current_pos ./= current_pos[N]

    desired_pos = 0:1/(N-1):1

    newPstring = zeros(size(Pstring))

    newPstring .= Pstring


    for b in 2:N-1
        p = desired_pos[b]
        #println(p)
        shallweSTOP = false
        a = 1

       #println(current_pos[a])

            while shallweSTOP == false

                if a > N
                    shallweSTOP = true
                end

                if current_pos[a] < p
                    a += 1
                else
                    shallweSTOP = true
                end
            end


            println(p," is between ", current_pos[a-1], " and ", current_pos[a] )

            println( (p - current_pos[a-1])/(current_pos[a]-current_pos[a-1]) )
            println( (current_pos[a]-p)/(current_pos[a]-current_pos[a-1]) )

            newPstring[b,:,:,:] .= (p - current_pos[a-1])/(current_pos[a]-current_pos[a-1]).*Pstring[a,:,:,:] + (current_pos[a]-p)/(current_pos[a]-current_pos[a-1]).*Pstring[a-1,:,:,:]

    end

    return newPstring
end

function get_string_energy(Pstring,G4, A2, A4b,a111,a112,a123,V0,R2,Nx,Ny,dx,dy)

    ED = zeros(Nx,Ny)
    N = size(Pstring)[1]
    engList = zeros(N)

    for a in 1:N
        engList[a] = energy2D(Pstring[a,:,:,:], ED, G4, A2, A4b,a111,a112,a123,0.0,R2,Nx,Ny,dx,dy)
    end

    return engList

end






