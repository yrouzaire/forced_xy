"Copyright (c) 2021 Y.Rouzaire All Rights Reserved."
logspace(x1, x2, n) = [10.0 ^y for y in range(log10(x1), log10(x2), length=n)]
prinz(z) = println("Runtime : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes = $(round(z/3600,digits=2)) hours.")
each = eachindex
nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,x,dims=y)
nanstd(x) = std(filter(!isnan,x))
nanstd(x,y) = mapslices(nanstd,x,dims=y)

## Methods to deal with a square lattice LxL with periodic BCs
function dist(a::Tuple{Int,Int},b::Tuple{Int,Int},L::Int)::Float64  # euclidian distance with Periodic BCs
    dx = abs(a[1] - b[1]) ; dx = min(dx,L-dx)
    dy = abs(a[2] - b[2]) ; dy = min(dy,L-dy)
    return sqrt(dx^2 + dy^2)
end

function dist2(a::Tuple{Int,Int},b::Tuple{Int,Int},L::Int)::Float64  # euclidian distance with Periodic BCs
    dx = abs(a[1] - b[1]) ; dx = min(dx,L-dx)
    dy = abs(a[2] - b[2]) ; dy = min(dy,L-dy)
    return dx^2 + dy^2
end

function dist(a::Tuple{Number,Number},b::Tuple{Number,Number},L::Int)::Float64  # euclidian distance with Periodic BCs
    dx = abs(a[1] - b[1]) ; dx = min(dx,L-dx)
    dy = abs(a[2] - b[2]) ; dy = min(dy,L-dy)
    return sqrt(dx^2 + dy^2)
end

function dist(a::Missing,b::Tuple{Number,Number},L::Int)  # euclidian distance with Periodic BCs
    return missing
end

function manhattan_dist(a::Tuple{Int,Int},b::Tuple{Int,Int},L::Int)::Float64  # euclidian distance with Periodic BCs
    dx = abs(a[1] - b[1]) ; dx = min(dx,L-dx)
    dy = abs(a[2] - b[2]) ; dy = min(dy,L-dy)
    return dx + dy
end

function linear_to_square_index(n::Int,L::Int)
    i = div(n-1,L) + 1  # 1 ≤ i ≤ L
    j = mod1(n,L)       # 1 ≤ j ≤ L
    #= formula for reversing the linear indexing. i is the quotient and j
    the reminder of the euclidian division of n by L and the corrections
    deal with the 1-indexing of Julia =#
    return i,j
end

function square_to_linear_index(i::Int,j::Int,L::Int)::Int
    return L*(i-1) + j # 1 ≤ n ≤ L^2
end

function isneighbour(a,b,L) # a and b are both position Tuples
    if     mod1.(a .+ (+1, 0),L) == b return true # b est a droite de a
    elseif mod1.(a .+ (-1, 0),L) == b return true # b est à gauche de a
    elseif mod1.(a .+ ( 0,-1),L) == b return true # b est en haut à gauche de a
    elseif mod1.(a .+ ( 0,+1),L) == b return true # b est en bas à droite de a
    else return false
    end
end


## Creation and Saving Methods
function initialize(L,T,Var,init,q,r0)
    # Theta field
    if     init == "lowtemp"   thetas = zeros(L,L)
    elseif init == "hightemp"  thetas = rand(L,L)*2π
    elseif init == "isolated"  thetas = create_isolated_vortex(L,L,q)
    elseif init == "pair"      thetas = create_pair_vortices(L,r0)
    elseif init == "2pair"    thetas = create_2pair_vortices(L,r0)
    else error("ERROR : Type of initialisation unknown. Choose among \"hightemp\",\"lowtemp\",\"isolated\" or \"pair\" .")
    end

    omegas = rand(Normal(0,sqrt(Var)),L,L)
    dt = determine_dt(T,Var,L)

    return thetas,omegas,dt
end

function initialize(L,T,Var,init,q,r0,m)
    # Theta field
    if     init == "lowtemp"   thetas = zeros(L,L)
    elseif init == "hightemp"  thetas = rand(L,L)*2π
    elseif init == "isolated"  thetas = create_isolated_vortex(L,L,q)
    elseif init == "pair"      thetas = create_pair_vortices(L,r0)
    elseif init == "2pair"    thetas = create_2pair_vortices(L,r0)
    else error("ERROR : Type of initialisation unknown. Choose among \"hightemp\",\"lowtemp\",\"isolated\" or \"pair\" .")
    end

    omegas = rand(Normal(0,sqrt(Var)),L,L)
    dt = determine_dt(T,Var,m,L)

    return thetas,omegas,dt
end

function determine_dt(T,Var,m,L)
    estimation_max = sqrt(4log(L) * Var) # from math.stackexchange.com/questions/89030, simple answer of Sivaraman
    arbitrary_coeff =  π/10
    dt = min(arbitrary_coeff/estimation_max , arbitrary_coeff/2 , arbitrary_coeff^2*π/4/T , sqrt(m)/5)
    #= A few remarks on the previous line :
        1. Julia can deal with 1/0 = Inf and any number is smaller than Inf
            so no problem with the divisions
        2. The first term (π/10max) is here to ensure that the oscillator
            with the max frequency will be resolved correctly. It grows as 1/sqrt(logN).
        3. The second term should have been π/40 to ensure that the greater
            possible coupling (when each one of all 4 neighbours produces a
            unit coupling). Since this situation only happens at very high T°,
            that it will basically never occur, I chose to relax this upper bound
            to π/20. It's more of a guardrail in case of very small σ AND very
            small T.
        4. The last term is obtained by getting the expectation value of the
            absolute value of the noise ξ, leading to sqrt(π/T)/2. Remember that
            it is sqrt(dt) and not dt as is that enters the Euler Maryuama integrator.
            It is the more constraining term for usual ranges of \sigma and T.
        5. The overall coefficient π/10 is arbitrary and can be modulated as desired.
    =#
    return dt
end

function determine_dt(T,Var,L)
    estimation_max = sqrt(4log(L) * Var) # from math.stackexchange.com/questions/89030, simple answer of Sivaraman
    arbitrary_coeff =  π/10
    dt = min(arbitrary_coeff/estimation_max , arbitrary_coeff/2 , arbitrary_coeff^2*π/4/T)
    return dt
end

function create_isolated_vortex(H,W,q)
    @assert abs(q) == 1
    thetas = zeros(H,W)
    H2 = round(Int,H/2)
    W2 = round(Int,W/2)
    for j in 1:W
        for i in 1:H
            y = H2 - i
            x = j - W2
            thetas[i,j] = q * atan(y,x)
        end
    end
    return thetas
end

function create_2pair_vortices(L,r0)

    return create_pair_vortices(L,r0) + create_pair_vortices(L,min(3r0,Int(L/2)))'
end

function create_pair_vortices(L,r0)
    #= Check for meaningfulness of the defaults separation,
    otherwise the defaults will annihilate with relaxation =#
    @assert r0 ≤ 0.5L  "Error : r0 > L/2. "
    @assert iseven(r0) "Error : r0 has to be even. "

    L2  = round(Int,L/2)
    r02 = floor(Int,r0/2)

    thetas = zeros(L,L)
    for i in 1:L
        for j in 1:L
            thetas[i,j] = atan(j-L2-r02,i-L2) - atan(j-L2+r02,i-L2)
        end
    end
    # Smooth domain walls in order not to create news vortices when relaxing
    thetas[1:3,:] = thetas[end-2:end,:] = zeros(3,L)

    # Let the system relax, while enforcing periodic BCs
    thetas = relax(thetas)

    return thetas
end

function relax(thetas) # Let the system relax at T = Var = 0
    L = size(thetas)[1]

    t = 0.0 ; dt = determine_dt(0,0,L) ; t_relax = 2 # t_relax = 2 is enough to smoothen the manually imposed PBCs
    while t<t_relax
        thetas_new = thetas
        for j in 1:L
            for i in 1:L
                angle_neighbours = get_angle_neighbours(thetas,i,j,L) # O(1)
                θ = thetas[i,j]
                sin_angle_neighbours = sin.(angle_neighbours .- θ)
                thetas_new[i,j] =  θ + dt*sum(sin_angle_neighbours)
            end
        end
        t += dt
        thetas = thetas_new
    end
    return thetas
end

## Get info about nearest spins
function get_angle_neighbours(thetas::Matrix{Float64},i::Int,j::Int,L::Int)::Vector{Float64}
    return [thetas[i,mod1(j+1,L)],
            thetas[mod1(i-1,L),j],
            thetas[i,mod1(j-1,L)],
            thetas[mod1(i+1,L),j]] ## Note : mod1 instead of mod to deal with the indexing at from 1
end

function get_angle_neighbours_borders(thetas::Matrix{Float64},i::Int,j::Int,L::Int)::Vector{Float64}
    neighbours = []
    if     i == 1 && j == 1 return [thetas[1,2],thetas[2,1]]     # topleft corner
    elseif i == 1 && j == L return [thetas[1,L-1],thetas[2,L]]   # topright corner
    elseif i == L && j == 1 return [thetas[L-1,1],thetas[L,2]]   # bottomleft corner
    elseif i == L && j == L return [thetas[L,L-1],thetas[L-1,L]] # bottomright corner

    elseif i == 1 return [thetas[1,j-1],thetas[1,j+1],thetas[2,j]]   # top
    elseif i == L return [thetas[L,j-1],thetas[L,j+1],thetas[L-1,j]] # bottom
    elseif j == 1 return [thetas[i-1,1],thetas[i+1,1],thetas[i,2]]   # left
    elseif j == L return [thetas[i-1,L],thetas[i+1,L],thetas[i,L-1]] # right

    end
end

function get_omega_nearest_neighbours(omegas::Matrix{Float64},i::Int,j::Int,L::Int)::Vector{Float64}
    return [omegas[i,mod1(j+1,L)],
            omegas[mod1(i-1,L),j],
            omegas[i,mod1(j-1,L)],
            omegas[mod1(i+1,L),j]] ## Note : mod1 instead of mod to deal with the indexing at from 1
    end

## "Get" methods
function get_Delta(thetas::Matrix{Float64})::Vector{Float64}
    L = size(thetas)[1]
    Δ = zero(ComplexF64)
    Δ_proj = 0.0 # Δ projected on x-axis
    for theta in thetas # scan all spin values, no matter the order
        Δ      += exp(im*theta) # im is the imaginary unit : im² = -1
        Δ_proj += cos(theta)
    end
    Δ = abs(Δ)/L^2
    Δ_proj = Δ_proj/L^2
    return [Δ,Δ_proj]
end


function get_delta(thetas::Matrix{Float64},L::Int,l::Int)::Matrix{Float64}
    # L is the size of the system
    # l is the size of the local box for averaging the magnetisation
    @assert isinteger(L/l)
    c = Int(L/l)
    δ = zeros(2,c^2)
    for j in 1:c
        for i in 1:c
            δ[:,square_to_linear_index(i,j,c)] = get_Delta(thetas[1+l*(i-1):i*l, 1+l*(j-1):j*l])
        end
    end
    return δ
end

function get_C_perp(thetas::Matrix{Float64},i::Int,j::Int,n::Int,L::Int)::Float64 # used in get_C below
    angle_neighbours_at_distance_n = [thetas[mod1(i+n,L),j], ## Note : mod1 instead of mod to deal with the indexing at from 1
                        thetas[mod1(i-n,L),j],
                        thetas[i,mod1(j+n,L)],
                        thetas[i,mod1(j-n,L)] ]
    return mean(Float64[cos(angle - thetas[i,j]) for angle in angle_neighbours_at_distance_n])
end

function get_C_diag(thetas::Matrix{Float64},i::Int,j::Int,n::Int,L::Int)::Float64 # used in get_C below
    angle_neighbours_at_distance_sqrt2_n = [thetas[mod1(i+n,L),mod1(j+n,L)], ## Note : mod1 instead of mod to deal with the indexing at from 1
                        thetas[mod1(i+n,L),mod1(j-n,L)],
                        thetas[mod1(i-n,L),mod1(j+n,L)],
                        thetas[mod1(i-n,L),mod1(j-n,L)] ]
    return mean(Float64[cos(angle - thetas[i,j]) for angle in angle_neighbours_at_distance_sqrt2_n])
end

# function get_C(thetas::Matrix{Float64})::Vector{Float64} # returns C(r) for angles
#     L = size(thetas)[1] ; Lover2 = round(Int,L/2,RoundDown)
#     matrix = Matrix{Float64}(undef,L,L^2)
#
#     # r_vector = vcat([n for n in 1:round(Int,L/2,RoundDown)],[n*sqrt(2) for n in 1:round(Int,L/2,RoundDown)])
#
#     # We fill the matrix with perpendicular directions first (vertical and horizontal)
#     for i in 1:L
#         for j in 1:L
#             for n in 1:Lover2
#                 matrix[n,L*(i-1) + j] = get_C_perp(thetas,i,j,n,L)
#                 # (i,j) -> L*(i-1) + j in linear indexing
#             end
#         end
#     end
#
#     # We then fill the matrix with diagonal directions
#     for i in 1:L
#         for j in 1:L
#             for n in 1:Lover2
#                 matrix[n+Lover2,L*(i-1) + j] = get_C_diag(thetas,i,j,n,L)
#             end
#         end
#     end
#
#     result_unsorted = Statistics.mean(matrix,dims=2) # average over spins = over columns
#         # discretization of the distance r, in 2 parts for simplicity
#         # in the matrix construction
#     r_vector = vcat([n for n in 1:Lover2],[n*sqrt(2) for n in 1:Lover2])
#     sortpermm = sortperm(r_vector)
#     result_sorted = result_unsorted[sortpermm]
#
#     return result_sorted # finally, one gets C(t,r)
#
#     #= Note 1 : Since the work needed to fill each cell of the Matrix "matrix"
#     is consequent, I believe the ordering row vs column will not be limiting
#     in this case. =#
#
#     #= Note 2 : This part of the code is of complexity O(L^3). It is already
#     one order of complexity above the other methods, which are usually O(L^2).
#     This is why I chose not to implement other approximations for the euclidian
#     distance r, which would have pleasantly refined its discretization.
#     There are a lot of other possibilities : i.e. the displacement of a cavalier
#     in a chess game. But already there, the combinatorics begin to increase
#     (8 possibilities). Implementing more and more of these would let the complexity
#     go up to  O(L^4). =#
#
#     #= Note 3 : One may remark that the output is an average of averages. In this
#     particular case, it is totally correct because all the sub-samples are of the
#     same size. Note that this is another argument is favour of not computing other
#     types of discretized distances. =#
# end

function get_C(thetas::Matrix{Float64})::Vector{Float64} # returns C(r) for angles
    L = size(thetas)[1] ; Lover2 = round(Int,L/2,RoundDown)
    matrix = Matrix{Float64}(undef,Lover2,L^2)

    # We fill the matrix with perpendicular directions first (vertical and horizontal)
    for i in 1:L
        for j in 1:L
            for n in 1:Lover2
                matrix[n,L*(i-1) + j] = get_C_perp(thetas,i,j,n,L)
                # (i,j) -> L*(i-1) + j in linear indexing
            end
        end
    end
    return Statistics.mean(matrix,dims=2)[:,1]
end

function C_all_directions(lattice) # to make sure there is no anisotropy
    L = size(lattice)[1] ; L2 = round(Int,L/2)
    C_vertical = zeros(L2,L^2)
    C_horizontal = zeros(L2,L^2)
    C_diag = zeros(L2,L^2)
    C_anti_diag = zeros(L2,L^2)
    for n in 1:1:L^2
        i,j = linear_to_square_index(n,L)
        angle = lattice[i,j,1]
        for k in 1:L2
            C_vertical[k,n]   = cos(angle - lattice[mod1(i+k,L),j,1])
            C_horizontal[k,n] = cos(angle - lattice[i,mod1(j+k,L),1])
            C_diag[k,n]       = cos(angle - lattice[mod1(i+k,L),mod1(j+k,L),1])
            C_anti_diag[k,n]  = cos(angle - lattice[mod1(i-k,L),mod1(j+k,L),1])
        end
    end
    C_vertical = mean(C_vertical,dims=2)
    C_horizontal = mean(C_horizontal,dims=2)
    C_diag = mean(C_diag,dims=2)
    C_anti_diag = mean(C_anti_diag,dims=2)

    return C_vertical , C_horizontal , C_diag , C_anti_diag
end


function get_Ctw(thetas_t::Matrix{Float64},t::Float64,token_time_t::Int,lattices_tw::Array{Matrix{Float32},1},tws::Vector{Int64},Ctw::Matrix{Float64}) ## remplit une colonne de la matrice Ctw
    for i in eachindex(tws)
        tw = tws[i]
        if t > tw # on peut remplir la colonne 'token_time_t' de la matrice Ctw
            Ctw[i,token_time_t] = mean(cos.(thetas_t .- lattices_tw[i])) # returns ⟨∑ cos(θ_i(t) - θ_i(t_w))⟩/L²
        else
            #= On ne peut pas encore, donc on écrit des NaN.
            Missing aurait été plus approprié mais ce n'est pas un Float, et on obtient
            l'erreur suivante : ArgumentError: type of SharedArray elements must
            be bits types, got Union{Missing, Float64} =#
            Ctw[i,token_time_t] = NaN
        end
    end
    return nothing # one writes the measurement in place, in the matrix Ctw
end

function get_fastest_oscillators(omegas::Matrix{Float64} ; threshold::Number = NaN , number::Int = 0 )
    threshold_is_modified = threshold ≠ NaN && threshold > 0
    number_is_modified    = threshold ≠ 0   && number    > 0
    @assert xor(threshold_is_modified,number_is_modified) "Either the entered value is incorrect or both keyargs were used at the same time."
    L = size(omegas)[1]
    if threshold_is_modified
        result = Tuple{Int,Int}[]
        for j in 1:L
            for i in 1:L
                if abs.(omegas[i,j]) > threshold push!(result,(i,j)) end
            end
        end
        return result
    elseif number_is_modified
        omegas_abs = abs.(omegas)
        return linear_to_square_index.(Array(partialsortperm(vec(omegas_abs), 1:number, rev=true)),L)
    end
end

## Gestion des Runaways
function spot_runaways(lattice_omega::Matrix{Float64},threshold::Number,nbroken::Int)::Vector{Tuple{Int,Int}} # return their location (i,j) in a Array of Tuple
    runaway_locations = Tuple[]
    L = size(lattice_omega)[1]

    for i in 1:L
        for j in 1:L
            omega_spin = lattice_omega[i,j]
            omega_neighbours = get_omega_nearest_neighbours(lattice_omega,i,j,L)
            Δω = abs.(omega_neighbours .- omega_spin)

            # if sum(Δω .> threshold) ≥ nbroken
            if sum(Δω .> threshold) == nbroken
                push!(runaway_locations,(i,j))
            end
        end
    end
    return runaway_locations
end

function create_lattice_runaway(L,T,Var,init,threshold,nbroken)::Tuple{Array{Float64,3},Float64,Int,Int} # return lattice,dt,i_runaway,j_runaway
    if Var == 0 error("σ = 0 , no runaway configuration is possible !") end
    ntries_max = 10000 ; ntries = 1
    while ntries < ntries_max
        try
        # Create a lattice
        lattice,dt = create_Lattice(L,T,Var,init)
        current_nbroken = nbroken # the nbroken at which one currently works (will decrease in the SEARCH PHASE and increase in the SWAP PHASE)

        ## SEARCH PHASE : Find the weakest spin (the one that has the most broken links)
        iweak,jweak = -1,-1 # will contain the location of the runaway
        omegaweak = NaN      # will contain the \omega of the runaway
        list_runaways = spot_runaways(lattice[:,:,2],threshold,current_nbroken) # array of tuples (i,j)
        if length(list_runaways) > 0
            iweak,jweak = list_runaways[1][1],list_runaways[1][2]
            return lattice,dt,iweak,jweak # because the lattice has already what you want, i.e at least one runaway oscillator (thus return the first element of the list)
        else # there is currently no runaways with the demanded threshold/nbroken
            bool_found = false
            current_nbroken = current_nbroken - 1 # now we look for runaways in the unchanged lattice but with a looser constraint : nbroken-1
            while current_nbroken > 0
                list_runaways = spot_runaways(lattice[:,:,2],threshold,current_nbroken)
                if length(list_runaways) > 0
                    #= we have a list of tuples of weak spins. We want the most fragile,
                    i.e. the one that has the largest intrinsic frequency difference
                    with its neighbours involved in currently existing links =#
                    tmp = 0 # will contain the index of the weakest runaway found so far
                    value_tmp = -1 # will contain the biggest frequency difference found so far
                    for i in eachindex(list_runaways) # scan over all runaways
                        # for each runaway, get the sum of the frequency difference with their existing neighbours
                        omega_neighbours = [ lattice[neighbour[1],neighbour[2],2] for neighbour in get_existing_neighbours(lattice,list_runaways[i][1],list_runaways[i][2],threshold) ]
                        diff_omega_neighbours = abs.( omega_neighbours .- lattice[list_runaways[i][1],list_runaways[i][2],2] )
                        sum_diff_frequency = sum(diff_omega_neighbours)
                        if sum_diff_frequency > value_tmp
                            tmp = i
                        end
                    end
                    # We now have the most fragile of the spins :
                    iweak,jweak = list_runaways[tmp][1] , list_runaways[tmp][2]
                    omegaweak   = lattice[iweak,jweak,2]
                    bool_found  = true
                    break # we found the weakest spin -> phase 1 (SEARCH) is thus over. One goes directly to PHASE 2 : the swap procedure
                else
                    current_nbroken = current_nbroken - 1 # now we look again for runaways in the unchanged lattice but with a looser constraint : nbroken-1
                end # if
            end # while current_nbroken > 0
            if !bool_found error("No runaway with nbroken ≥ 1 ! Try again.") end # If one gets here, it means that there is no runaway with nbroken ≥ 1
        end # if length() > 0

            ## SWAP PHASE
            #= at this stage, "weak_spin" has only current_nbroken links broken.
            We aim for nbroken links. Thus, we are going to swap spins involved
            in existing links to weaken a link so that it breaks. Therefore, the
            number of actually broken links will now increase. =#
            broken_neighbours_list = get_broken_neighbours(lattice,iweak,jweak,threshold) # list of tuples (i,j) of neighbours of the weak spin that are involved in a broken link
            while current_nbroken < nbroken
                i_existing_neighbour,j_existing_neighbour = get_existing_neighbours(lattice,iweak,jweak,threshold)[1] # one gets the/one neighbour (of our weak spin) that is involved in a still existing link

                # Let's find the spin most distant to the weak spin's omega
                mdsf = (0,0,0) # will contain both the coordinates and Δω of the Most Distant (wrt ω) So Far (=mdsf) to our weak spin throughout the search : (i,j,Δω)
                for i in 1:L
                    for j in 1:L
                        Δω = abs(lattice[i,j,2] - omegaweak)
                        if Δω > mdsf[3] && (i,j) ∉ broken_neighbours_list # check whether its not one of the already broken neighbours
                            mdsf = (i,j,Δω) # update the Most Distant spin So Far in our search
                        end
                    end
                end

                # Check whether a swap would break the link. If so, proceed to the swap. If not, throw an error and abort mission.
                if abs(omegaweak - lattice[mdsf[1],mdsf[2],2]) ≥ threshold # a swap would break the new link
                    # swap both spins
                        tmp = lattice[mdsf[1],mdsf[2],:]
                        lattice[mdsf[1],mdsf[2],:] = lattice[i_existing_neighbour,j_existing_neighbour,:]
                        lattice[i_existing_neighbour,j_existing_neighbour,:] = tmp
                    current_nbroken = current_nbroken + 1 # the swap has broken the link ...
                    push!(broken_neighbours_list,(i_existing_neighbour,j_existing_neighbour)) # ... so one has to add it to the list of broken neighbours
                    # from there, the while loop will continue and if the required number of broken neighbours is met, job is done. Otherwise, loop again
                else
                    # even a swap with the furthest spin would not break the link, our quest is hopeless
                    error("There is no spin in this lattice distant enough to meet the requierements in terms of threshold/nbroken. Try again")
                end
            end # While current_nbroken < nbroken

            # if you made it until here, it's because the modified lattice now has what you want : a runaway oscillator at location (iweak,jweak)
            return lattice,dt,iweak,jweak
        catch
            ntries = ntries + 1
        end # end try
    end # while try delivers an error
    error("After $ntries_max attempts, it has been impossible to create the requiered lattice.")
end

function get_broken_neighbours(lattice,iweak,jweak,threshold)::Vector{Tuple{Int,Int}} # returns a list of tuples (i,j) of neighbours of the weak spin that are involved in a broken link
    L = size(lattice)[1]
    list = Tuple[]
    neighbours = get_neighbours(L,iweak,jweak) # neighbours is a list of tuples (i,j)
    for neighbour in neighbours
        i_neighbour = neighbour[1] ; j_neighbour = neighbour[2]
        if abs(lattice[iweak,jweak,2] - lattice[i_neighbour,j_neighbour,2]) ≥ threshold # that link is broken
            push!(list,neighbour)
        end
    end
    # return is hidden in next line
    if length(list) > 0 return list else error("No broken neighbours. Something must have gone wrong.") end
end

function get_existing_neighbours(lattice,iweak,jweak,threshold)::Vector{Tuple{Int,Int}}  # returns all neighbours of the weak spin that are involved in an existing link.
    L = size(lattice)[1]
    list = Tuple[]
    neighbours = get_neighbours(L,iweak,jweak) # neighbours is a list of tuples (i,j)
    for neighbour in neighbours
        i_neighbour = neighbour[1] ; j_neighbour = neighbour[2]
        if abs(lattice[iweak,jweak,2] - lattice[i_neighbour,j_neighbour,2]) < threshold # that link does exists
            push!(list,neighbour)
        end
    end

    # return is hidden in next line
    if length(list) > 0 return list else error("No effective (contrary of 'broken') neighbours. Something must have gone wrong.") end
end

## Gestion des Vortex
function arclength(theta1::Float64,theta2::Float64)::Float64
    #= This function returns the signed arclength on the unit trigonometric circle .
    Clockwise        >> sign -
    Counterclockwise >> sign +
    Note that the inputs thetas are within [0,2π] =#
    dtheta = theta2 - theta1
    dtheta_abs = abs(theta2 - theta1)

    shortest_unsigned_arclength = min(2π-dtheta_abs,dtheta_abs)
    if dtheta_abs < π
        signe = sign(dtheta)
    else
        signe = -sign(dtheta)
    end
    return signe*shortest_unsigned_arclength
end

function get_vorticity(thetas::Matrix{Float64},i::Int,j::Int,L::Int)::Int
    #= Note : here, i and j are the coordinates of the plaquette one considers.
        The top-left corner spin of a plaquette (i,j) also has coordinates (i,j).
        By convention, a plaquette has an orientation : we start with the topleft
        corner and rotate in the Counterclockwise direction, hence the particular
        ordering of the spins in the following line of code
        1---4
        |   |
        2---3
        =#
    angles_corners = mod.([thetas[i,j],thetas[mod1(i+1,L),j],thetas[mod1(i+1,L),mod1(j+1,L)],thetas[i,mod1(j+1,L)]],2π)
    perimeter_covered  = arclength(angles_corners[1],angles_corners[2])
    perimeter_covered += arclength(angles_corners[2],angles_corners[3])
    perimeter_covered += arclength(angles_corners[3],angles_corners[4])
    perimeter_covered += arclength(angles_corners[4],angles_corners[1])
    # if isnan(perimeter_covered) println("i $i j $j lattice $lattice_angles") end
    return round(Int,perimeter_covered/2π)
end

function spot_defects(thetas::Matrix{Float64},BC)::Vector{Tuple{Float64,Float64,Float64}}
    L = size(thetas)[1]
    if BC == "free" range_bc = 2:L-1
    elseif BC == "periodic" range_bc = 1:L
    else error("Boundary conditions not known. Should be either 'periodic' or 'free'.")
    end
    list_vortices = Tuple{Float64,Float64,Float64}[] # will contain a tuple for each defect : (i,j,q), where q is its charge
    for i in range_bc
        for j in range_bc
            q = get_vorticity(thetas,i,j,L)
            if abs(q) > 0.1 # we want to keep ±½ defects, and not rounding errors
                push!(list_vortices,(i,j,round(q,digits=2)))
                # store rounded values because otherwise a q = 0.499 is considered different from a q = 0.5
            end
        end
    end
    #= In this list, there might be doubles/triples (2/3 locations for the
    same physical vortex). We thus seek for numerically identified vortices
    which are neighbours and with the same charge to delete them. =#
    # elements_to_keep = BitVector(ones(length(list_vortices)))
    # for i in 2:length(list_vortices)
    #     for j in 1:i-1
    #         if list_vortices[i][3] == list_vortices[j][3] && isneighbour(list_vortices[i][1:2],list_vortices[j][1:2],L)
    #             # if same charge and isneighbour(vortex_i,vortex_j) == true (for at least one j), we have a double, so we don't keep it
    #             elements_to_keep[i] = 0
    #             break # will only break the inner for loop and continue the outer one
    #         end
    #     end
    # end
    return list_vortices
    # return list_vortices[elements_to_keep]
end

function spot_defects_separated(thetas::Matrix{Float64},BC)
    L = size(thetas)[1]
    if BC == "free" range_bc = 2:L-1
    elseif BC == "periodic" range_bc = 1:L
    else error("Boundary conditions not known. Should be either 'periodic' or 'free'.")
    end
    list_vortices     = Tuple{Int,Int}[] # will contain a tuple for each defect : (i,j)
    list_antivortices = Tuple{Int,Int}[] # will contain a tuple for each defect : (i,j)
    for i in range_bc
        for j in range_bc
            q = get_vorticity(thetas,i,j,L)
            if     q > +0.1 push!(list_vortices,(i,j))
            elseif q < -0.1 push!(list_antivortices,(i,j))
            end
        end
    end
    return list_vortices,list_antivortices
end

function dist_pairs(u::Vector{Tuple{Int,Int}},v::Vector{Tuple{Int,Int}},L)
    m = length(u)
    distances = Float64[]
    if u == v
        for i in 1:m , j in i+1:m  push!(distances,dist(u[i],u[j],L)) end
    else
        for i in 1:m , j in 1:m    push!(distances,dist(u[i],v[j],L)) end
    end
    return filter!(x->x≤L/2 && x>0,distances)
end

function grt(u::Vector{Tuple{Int,Int}},v::Vector{Tuple{Int,Int}},L)
    distances = dist_pairs(u,v)
    n = length(distances)
    h = fit(Histogram,distances,0:0.25:floor(Int,L/2))
    r = Array(h.edges[1]) ; dr = r[2]-r[1]
    g = h.weights/π/dr ./ r[2:end]
    g_gas = length(distances)/π/(L/2)^2

    return g/g_gas
end


# function spot_vortices(thetas::Matrix{Float64})::Vector{Tuple{Int,Int,Int}}
#     L = size(thetas)[1]
#     list_vortices = Tuple{Int,Int,Int}[] # will contain a tuple for each vortex : (i,j,q), where q is its charge
#     for i in 1:L
#         for j in 1:L
#             q = get_vorticity(thetas,i,j,L)
#             if q ≠ 0
#                 push!(list_vortices,(i,j,q))
#             end
#         end
#     end
#     # @assert sum_q == 0
#     return list_vortices
# end

function classify_vortices(thetas::Matrix{Float64})
    #= Let's sum up what this function does. Characterising a vortex pair as free
    or bounded is not a trivial question, be it numerically or conceptually.
    There currently exists no precise definition in the litterature and what follows
    does not pretend to be a definitive one. The underlying idea is to OPTIMALLY
    pair the vortices with their closest opposite-polarity neighbour, optimal in
    the sense that the total distance between pairs should be minimized.
    Then, we define a pair as free if its distance is above a given threshold,
    and as bounded if below that threshold. If one decides to reason in terms of
    distances to classify pairs of vortices, I reckon this optimal bipartite matching
    is the only valid solution. Any of the other simpler procedures I've tried led to
    incoherences, at least when visually checking the results.

    To carry out the bipartite optimal matching, we use the so-called Hungarian Matching algoritm,
    also called Munkres' algoritm after James Munkres, the first to have proved that the
    complexity of this algo was not actually exponential as thought until then, but rather
    polynomial : O(m^3) where m is the number of vortex pairs. (for the record, he showed
    it was O(m^4) and the bound was lowered in the 70's).
    Note 1 : m = O(L²) is the worst-case (high T and/or high σ²), leading the algorithm
        to be O(L^6) in the worst-case. =#

    L = size(thetas)[1]

    list_vortex,list_antivortex  = spot_vortex_antivortex(thetas) # Vector{Tuple{Int,Int}}
    m = length(list_vortex) # number of pairs vortex/antivortex

    # Labels
        # label = 0  : free vortex
        # label ≥ 1  : bounded vortex, each pair has a unique ID (the token) ≥ 1
    label_vortex = fill(-1,m) ; label_antivortex = fill(-1,m)
    token = 1

    #= Now comes the tricky question of what distance should we choose as
    threshold to determine whether a vortex pair is bounded or not. Let's dare the
    crude approximation that the pairs of vortices are uniformly distributed on
    the square.
    Let's consider the behaviour only at large L (actually, I reckon that this second
    approximation should be more controlled in a periodic lattice).
    Based on the first answer by Douglas Zare on mathoverflow.net/questions/124579
    the average minimal distance between each one of the m/2 pairs is of the order of 2/m.
    Since if the threshold is too large, one will never end up with free vortices
    classified as such, let's divide for security this threshold by a factor two,
    leading to threshold_dist = 1/m (here for a unit square -> L/m for a LxL lattice.)

    Note 2 : the threshold distance decreases exponentially with T and/or σ² since
        the number of vortices increases exponentially.
    Note 3 : since m = O(L²) is the worst-case, one has to manually ensure that
        the threshold does not decrease too much. Here comes in the first arbitrary
        parameter, namely that constant in the max. Indeed, we observe in the simulations
        that vortices that one would intuitively consider as bounded can very well diffuse
        (on distance of the order of a few lattice-spacings) due to random thermal
        fluctuations.
        The value 3 is justified thanks to the benchmarks, showing that for ≤ 2 and ≥ 5
        it produces poor results (L=200 , R=30 , tmax=500) =#
    threshold_dist = max(3,L/m)

    dist_matrix = Array{Float64,2}(undef,m,m)
    for j in 1:m # j is an antivortex
        for i in 1:m # i is a vortex
            dist_matrix[i,j] = dist(list_vortex[i],list_antivortex[j],L)
        end
    end

    try # if there is no vortex at all, hungarian rises "ArgumentError: reducing over an empty collection is not allowed"
        matching,~ = Hungarian.hungarian(dist_matrix) # here is done the optimal bipartite matching

        for i in 1:m # i is a vortex
            j = matching[i] # j is an antivortex
            if dist_matrix[i,j] ≤ threshold_dist # you just found its pair, so label it with a unique token !
                label_vortex[i] = token
                label_antivortex[j] = token
                token += 1
            else # its closest pair is further than the threshold -> free vortex
                label_vortex[i] = label_antivortex[j] = 0
            end
        end
        number_free_vortices  = 2sum(label_vortex .== 0)
        number_all_vortices   = 2m

        return label_vortex,label_antivortex,number_all_vortices,number_free_vortices
    catch e
        # println("Caught Error ",sprint(showerror, e))
        return Tuple{Int,Int}[],Tuple{Int,Int}[],0,0  # if there was no vortex at all in the first place, output this
    end
end


function initialize_fiches_defects(thetas,BC,t)
    fiche_vortices     = Vector{Vector{Any}}()
    fiche_antivortices = Vector{Vector{Any}}()
    #= Explanation of the data structure
        - each fiche_vortices[i] is the data corresponding to a vortex
        - this data is a vector composed of
            1 * [the vortex ID ,
            2 * the last known location (i,j)
            3 * the creation time ,
            4 * the location of the creation
            5 * the annihilation time (nothing until known) ,
            6 * the location of the annihilation (nothing until known)
            7 * the ID of the antivortex with which it annihilated (nothing until known)]
    =#

    vortices_old,antivortices_old = spot_defects_separated(thetas,BC)
    m_old = length(vortices_old)
    for i in 1:m_old
        push!(fiche_vortices,[i,vortices_old[i],t,vortices_old[i],nothing,nothing,nothing])
        push!(fiche_antivortices,[i,antivortices_old[i],t,antivortices_old[i],nothing,nothing,nothing])
    end

        return fiche_vortices,fiche_antivortices,vortices_old,antivortices_old
end

function location_creation(fiche)
    locations = Vector{Tuple{Int,Int}}(undef,length(fiche))
    for i in eachindex(fiche)
        locations[i] = fiche[i][4]
    end
    return locations
end

function location_annihiliation(fiche)
    locations = Tuple{Float64,Float64}[]
    for i in eachindex(fiche)
        tmp = fiche[i][6]
        if tmp ≠ nothing push!(locations,tmp) end
    end
    return locations
end

# function defect_trackerv1(thetas::Matrix{Float64},BC,fiche_vortices::Vector{Vector{Any}},fiche_antivortices::Vector{Vector{Any}},vortices_old::Vector{Tuple{Int,Int}},antivortices_old::Vector{Tuple{Int,Int}},t::Float64)
#     vortices_new,antivortices_new = spot_defects_separated(thetas,BC)
#     m_new,m_old = length(vortices_new),length(vortices_old)
#
#     # Special simple cases to deal with upstream : creation of a pair when there was no defects before
#     if m_new == m_old == 0
#         # do nothing, otherwise, "reducing over empty collection blablabla"
#     elseif m_new == 1 && m_old == 0
#         push!(fiche_vortices,[length(fiche_vortices)+1,vortices_new[1],t,nothing,nothing])
#         push!(fiche_antivortices,[length(fiche_antivortices)+1,antivortices_new[1],t,nothing,nothing])
#     elseif m_new == 0 && m_old > 0
#         for i in eachindex(fiche_vortices)
#             # ici je triche, je laisse par flemme l'ID de l'autre defaut = nothing
#             if fiche_vortices[i][4] == nothing fiche_vortices[i][4] = t end
#             if fiche_antivortices[i][4] == nothing fiche_antivortices[i][4] = t end
#         end
#     elseif abs(m_new - m_old) ≤ 1 # allows for at most one creation or annihilation process
#         distance_matrix_vortices = get_distance_matrix(vortices_new,vortices_old,size(thetas,1)) # m_new lignes, m_old colonnes
#
#         assignment_vortices = proposal_vortices = hungarian(distance_matrix_vortices)[1] # length = m_new
#         distance_matrix_antivortices = get_distance_matrix(antivortices_new,antivortices_old,size(thetas,1))
#         assignment_antivortices = proposal_antivortices = hungarian(distance_matrix_antivortices)[1] # length = m_new
#         #= NB : the proposal is the matching between indices of the vectors
#         vortices_new,vortices_old while the assignment matches actual IDs. =#
#         for i in eachindex(assignment_vortices) , j in eachindex(fiche_vortices)
#             if assignment_vortices[i] ≠ 0 && fiche_vortices[j][4] == nothing && fiche_vortices[j][2] == vortices_old[proposal_vortices[i]]
#                 # the first condition is a protection against the creation case, where the Hungarian algo matches the newly created vortex to 0
#                 # the second condition ensures that one only considers currently living vortices and not vortices now annihilated
#                 assignment_vortices[i] = j
#                 break
#             end
#         end
#         for i in eachindex(assignment_antivortices) , j in eachindex(fiche_antivortices)
#             if assignment_antivortices[i] ≠ 0 && fiche_antivortices[j][4] == nothing && fiche_antivortices[j][2] == antivortices_old[proposal_antivortices[i]]
#                 assignment_antivortices[i] = j
#                 break
#             end
#         end
#
#         #= =#
#
#
#     # Case 1 : no creation, no annihilation : simply update the data structure
#         if m_new == m_old
#             for i in 1:m_new # scan assignment
#                 fiche_vortices[assignment_vortices[i]][2] = vortices_new[i]
#                 fiche_antivortices[assignment_antivortices[i]][2] = antivortices_new[i]
#             end
#
#     # Case 2 : creation !
#         elseif m_new > m_old
#             #= locate the newly created defects :
#             the assignment vector contains a 0, corresponding to the newly created vortex =#
#             ind_created_vortex = findfirst(iszero,assignment_vortices)
#             loc_created_vortex = vortices_new[ind_created_vortex]
#             ind_created_antivortex = findfirst(iszero,assignment_antivortices)
#             loc_created_antivortex = antivortices_new[ind_created_antivortex]
#
#             # push them to the data structure
#             id_new = length(fiche_vortices) + 1
#             push!(fiche_vortices,[id_new,loc_created_vortex,t,nothing,nothing])
#             push!(fiche_antivortices,[id_new,loc_created_antivortex,t,nothing,nothing])
#
#             # erase them from the assignment vector and from vortices_new (unused after)
#             deleteat!(assignment_vortices,ind_created_vortex)
#             deleteat!(vortices_new,ind_created_vortex)
#             deleteat!(assignment_antivortices,ind_created_antivortex)
#             deleteat!(antivortices_new,ind_created_antivortex)
#
#             # update via the remaining assignment vector
#             for i in 1:m_new-1 # scan assignment modified, hence the -1
#                 fiche_vortices[assignment_vortices[i]][2] = vortices_new[i]
#                 fiche_antivortices[assignment_antivortices[i]][2] = antivortices_new[i]
#             end
#
#     # Case 3 : annihilation !
#         elseif m_new < m_old
#              # the annihilated vortex is absent from the assignment vector, so let's proceed
#             for i in 1:m_new # scan assignment
#                 fiche_vortices[assignment_vortices[i]][2] = vortices_new[i]
#                 fiche_antivortices[assignment_antivortices[i]][2] = antivortices_new[i]
#             end
#
#             # now deal with the annihilation stuff
#             # first, find the missing index in the assignment vector
#             annihilated_vortex = -1 ; annihilated_antivortex = -1
#             for i in eachindex(fiche_vortices)
#                 if (i ∉ assignment_vortices ) && (fiche_vortices[i][4] == nothing) # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
#                     annihilated_vortex = i
#                 end
#             end
#             for i in eachindex(fiche_antivortices)
#                 if (i ∉ assignment_antivortices ) && (fiche_antivortices[i][4] == nothing) # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
#                     annihilated_antivortex = i
#                 end
#             end
#             fiche_vortices[annihilated_vortex][4] = t # fiche_vortices[i][4] is the time of annihilation
#             fiche_vortices[annihilated_vortex][5] = annihilated_antivortex # fiche_vortices[i][5] is the ID of the colliding antivortex
#
#             fiche_antivortices[annihilated_antivortex][4] = t
#             fiche_antivortices[annihilated_antivortex][5] = annihilated_vortex
#         end
#     else # if there is more than one creation or annihilation process
#         error("There is more than one creation/annihilation process. Decrease track_every !")
#     end
#
#     return fiche_vortices,fiche_antivortices
# end


function defect_tracker(thetas::Matrix{Float64},BC,fiche_vortices::Vector{Vector{Any}},fiche_antivortices::Vector{Vector{Any}},vortices_old::Vector{Tuple{Int,Int}},antivortices_old::Vector{Tuple{Int,Int}},t::Float64)
    vortices_new,antivortices_new = spot_defects_separated(thetas,BC)
    m_new,m_old = length(vortices_new),length(vortices_old)

    # Special simple cases to deal with upstream
    if m_new == m_old == 0
        # do nothing, otherwise, "reducing over empty collection blablabla"
    elseif m_new > 0 && m_old == 0
        ll = length(fiche_vortices)
        for i in 1:m_new
            push!(fiche_vortices,[ll+i,vortices_new[i],t,vortices_new[i],nothing,nothing,nothing])
            push!(fiche_antivortices,[ll+i,antivortices_new[i],t,antivortices_new[i],nothing,nothing,nothing])
        end
    elseif m_new == 0 && m_old > 0
        # for i in eachindex(fiche_vortices)
        #     # ici je laisse par flemme l'ID de l'autre defaut = nothing
        #     if fiche_vortices[i][5] == nothing fiche_vortices[i][5] = t end
        #     if fiche_antivortices[i][5] == nothing fiche_antivortices[i][5] = t end
        # end
        ID_annihilated_vortices = []
        ID_annihilated_antivortices = []
        for i in eachindex(fiche_vortices)
            if fiche_vortices[i][5] == nothing
                # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                push!(ID_annihilated_vortices,i)
            end
        end
        for i in eachindex(fiche_antivortices)
            if fiche_antivortices[i][5] == nothing
                push!(ID_annihilated_antivortices,i)
            end
        end
        for i in eachindex(ID_annihilated_vortices)
            # Save the annihilation time
            fiche_vortices[ID_annihilated_vortices[i]][5] = t
            fiche_antivortices[ID_annihilated_antivortices[i]][5] = t # fiche_vortices[i][4] is the time of annihilation

            # Pair up each annihilated vortex with the closest antivortex
            old_loc_vortex = fiche_vortices[ID_annihilated_vortices[i]][2]
            old_loc_antivortex = fiche_antivortices[ID_annihilated_antivortices[i]][2]
            distance = Inf ; ID_antivortex = -1 # dummy
            for j in eachindex(ID_annihilated_antivortices)
                tmp = dist(old_loc_vortex,old_loc_antivortex,L)
                if tmp < distance
                    distance = tmp
                    ID_antivortex = findfirst(x->x[2] == old_loc_antivortex,fiche_antivortices)
                end
            end
            ID_vortex = findfirst(x->x[2] == old_loc_vortex,fiche_vortices)
            fiche_vortices[ID_vortex][7] = ID_antivortex
            fiche_antivortices[ID_antivortex][7] = ID_vortex

            # Estimation of the collision location (midpoint between old vortex and old antivortex)
            dx = abs(old_loc_antivortex[1] - old_loc_vortex[1]) ; dx = min(dx,L-dx)
            dy = abs(old_loc_antivortex[2] - old_loc_vortex[2]) ; dy = min(dy,L-dy)
            estimate_loc_collision = rem.(old_loc_vortex .+ 0.5.*(dx,dy),L)

            fiche_vortices[ID_vortex][6] = estimate_loc_collision
            fiche_antivortices[ID_antivortex][6] = estimate_loc_collision
        end

        for i in eachindex(ID_annihilated_vortices)
            fiche_vortices[ID_annihilated_vortices[i]][5] = fiche_antivortices[ID_annihilated_antivortices[i]][5] = t # fiche_vortices[i][4] is the time of annihilation
            fiche_vortices[ID_annihilated_vortices[i]][7] = ID_annihilated_antivortices[i] # fiche_vortices[i][5] is the ID of the colliding antivortex
            fiche_antivortices[ID_annihilated_antivortices[i]][7] = ID_annihilated_vortices[i]
        end
         # end ???
    else
        distance_matrix_vortices = get_distance_matrix(vortices_new,vortices_old,size(thetas,1)) # m_new lignes, m_old colonnes
        proposal_vortices = hungarian(distance_matrix_vortices)[1] # length = m_new
        assignment_vortices = copy(proposal_vortices)
        distance_matrix_antivortices = get_distance_matrix(antivortices_new,antivortices_old,size(thetas,1))
        proposal_antivortices = hungarian(distance_matrix_antivortices)[1] # length = m_new
        assignment_antivortices = copy(proposal_antivortices)
        #= NB : the proposal is the matching between indices of the vectors
        vortices_new,vortices_old while the assignment matches actual IDs. =#
        for i in eachindex(assignment_vortices)
            for j in eachindex(fiche_vortices)
                if assignment_vortices[i] ≠ 0 && fiche_vortices[j][5] == nothing && fiche_vortices[j][2] == vortices_old[proposal_vortices[i]]
                    # the first condition is a protection against the creation case, where the Hungarian algo matches the newly created vortex to 0
                    # the second condition ensures that one only considers currently living vortices and not vortices now annihilated
                    assignment_vortices[i] = j
                    break
                end
            end
        end
        for i in eachindex(assignment_antivortices)
            for j in eachindex(fiche_antivortices)
                if assignment_antivortices[i] ≠ 0 && fiche_antivortices[j][5] == nothing && fiche_antivortices[j][2] == antivortices_old[proposal_antivortices[i]]
                    assignment_antivortices[i] = j
                    break
                end
            end
        end


    # Case 1 : no creation, no annihilation : simply update the data structure
        if m_new == m_old
            for i in 1:m_new # scan assignment
                fiche_vortices[assignment_vortices[i]][2] = vortices_new[i]
                fiche_antivortices[assignment_antivortices[i]][2] = antivortices_new[i]
            end

    # Case 2 : creation !
        elseif m_new > m_old
            #= locate the newly created defects :
            the assignment vector contains a 0, corresponding to the newly created vortex =#
            ind_created_vortex = findall(iszero,assignment_vortices)
            loc_created_vortex = vortices_new[ind_created_vortex]
            ind_created_antivortex = findall(iszero,assignment_antivortices)
            loc_created_antivortex = antivortices_new[ind_created_antivortex]

            # push them to the data structure
            ll = length(fiche_vortices)
                for i in 1:length(loc_created_vortex)
                    push!(fiche_vortices,[ll + i,loc_created_vortex[i],t,loc_created_vortex[i],nothing,nothing,nothing])
                    push!(fiche_antivortices,[ll + i,loc_created_antivortex[i],t,loc_created_antivortex[i],nothing,nothing,nothing])
                end

            # erase them from the assignment vector and from vortices_new (unused after)
            deleteat!(assignment_vortices,ind_created_vortex)
            deleteat!(vortices_new,ind_created_vortex)
            deleteat!(assignment_antivortices,ind_created_antivortex)
            deleteat!(antivortices_new,ind_created_antivortex)

            # update via the remaining assignment vector
            for i in eachindex(assignment_vortices) # scan assignment modified
                fiche_vortices[assignment_vortices[i]][2] = vortices_new[i]
                fiche_antivortices[assignment_antivortices[i]][2] = antivortices_new[i]
            end

    # Case 3 : annihilation !
        elseif m_new < m_old
             # the annihilated vortex is absent from the assignment vector, so let's proceed
            for i in 1:m_new # scan assignment
                fiche_vortices[assignment_vortices[i]][2] = vortices_new[i]
                fiche_antivortices[assignment_antivortices[i]][2] = antivortices_new[i]
            end

            # now deal with the annihilation stuff
            ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
            for i in eachindex(fiche_vortices)
                if i ∉ assignment_vortices && fiche_vortices[i][5] == nothing
                    # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                    push!(ID_annihilated_vortices,i)
                end
            end
            for i in eachindex(fiche_antivortices)
                if i ∉ assignment_antivortices && fiche_antivortices[i][5] == nothing
                    push!(ID_annihilated_antivortices,i)
                end
            end
            # println()
            # println("m_new = $m_new  m_old = $m_old")
            # println(ID_annihilated_vortices)
            # println(ID_annihilated_antivortices)

            # Now find ID defect of opposite charge closest to each defect
            for i in eachindex(ID_annihilated_vortices)
                # Save the annihilation time
                fiche_vortices[ID_annihilated_vortices[i]][5] = t
                fiche_antivortices[ID_annihilated_antivortices[i]][5] = t # fiche_vortices[i][4] is the time of annihilation

                # Pair up each annihilated vortex with the closest antivortex
                old_loc_vortex = fiche_vortices[ID_annihilated_vortices[i]][2]
                old_loc_antivortex = fiche_antivortices[ID_annihilated_antivortices[i]][2]
                distance = Inf ; ID_antivortex = -1 # dummy
                for j in eachindex(ID_annihilated_antivortices)
                    tmp = dist(old_loc_vortex,old_loc_antivortex,L)
                    if tmp < distance
                        distance = tmp
                        ID_antivortex = findfirst(x->x[2] == old_loc_antivortex,fiche_antivortices)
                    end
                end
                ID_vortex = findfirst(x->x[2] == old_loc_vortex,fiche_vortices)
                fiche_vortices[ID_vortex][7] = ID_antivortex
                fiche_antivortices[ID_antivortex][7] = ID_vortex

                # Estimation of the collision location (midpoint between old vortex and old antivortex)
                dx = abs(old_loc_antivortex[1] - old_loc_vortex[1]) ; dx = min(dx,L-dx)
                dy = abs(old_loc_antivortex[2] - old_loc_vortex[2]) ; dy = min(dy,L-dy)
                estimate_loc_collision = rem.(old_loc_vortex .+ 0.5.*(dx,dy),L)

                fiche_vortices[ID_vortex][6] = estimate_loc_collision
                fiche_antivortices[ID_antivortex][6] = estimate_loc_collision
            end

            for i in eachindex(ID_annihilated_vortices)
                fiche_vortices[ID_annihilated_vortices[i]][5] = fiche_antivortices[ID_annihilated_antivortices[i]][5] = t # fiche_vortices[i][4] is the time of annihilation
                fiche_vortices[ID_annihilated_vortices[i]][7] = ID_annihilated_antivortices[i] # fiche_vortices[i][5] is the ID of the colliding antivortex
                fiche_antivortices[ID_annihilated_antivortices[i]][7] = ID_annihilated_vortices[i]
            end
        end
    end

    return fiche_vortices,fiche_antivortices
end

function get_distance_matrix(vortices_new,vortices_old,L)
    m_new,m_old = length(vortices_new),length(vortices_old)
    distance_matrix = zeros(m_new,m_old)
    for j in 1:m_old
        for i in 1:m_new
            distance_matrix[i,j] = dist(vortices_new[i],vortices_old[j],L)
        end
    end
    return distance_matrix
end


function lifetimes(fiche,tmax)
    lifetimes = zeros(length(fiche))
    for i in eachindex(fiche)
        if fiche[i][4] ≠ nothing
            lifetimes[i] = fiche[i][4] - fiche[i][3]
        else
            lifetimes[i] = tmax - fiche[i][3]
        end
    end
    return lifetimes
end

## Time evolution methods
function update(thetas::Matrix{Float64},omegas::Matrix{Float64},L::Int,T::Number,dt::Float64)::Array{Float64}
    thetas_old = copy(thetas)
    for j in 1:L
        for i in 1:L
            θ = thetas_old[i,j]
            angle_neighbours = get_angle_neighbours(thetas_old,i,j,L) # O(1)
            sin_angle_neighbours = sin.(angle_neighbours .- θ)
            thetas[i,j] =  θ + omegas[i,j]*dt + dt*sum(sin_angle_neighbours) + sqrt(2T*dt)*randn() # O(1)
        end
    end
    return thetas
end

function update_FBC(thetas::Matrix{Float64},omegas::Matrix{Float64},L::Int,T::Number,dt::Float64)::Array{Float64}
    thetas_old = copy(thetas)
    for i in 1:L # j =1
        θ = thetas_old[i,1]
        thetas[i,1] =  θ + omegas[i,1]*dt + dt*sum(sin.(get_angle_neighbours_borders(thetas_old,i,1,L) .- θ)) + sqrt(2T*dt)*randn() # O(1)
    end

    for i in 1:L # j = L
        θ = thetas_old[i,L]
        thetas[i,L] =  θ + omegas[i,L]*dt + dt*sum(sin.(get_angle_neighbours_borders(thetas_old,i,L,L) .- θ)) + sqrt(2T*dt)*randn() # O(1)
    end

    for j in 1:L # i = 1
        θ = thetas_old[1,j]
        thetas[1,j] =  θ + omegas[1,j]*dt + dt*sum(sin.(get_angle_neighbours_borders(thetas_old,1,j,L) .- θ)) + sqrt(2T*dt)*randn() # O(1)
    end

    for j in 1:L # i = L
        θ = thetas_old[L,j]
        thetas[L,j] =  θ + omegas[L,j]*dt + dt*sum(sin.(get_angle_neighbours_borders(thetas_old,L,j,L) .- θ)) + sqrt(2T*dt)*randn() # O(1)
    end

    for j in 2:L-1, i in 2:L-1
        θ = thetas_old[i,j]
        thetas[i,j] =  θ + omegas[i,j]*dt + dt*sum(sin.(get_angle_neighbours(thetas_old,i,j,L) .- θ)) + sqrt(2T*dt)*randn() # O(1)
    end
    return thetas
end

# function update_underdamped(thetas::Matrix{Float64},thetas_old::Matrix{Float64},omegas::Matrix{Float64},L::Int,T::Number,dt::Float64,m)::Array{Float64}
#     thetas_new = zeros(L,L)
#     for j in 1:L
#         for i in 1:L
#             angle_neighbours = get_angle_neighbours(thetas,i,j,L) # O(1)
#             sin_angle_neighbours = sin.(angle_neighbours .- thetas[i,j])
#             thetas_new[i,j] =  (omegas[i,j] + sum(sin_angle_neighbours) + sqrt(2T/dt)*randn() + thetas[i,j]*(1/dt+2m/dt^2) - m/dt^2*thetas_old[i,j]) / (m/dt^2 + 1/dt)
#         end
#     end
#     return thetas_new
# end

function update_underdamped_FBC(thetas::Matrix{Float64},thetas_old::Matrix{Float64},omegas::Matrix{Float64},L::Int,T::Number,dt::Float64,m)::Array{Float64}
    thetas_new = zeros(L,L)
    for i in 1:L # j = 1
        angle_neighbours = get_angle_neighbours_borders(thetas,i,1,L) # O(1)
        sin_angle_neighbours = sin.(angle_neighbours .- thetas[i,1])
        thetas_new[i,1] =  (omegas[i,1] + sum(sin_angle_neighbours) + sqrt(2T/dt)*randn() + thetas[i,1]*(1/dt+2m/dt^2) - m/dt^2*thetas_old[i,1]) / (m/dt^2 + 1/dt)
    end
    for i in 1:L # j = L
        angle_neighbours = get_angle_neighbours_borders(thetas,i,L,L) # O(1)
        sin_angle_neighbours = sin.(angle_neighbours .- thetas[i,L])
        thetas_new[i,L] =  (omegas[i,L] + sum(sin_angle_neighbours) + sqrt(2T/dt)*randn() + thetas[i,L]*(1/dt+2m/dt^2) - m/dt^2*thetas_old[i,L]) / (m/dt^2 + 1/dt)
    end
    for j in 1:L # i = 1
        angle_neighbours = get_angle_neighbours_borders(thetas,1,j,L) # O(1)
        sin_angle_neighbours = sin.(angle_neighbours .- thetas[1,j])
        thetas_new[1,j] =  (omegas[1,j] + sum(sin_angle_neighbours) + sqrt(2T/dt)*randn() + thetas[1,j]*(1/dt+2m/dt^2) - m/dt^2*thetas_old[1,j]) / (m/dt^2 + 1/dt)
    end
    for j in 1:L # i = L
        angle_neighbours = get_angle_neighbours_borders(thetas,L,j,L) # O(1)
        sin_angle_neighbours = sin.(angle_neighbours .- thetas[L,j])
        thetas_new[L,j] =  (omegas[L,j] + sum(sin_angle_neighbours) + sqrt(2T/dt)*randn() + thetas[L,j]*(1/dt+2m/dt^2) - m/dt^2*thetas_old[L,j]) / (m/dt^2 + 1/dt)
    end

    for j in 2:L-1 , i in 2:L-1
        angle_neighbours = get_angle_neighbours(thetas,i,j,L) # O(1)
        sin_angle_neighbours = sin.(angle_neighbours .- thetas[i,j])
        thetas_new[i,j] =  (omegas[i,j] + sum(sin_angle_neighbours) + sqrt(2T/dt)*randn() + thetas[i,j]*(1/dt+2m/dt^2) - m/dt^2*thetas_old[i,j]) / (m/dt^2 + 1/dt)
    end
    return thetas_new
end


function update_underdamped_leapfrog(Q::Matrix{Float64},P::Matrix{Float64},omegas::Matrix{Float64},L::Int,T::Number,dt::Float64,m)
    #= Demarche :
        * the equation is q' = p && p' + p/m = f + noise , with p = dtheta/dt and q = theta
        * first, one integrates the Langevin equation with \gamma = 0, ie p'(t) = f(t)
        * then, one introduces the damping and the noise p(t+dt) = p(t) (1-dt/sqrt(m)) + sqrt(2T*dt/sqrt(m))*randn()
    =#
    Q_new = zeros(L,L)
    for j in 1:L , i in 1:L
        angle_neighbours_old = get_angle_neighbours(Q,i,j,L) # O(1)
        sin_angle_neighbours_old = sin.(angle_neighbours_old .- Q[i,j])
        Q_new[i,j] = Q[i,j] + dt * P[i,j] + (dt^2) / 2 * (omegas[i,j] + sum(sin_angle_neighbours_old))
    end
    P_new = zeros(L,L)
    for j in 1:L , i in 1:L
        P_new[i,j] = P[i,j] + dt/2 * (2omegas[i,j] + sum(sin.(get_angle_neighbours(Q,i,j,L) .- Q[i,j])) + sum(sin.(get_angle_neighbours(Q_new,i,j,L) .- Q_new[i,j])))
    end
    P_new = P_new*(1-dt/sqrt(m)) + sqrt(2T*dt/sqrt(m))*randn(L,L)
    return Q_new,P_new
end

function update_underdamped_leapfrog_FBC(Q::Matrix{Float64},P::Matrix{Float64},omegas::Matrix{Float64},L::Int,T::Number,dt::Float64,m)
    Q_new = zeros(L,L)
    for i in 1:L Q_new[i,1] = Q[i,1] + dt * P[i,1] +  dt^2 / 2 * (omegas[i,1] + sum( sin.(get_angle_neighbours_borders(Q,i,1,L) .- Q[i,1]) )) end
    for i in 1:L Q_new[i,L] = Q[i,L] + dt * P[i,L] +  dt^2 / 2 * (omegas[i,L] + sum( sin.(get_angle_neighbours_borders(Q,i,L,L) .- Q[i,L]) )) end
    for j in 1:L Q_new[1,j] = Q[1,j] + dt * P[1,j] +  dt^2 / 2 * (omegas[1,j] + sum( sin.(get_angle_neighbours_borders(Q,1,j,L) .- Q[1,j]) )) end
    for j in 1:L Q_new[L,j] = Q[L,j] + dt * P[L,j] +  dt^2 / 2 * (omegas[L,j] + sum( sin.(get_angle_neighbours_borders(Q,L,j,L) .- Q[L,j]) )) end
    for j in 2:L-1 , i in 2:L-1 Q_new[i,j] = Q[i,j] + dt * P[i,j] +  dt^2 / 2 * (omegas[i,j] + sum( sin.(get_angle_neighbours(Q,i,j,L) .- Q[i,j]) )) end

    P_new = zeros(L,L)
    for i in 1:L  P_new[i,1] = P[i,1] + dt/2 * (2omegas[i,1] + sum(sin.(get_angle_neighbours_borders(Q,i,1,L) .- Q[i,1])) + sum(sin.(get_angle_neighbours_borders(Q_new,i,1,L) .- Q_new[i,1])) ) end
    for i in 1:L  P_new[i,L] = P[i,L] + dt/2 * (2omegas[i,L] + sum(sin.(get_angle_neighbours_borders(Q,i,L,L) .- Q[i,L])) + sum(sin.(get_angle_neighbours_borders(Q_new,i,L,L) .- Q_new[i,L])) ) end
    for j in 1:L  P_new[1,j] = P[1,j] + dt/2 * (2omegas[1,j] + sum(sin.(get_angle_neighbours_borders(Q,1,j,L) .- Q[1,j])) + sum(sin.(get_angle_neighbours_borders(Q_new,1,j,L) .- Q_new[1,j])) ) end
    for j in 1:L  P_new[L,j] = P[L,j] + dt/2 * (2omegas[L,j] + sum(sin.(get_angle_neighbours_borders(Q,L,j,L) .- Q[L,j])) + sum(sin.(get_angle_neighbours_borders(Q_new,L,j,L) .- Q_new[L,j])) ) end
    for j in 2:L-1 , i in 2:L-1 P_new[i,j] = P[i,j] + dt/2 * (2omegas[i,j] + sum(sin.(get_angle_neighbours(Q,i,j,L) .- Q[i,j])) + sum(sin.(get_angle_neighbours(Q_new,i,j,L) .- Q_new[i,j])) ) end

    P_new = P_new*(1-dt/sqrt(m)) + sqrt(2T*dt/sqrt(m))*randn(L,L)
    return Q_new,P_new
end


function evolve_FBC(Q::Int,thetas::Matrix{Float64},omegas::Matrix{Float64},T::Number,dt::Float64,tmax::Number,every)
    L = size(thetas)[1]
    nsteps_max = floor(Int,tmax/dt)
    history_locations = Vector{Union{Tuple{Int16,Int16},Missing}}(missing,1+floor(Int,nsteps_max/every))
    # Initial location
    last_location,~ = spot_single_default_global(thetas)
    history_locations[1] = last_location

    # Time Evolution
    for i in 1:nsteps_max
        thetas = update_FBC(thetas,omegas,L,T,dt) # O(L²)
        if i%every == 0
            try
                lastknown = (Int16(-1),Int16(-1))
                last_location,alone = spot_single_default_local(thetas,last_location,lastknown,Q)
                if alone # there are no other vortices in the surroundings (the '2margin' square)
                    saved_location = last_location
                else # there are other vortices in the surroundings that could pollute the data
                    saved_location = missing
                end
                history_locations[1+div(i,every)] = saved_location  # O(1)
            catch e  # we lost our vortex (out of lattice, annihilation)
                println(e)
                printstyled("Warning : Vortex lost, simulation stopped at t = $(dt*i*every). \n"; color = :yellow)
                break # we leave all the remaining values as they are : missing
            end
        end

    end

    return history_locations
end

function spot_single_default_global(thetas::Matrix{Float64})::Tuple{Tuple{Int16,Int16},Int}
    L = size(thetas)[1]
    list_defaults = Tuple{Int16,Int16,Int}[]
    for i in 2:L-1
        for j in 2:L-1
            q = get_vorticity(thetas,i,j,L)
            if q ≠ 0
                push!(list_defaults,(i,j,q))
            end
        end
    end
    return list_defaults[1][1:2],list_defaults[1][3]
end


function spot_single_default_local(thetas::Matrix{Float64},last_loc::Tuple{Int16,Int16},known_loc::Tuple{Int16,Int16},Q::Int,margin::Int=6)::Tuple{Tuple{Int16,Int16},Bool}
    # V    is the location (x,y,q) of the default
    #= The search for the new location of the default will scan a square of
    dimension 2margin x 2margin (within the limits of the LxL lattice) ,
    with the last location as center. The whole method is O(margin²). =#
    L = size(thetas)[1]
    positions = []


    # Dimensions of the square
    j_low,j_high = max(1,last_loc[2]-margin) , min(L,last_loc[2]+margin)
    i_low,i_high = max(1,last_loc[1]-margin) , min(L,last_loc[1]+margin)
    for j in j_low:j_high
        for i in i_low:i_high
            q = get_vorticity(thetas,i,j,L)
            if q ≠ 0
                push!(positions,(i,j,q))
            end
        end
    end


    ℓ = length(positions)
    #= Summary of what follows :
    ℓ = 0 : (i)   Controled Annihilation (encounter of the vortex with its antivortex)
            (ii)  Unexpected Annihilation (encounter of the vortex with another antivortex
            (iii) The (single) vortex left the lattice.
            In all three case, throw an explicit error and treat it later on.

    ℓ = 1 : All good ! One is pretty sure that the only detected vortex is the same
    than the one of the previous timestep. The only case where we could be mistaken
    is if a pair of vortex enters the square and that our previous vortex annihilates
    with the newcomer antivortex. In this case, the only default remaining would be
    the newcomer vortex, and we would be clueless about it. The only signature would be a
    possible jump in space in the r(t) graph.

    ℓ = 2 : If the other default is known/authorized, meaning that it is the former
    antivortex of the default we currently work on, that's not important and we keep
    going as if the vortex was alone. If the other defaut is NOT known/authorized,
    it will perturb the displacment of our vortex : don't save the current position
    by signalling it (alone = false)

    ℓ ≥ 3 : We are sure to be in the bad case #2.  =#
    if ℓ == 1
        alone = true # no other default in the searched region
        most_probable = positions[1][1:2]
        #= Actually,  =#
    elseif ℓ > 1
        alone = false # by default. We deal we the special case ℓ = 2 just below.

        if ℓ == 2 && ((known_loc[1],known_loc[2],-Q) in positions) ; alone = true ; end
        #= If ℓ == 2, one has 2 possibilities :
            - either the extra default is "known" (its the antidefault of
            the one we currently consider), in which case turn alone = true.
            - or it's not, in which case we leave alone = false.
            In any case, the following routine to determine the location of
            the default currenlty tracked is still completely valid.
        Note that when using this function for the study of a single vortex,
        one needs to provide an impossible known_loc, such as (-1,-1). =#

        distances_to_last = [dist(positions[i][1:2],last_loc,L) for i in 1:ℓ]
        # Kill candidates with opposite polarity
        for i in 1:ℓ
            element = positions[i]
            if element[3] ≠ Q distances_to_last[i] = Inf end
        end
        most_probable = positions[sortperm(distances_to_last)][1]
        #= Returns the position with the smallest distance to the last location,
        hence we choose that one as the most probable candidate for the vortex we
        are considering. =#

    else # dealing with the case ℓ = 0

        close_from_boundary = last_loc[1] < 4 || last_loc[1] > L - 4 || last_loc[2] < 4 || last_loc[2] > L - 4
        if close_from_boundary                      error("The vortex left the lattice.")
        elseif dist(last_loc,known_loc,L) ≤ sqrt(5) error("Controlled annihilation.") # sqrt(5) is arbitrary
        else                                        error("Unexpected annihilation.") end
    end
    # println(1)
    # println(typeof(Int16.(most_probable[1:2])))
    return most_probable[1:2],alone
end

## Meta Methods
function scan(configfile_name="configuration.jl")
    include(configfile_name)
    include("IDrealisation.jl")
    len = length(TV)

    # Declarations of the output containers
    Δ   = SharedArray{Float64,3}(len,2,nsave) # global magnetisation on the LXL system (genuine and projected on x-axis)
    δ   = SharedArray{Float64,4}(len,2,Int(L/l)^2,nsave) # local magnetisation on a lxl domain (genuine and projected on x-axis)
    C   = SharedArray{Float16,3}(len,L,nsave)
    Ctw = SharedArray{Float64,3}(len,length(tws),nsave)
    number_all_vortex  = SharedArray{Int,2}(len,nsave)
    number_free_vortex = SharedArray{Int,2}(len,nsave)
    # thetas_t           = SharedArray{Float16,3}(len,L^2,nsave) # vec(thetas) for each time saved

    LS = fill(L,len) ; lS = fill(l,len) ; initS = fill(init,len) ; tsaveS = fill(tsave,len) ; tmaxS = fill(tmax,len) ; twsS = fill(tws,len)
    ΔS = fill(Δ,len) ; δs = fill(δ,len) ; CS = fill(C,len) ; CtwS = fill(Ctw,len) ; number_all_vortexS = fill(number_all_vortex,len) ; number_free_vortexS = fill(number_free_vortex,len)
    pmap(simulation,LS,lS,TV,initS,tmaxS,tsaveS,twsS,Array(1:len),ΔS,δs,CS,CtwS,number_all_vortexS,number_free_vortexS)

    nT,nVar = length(Ts),length(Vars)
    # Convert the containers to the right format
    Δ   = convert(Array,reshape(Δ,nT,nVar,2,nsave))
    δ   = convert(Array,reshape(δ,nT,nVar,2,Int(L/l)^2,nsave))
    C   = convert(Array,reshape(C,nT,nVar,L,nsave))
    Ctw = convert(Array,reshape(Ctw,nT,nVar,length(tws),nsave))
    number_all_vortex  = convert(Array,reshape(number_all_vortex,nT,nVar,nsave))
    number_free_vortex = convert(Array,reshape(number_free_vortex,nT,nVar,nsave))
    # thetas_t           = convert(Array,reshape(thetas_t,nT,nVar,L^2,nsave))

    ## Save Outputs
    # println("On y arrive.")
    filename = "data/"*fname*"_r$r.jld"
    JLD.save(filename,"Δ",Δ,"δ",δ,"C",C,"Ctw",Ctw,"number_all_vortex",number_all_vortex,"number_free_vortex",number_free_vortex,"L",L,"l",l,"Vars",Vars,"Ts",Ts,"init",init,"tmax",tmax,"tsave",tsave,"tws",tws,"name",name,"comments",comments)
    println("Data saved in ",filename)
end

function simulation(L,l,TV,init,tmax,tsave,tws,i,Δ,δ,C,Ctw,number_all_vortex,number_free_vortex)
    T,Var = TV
    println("Simulation at T=$T & σ²=$Var launched at $(Dates.hour(now())):$(Dates.minute(now()))")
    lattice,dt = create_Lattice(L,T,Var,init)
    Δ[i,:,:],δ[i,:,:,:],C[i,:,:],Ctw[i,:,:],number_all_vortex[i,:],number_free_vortex[i,:] = evolve(lattice,l,T,dt,tmax,tsave,tws)
end

## Visualisation methods
function smooth(X) ## for smoother plots
    smoothed = copy(X)
    coeff = [1,2,1]
    coeff = coeff./sum(coeff)
    # @assert isodd(length(coeff))
    s = Int(floor(length(coeff)/2))
    for i in 1+s:length(smoothed)-s
        smoothed[i] = X[i-s:i+s]'*coeff
    end
    return smoothed
end

function smooth(X;over=3) ## for smoother plots
    smoothed = copy(X)
    coeff = [2^(i-1) for i in 1:over]
    coeff = coeff./sum(coeff)
    s = length(coeff)
    for i in 1+s:length(smoothed)
        smoothed[i] = X[i-s+1:i]'*coeff
    end
    return smoothed
end


function correct_location(x::Tuple{Int,Int},L::Int)::Tuple{Int,Int}
    #= This method deals with the difference in ordering between the lattice of
        angles and the way plots are layed out. Thetas spans from (1,1) in top left
        corner to (L,L) at bottom right. Plots layout goes from (1,1) in bottom left to (L,L) at
        top right. Thus, one needs to invert rows&columns and to reverse the order of the (new)
        columns (going down instead of going up).
        One can find these explanations with schemes at page 21 of the handwritten book. =#
    return (x[2],L-x[1])
end

# function twiny(sp::Plots.Subplot)
#     sp[:top_margin] = max(sp[:top_margin], 40Plots.px)
#     plot!(sp.plt, inset = (sp[:subplot_index], bbox(0,0,1,1)))
#     twinsp = sp.plt.subplots[end]
#     twinsp[:xaxis][:mirror] = true
#     twinsp[:background_color_inside] = RGBA{Float64}(0,0,0,0)
#     Plots.link_axes!(sp[:yaxis], twinsp[:yaxis])
#     return twinsp
# end
# twiny(plt::Plots.Plot = current()) = twiny(plt[1])
# #

function display_thetas(thetas;defects=false,title="",colorbar=true)
    cols=cgrad([:black,:blue,:green,:orange,:red,:black])
    if defects
        p = heatmap(mod.(thetas',2π),c=cols,title=title,clims=(0,2π),size=(512,512),colorbar=colorbar)
        highlight_defects(p,thetas)
    else
        p = heatmap(mod.(thetas',2π),c=cols,title=title,clims=(0,2π),size=(512,512),colorbar=colorbar)
    end
    return p
end

function highlight_defects(p,thetas,symbP=:circle,symbM=:utriangle)
    L = size(thetas,1)
    for defect in spot_defects(thetas,"periodic")
        if defect[3] > 0 symb = symbP else symb = symbM end
        scatter!((defect[1:2]), m = (8, 12.0, symb,:transparent, stroke(2, :white)))
    end
    xlims!((0,L))
    ylims!((0,L))
    return p
end

## Miscellaneous
function index_element_repeated(x)
    indices = []
    L = length(x)
    for i in 1:L
        count = 0
        for j in 1:L
            if x[j] == x[i] count += 1 end
        end
        if count > 1 push!(indices,i) end
    end
    return indices
end

# function remove_negative(series)
#     result = Vector{Union{Float64,Missing}}(undef,length(series)) # initializes by default to "missing"
#     for i in eachindex(series)
#         if series[i] > 0
#             result[i] = series[i]
#         end
#     end
#     return result
# end

function remove_negative(input)
    array = Float64.(copy(input))
    for i in 1:length(array)
        if array[i] ≤ 0 array[i] = NaN end
    end
    return array
end

function number_decimals(x::Number)::Int
    ndigits = 0
    while !isinteger(x) && ndigits < 4 # gardefou pour les nombres à écriture décimale infinie
        x = x * 10
        ndigits += 1
    end

    return ndigits
end

function grad2(thetas) # returns ||∇θ(x,t)||²
    L     = size(thetas)[1]
    grad2 = zeros(L,L)
    # thetas= mod.(thetas,2π)
    for j in 1:L
        for i in 1:L
            neighbours = get_angle_nearest_neighbours(thetas,i,j,L)
            dθdx       = (neighbours[1] - neighbours[2]) / 2
            dθdy       = (neighbours[3] - neighbours[4]) / 2
            grad2[i,j] = max(abs(dθdx), abs(dθdy))
            # grad2[i,j] = 1/min(abs(dθdx),abs(dθdy))
            # grad2[i,j] = dθdx^2 + dθdy^2
        end
    end
    return grad2
end

function walls(thetas;tol=π/4)
    L   = size(thetas)[1]
    walls = zeros(L,L)
    # tol_aligned = π/4
    # tol_perp    =
    for j in 1:L
        for i in 1:L
            bool_vertical   = abs(thetas[mod1(i+1,L),j] - thetas[mod1(i-1,L),j]) < tol && abs(π - abs(thetas[i,mod1(j+2,L)] - thetas[i,mod1(j+2,L)])) < tol
            bool_horizontal = abs(π - abs(thetas[mod1(i+2,L),j] - thetas[mod1(i-2,L),j])) < tol && abs(thetas[i,mod1(j+1,L)] - thetas[i,mod1(j+1,L)]) < tol
            # bool_diagonal   = abs(π - abs(thetas[mod1(i+1,L),j] - thetas[mod1(i-1,L),j])) < tol && abs(thetas[i,mod1(j+1,L)] - thetas[i,mod1(j+1,L)]) < tol
            walls[i,j] = bool_horizontal || bool_vertical
        end
    end
    return walls
end

function energy(thetas)
    L     = size(thetas)[1]
    energy_local = zeros(L,L)
    for j in 1:L
        for i in 1:L
            neighbours = get_angle_neighbours(thetas,i,j,L)
            θ = thetas[i,j]
            energy_local[i,j] = -sum(cos.(neighbours .- θ))
        end
    end
    energy_global = sum(energy_local)
    return energy_local,energy_global
end

function psi(thetas)
    L     = size(thetas)[1]
    m = [mean(cos.(thetas)),mean(sin.(thetas))]
    u = m / (sqrt(m[1]^2 + m[2]^2))

    psi = zeros(size(thetas))
    for i in eachindex(thetas)
        psi[i] = u[1]*cos(thetas[i]) + u[2]*sin(thetas[i])
    end
    return psi
end

function grad_sq(matrix)
    L     = size(thetas)[1]
    grad_sq = zeros(L-2,L-2)
    for i in 2:L-1
        for j in 2:L-1
            dx = (matrix[i+1] - matrix[i-1])/2
            dy = (matrix[j+1] - matrix[j-1])/2
            grad_sq[i-1,j-1] = sqrt(dx^2 + dy^2)
        end
    end
    return grad_sq
end

function compute_gr(thetas,dr=1)
    L = size(thetas,1)
    vortices,antivortices = spot_defects_separated(thetas,BC)
    m = length(vortices)

    rr = collect(0:dr:Int(L/2))
    aire_couronne = π*(2 .* rr[1:end-1] * dr .+ dr^2)

    # Compute distances between opposite charge defects
    distances_opposite = zeros(m,m)
    for i in 1:m , j in 1:m
        distances_opposite[i,j] = dist(vortices[i],antivortices[j],L)
    end

    # Compute distances between same charge defects
    distances_same = NaN*zeros(m,m)
    for i in 1:m , j in i+1:m
        distances_same[i,j] = dist(vortices[i],vortices[j],L)
        distances_same[j,i] = dist(antivortices[i],antivortices[j],L)
    end

    gr_opposite = fit(Histogram,vec(distances_opposite),rr).weights ./ aire_couronne / m^2 * L^2
    gr_same = fit(Histogram,vec(distances_same),rr).weights ./ aire_couronne / m^2 * L^2

    return gr_opposite,gr_same
end

function corr_length(rs,C::Vector{Float64},seuil=exp(-1))::Float64 # from a time series, returns the correlation length ()
    i_after = findfirst(x->x<seuil,C)
    if i_after ≠ nothing && i_after > 1
        # Linear interpolation
        i_befor = i_after - 1
        r_after = rs[i_after]
        r_befor = rs[i_befor]
        c_after = C[i_after]
        c_befor = C[i_befor]
        ξ = (seuil*(r_after-r_befor) -  (c_befor*r_after - r_befor*c_after))/(c_after-c_befor)
    else
    ξ = NaN
    end
    return ξ
end

function exp_moving_average(x,window)
    y = zeros(length(x))
    for t in 1:length(y)
        y[t] = sum([x[t-i]*((window-1)/window)^(i) for i in 0:t-1])/window
    end
    return y
end

function density_creation(fiche,L,cg_spacing)
    locations = location_creation(fiche)
    density_matrix = zeros(Int,Int(L/cg_spacing),Int(L/cg_spacing))
    for element in locations
        i = 1 + floor(Int,element[1]/cg_spacing)
        j = 1 + floor(Int,element[2]/cg_spacing)
        density_matrix[i,j] += 1
    end
    density_vec = vec(density_matrix)
    return density_matrix,density_vec
end

## Methods pour retrouver le 3/2
function get_neighbours(i::Int,j::Int,L::Int,topology)::Vector{Tuple{Int,Int}}
    if topology == "chebychev"
        indices_neighbours = [(mod1(i-1,L),mod1(j,L)) , (mod1(i,L),mod1(j-1,L)) , (mod1(i,L),mod1(j+1,L)) , (mod1(i+1,L),mod1(j,L)) , (mod1(i+1,L),mod1(j+1,L)) , (mod1(i-1,L),mod1(j-1,L)) , (mod1(i-1,L),mod1(j+1,L)) , (mod1(i+1,L),mod1(j-1,L))]
    elseif topology == "hexagonal"
        indices_neighbours = [(mod1(i-1,L),mod1(j,L)) , (mod1(i,L),mod1(j-1,L)) , (mod1(i,L),mod1(j+1,L)) , (mod1(i+1,L),mod1(j,L)) , (mod1(i+1,L),mod1(j+1,L)) , (mod1(i-1,L),mod1(j-1,L))] # hexagonal
    elseif topology == "square"
        indices_neighbours = [(mod1(i-1,L),mod1(j,L)) , (mod1(i,L),mod1(j-1,L)) , (mod1(i,L),mod1(j+1,L)) , (mod1(i+1,L),mod1(j,L))] # square
    else
        error("Please choose the topology among 'chebychev', 'hexagonal' and 'square'.")
    end
    return indices_neighbours
end

function step(grid,L,topology)
    grid_new = copy(grid)
    for j in 1:L , i in 1:L
        if grid[i,j] == 0
            indices_neighbours = get_neighbours(i,j,L,topology)
            ID_neighbours = [grid[indices_neighbours[a]...] for a in eachindex(indices_neighbours)]
            grid_new[i,j] = maximum(ID_neighbours)
        end
    end
    return grid_new
end

function init_grid(L,nb_nucleation)
    grid = zeros(L,L)
    nucleation_sites = [(rand(1:L),rand(1:L)) for i in 1:nb_nucleation]
    for i in eachindex(nucleation_sites)
        grid[nucleation_sites[i]...] = i
    end
    return grid,nucleation_sites
end

function DB(grid,L,k,thicken)
    db = zeros(L,L)
    for j in 1:L , i in 1:L
        indices_neighbours = get_neighbours(i,j,L,topology)
        ID_neighbours = [grid[indices_neighbours[a]...] for a in eachindex(indices_neighbours)]
        if sum(grid[i,j] .- ID_neighbours) > 0
            db[i,j] = 1
        end
    end

    # thicken a bit the domain boundaries
    if thicken
        for j in 1:L , i in 1:L
            indices_neighbours = get_neighbours(i,j,L,topology)
            db_neighbours = [db[indices_neighbours[a]...] for a in eachindex(indices_neighbours)]
            if sum(db_neighbours .== 1) ≥ 1 && db[i,j] == 0
                db[i,j] = .99
            end
        end
    end

    # expand artificially the matrix, ensuring PBC
    if k == 1
        return db
    elseif k > 1
        space = zeros(k*L,k*L)
        for n in 0:k-1 , m in 0:k-1
            space[1+n*L:(n+1)*L,1+m*L:(m+1)*L] = db
        end
        return space
    end
end

# fonction hyperlente
# function create_space(L;k=1,metric,topology="square",nb_nucleation=Int(round(log(L),RoundUp)),thicken=false)
#     grid,~ = init_grid(L,nb_nucleation)
#     while minimum(grid) == 0
#         grid = step(grid,L,metric,topology)
#     end
#     return DB(grid,L,k,thicken)
# end

# même fonction mais plus rapide
function create_space(L;k=1,metric="manhattan",nb_nucleation=Int(round(log(L),RoundUp)),thicken=false)
    grid,nucleation_sites = init_grid(L,nb_nucleation)
    if metric == "manhattan"
        for j in 1:L , i in 1:L
            if grid[i,j] == 0
                distances_to_nucl = [manhattan_dist((i,j),nuc,L) for nuc in nucleation_sites]
                grid[i,j] = argmin(distances_to_nucl)
            end
        end
    elseif metric == "euclidian"
        for j in 1:L , i in 1:L
            if grid[i,j] == 0
                distances_to_nucl = [dist2((i,j),nuc,L) for nuc in nucleation_sites]
                grid[i,j] = argmin(distances_to_nucl)
            end
        end
    end
    return DB(grid,L,k,thicken)
end


function direction2neighbour(dir::Vector{Float64},i::Int,j::Int,L::Int,topology)
    if norm(dir) > 0
        anglee = atan(dir[2],dir[1])
        if 0 ≤ anglee ≤ π/4 nn = (1,0) # 1er huitieme
        elseif π/4 < anglee < 3π/4 nn = (0,1) # 2e et 3e huitiemes
        elseif 3π/4 ≤ anglee ≤ π nn = (-1,0) # 4e huitieme
        elseif -3π/4 ≥ anglee ≥ -π nn = (-1,0) # 5e huitieme
        elseif -π/4 > anglee > -3π/4 nn = (0,-1) # 6e et 7e huitieme
        elseif 0 ≥ anglee ≥ -π/4 nn = (1,0) # 1er huitieme
        else error("J'ai du oublier un cas de figure, fonction a revoir !")
        end
        return mod1.((i,j) .+ nn,L)
    else
        return sample(get_neighbours(i,j,L,topology))
    end
end

function closetoborders(i::Int,j::Int,L::Int,c::Int=5)
    if     i ≤ c     return true
    elseif j ≤ c     return true
    elseif i + c ≥ L return true
    elseif j + c ≥ L return true
    else return false
    end
end

# function update((i,j)::Tuple{Int,Int},vision_radius::Int,space::Matrix{Float64},L::Int; alpha=0.8,metric,topology="square")
#     col = :red
#     if closetoborders(i,j,L,2vision_radius)
#         shiftedspace = circshift(space,-1 .* (i-Int(L/2),j-Int(L/2)))
#         vision = shiftedspace[Int(L/2)-vision_radius:Int(L/2)+vision_radius , Int(L/2)-vision_radius:Int(L/2)+vision_radius]
#
#     else
#         vision = space[i-vision_radius:i+vision_radius , j-vision_radius:j+vision_radius]
#     end
#
#     if sum(vision) == 0 # no visible boundary, perform simple RW
#         i,j = sample(get_neighbours(i,j,L,metric,topology))
#     else # there is a boundary within the vision radius
#         u = rand()
#         if u < alpha # with proba alpha, get attracted towards this boundary
#             vector_attraction_to_boundary = [0,0]
#             R = 2vision_radius + 1
#             for m in 1:R , n in 1:R
#                 vec_tmp = [n-(R+1)/2,m-(R+1)/2]
#                 if !iszero(vec_tmp) vector_attraction_to_boundary += vision[n,m]*vec_tmp/norm(vec_tmp)^2 * (1 .+ randn()*0.01) end
#             end
#             i,j = direction2neighbour(vector_attraction_to_boundary,i,j,L,metric,topology)
#             col = :green
#         else # with proba 1-alpha, simple RW
#             i,j = sample(get_neighbours(i,j,L,metric,topology))
#         end
#     end
#     space[i,j] = 0 # the defect devours the boundary
#     return (i,j),space, col
# end

function update(loc::Vector{Tuple{Int,Int}},vision_radius::Int,space::Matrix{Float64},L::Int; alpha=0.8,topology="square",erasure=0.0)
    N = length(loc)
    col = fill(:red,N)
    for n in 1:N
        i,j = loc[n]
        if closetoborders(i,j,L,2vision_radius)
            shiftedspace = circshift(space,-1 .* (i-Int(L/2),j-Int(L/2)))
            vision = shiftedspace[Int(L/2)-vision_radius:Int(L/2)+vision_radius , Int(L/2)-vision_radius:Int(L/2)+vision_radius]

        else
            vision = space[i-vision_radius:i+vision_radius , j-vision_radius:j+vision_radius]
        end

        if sum(vision .== 1) == 0 # no visible boundary, perform simple RW
            i,j = sample(get_neighbours(i,j,L,topology))
        else # there is a boundary within the vision radius
            u = rand()
            if u < alpha # with proba alpha, get attracted towards this boundary
                vector_attraction_to_boundary = vision2force(vision)
                i,j = direction2neighbour(vector_attraction_to_boundary,i,j,L,topology)
                col[n] = :green
            else # with proba 1-alpha, simple RW
                i,j = sample(get_neighbours(i,j,L,topology))
            end
        end

        space[max(i-1,1):min(i+1,L),max(j-1,1):min(j+1,L)] .= erasure # the defect devours the boundary
        # space[i,j] = erasure # the defect devours the boundary
        loc[n] = i,j
    end
    return loc,space,col
end

function vision2force(vision)
    R = size(vision,1)
    vector_attraction_to_boundary = [0,0]
    for m in 1:R , n in 1:R
        vec_tmp = [n-(R+1)/2,m-(R+1)/2]
        if !iszero(vec_tmp) vector_attraction_to_boundary += vision[n,m]*vec_tmp/norm(vec_tmp)^2 * (1 .+ randn()*0.001) end
    end
    return vector_attraction_to_boundary
end

#
# vision = zeros(5,5) ; vision[1,:] .= 1  ; vision[3,1] = -1
# vision
# vector_attraction_to_boundary = vision2force(vision)
# i,j = direction2neighbour(vector_attraction_to_boundary,3,3,5,metric,topology)
#



function remove_parts(space;c=25,r=0.6)
    L = size(space,1)
    space_removed = copy(space)
    N = sum(space .== 1)
    count = sum(space .== 1)
    gardefou = 1
    while count ≥ r*N && gardefou<50
        i,j = Tuple(rand(findall(isequal(1),space)))
        space_removed = circshift(space_removed,-1 .* (i-Int(L/2),j-Int(L/2)))
        space_removed[Int(L/2)-c:Int(L/2)+c , Int(L/2)-c:Int(L/2)+c] .= 0
        space_removed = circshift(space_removed,(i-Int(L/2),j-Int(L/2)))
        count = sum(space_removed .== 1)
        gardefou += 1
    end
    if gardefou ≥ 49 println("Ca mouline.") end
    return space_removed
end


function run_simu(Lbase,dilation,nbs_nucleation,metric,topology,R,times,N, pourcentage;remove=false,erasure=0.0)
    history_locations = Array{Tuple{Int,Int}}(undef,N,length(times),length(nbs_nucleation),length(pourcentage),R)
    Threads.@threads for r in 1:R
        for a in eachindex(nbs_nucleation) , b in eachindex(pourcentage)
            radius_vision = round(Int,10/sqrt(nbs_nucleation[a]))
            space = create_space(Lbase,k=dilation,metric=metric,nb_nucleation=nbs_nucleation[a],thicken=false)
            L = Lbase*dilation
            if remove space = remove_parts(space,c=25,r=pourcentage[b]) end
            loc = [Tuple(rand(1:L,2)) for n in 1:N]
            history_locations[:,1,a,b,r] = loc
            token = 2
            for t in 1:times[end]
                loc,space,~ = update(loc,radius_vision,space,L ; alpha=0.95,topology=topology,erasure=erasure)
                if t ≥ times[token]
                    history_locations[:,token,a,b,r] = loc
                    token = min(token+1,length(times))
                end

            end
        end
    end
    return history_locations
end

function shift_trajectory(locations,tup,L=200)
    for i in eachindex(locations)
        locations[i] = mod1.(locations[i] .+ tup,L)
    end
    return locations
end

function smooth(x,c=1)
    len = length(x)
    result = zeros(len)
    result[1:c] .= NaN
    result[end-c+1:end] .= NaN
    for i in 1+c:len-c
        result[i] = mean(x[i-c:i+c])
    end
    return result
end

## Self Avoiding Random Walks
# Source https://www.hiskp.uni-bonn.de/uploads/media/polymers.pdf
function SAW(Ntarget,params,method)
    if method == "BruteForce" || method == "BF"
        return SAW_brute_force(Ntarget)
    elseif method == "SimpleSampling" || method == "SS"
        return SAW_simple_sampling(Ntarget)
    elseif method == "BruteForcePersistence" || method == "BFPersistence"
        return SAW_brute_force_persistence(Ntarget,params)
    else println("The method should be chosen amongst 'SimpleSampling/SS' , 'BruteForce/BF' or 'BruteForcePersistence/BFPersistence'.")
    end
end

function SAW_brute_force_persistence(Ntarget,p)
    path = Tuple{Int, Int}[]
    push!(path, (0,0))
    push!(path, (0,1)) # choose first step arbitrarily

    displacement = (0,1)
    possibilities = [(+1, 0),(0, +1),(-1, 0),(0, -1)] # ordre = convention
    if displacement     == (1,0)  w = [p,(1-p)/2,0,(1-p)/2]
    elseif displacement == (0,1)  w = [(1-p)/2,p,(1-p)/2,0]
    elseif displacement == (-1,0) w = [0,(1-p)/2,p,(1-p)/2]
    elseif displacement == (0,-1) w = [(1-p)/2,0,(1-p)/2,p]
    end

    displacement = sample(possibilities,Weights(w))
    proposal = path[end] .+ displacement
    while !(proposal in path)
        push!(path,proposal) # accept proposal

        if displacement     == (1,0)  w = [p,(1-p)/2,0,(1-p)/2]
        elseif displacement == (0,1)  w = [(1-p)/2,p,(1-p)/2,0]
        elseif displacement == (-1,0) w = [0,(1-p)/2,p,(1-p)/2]
        elseif displacement == (0,-1) w = [(1-p)/2,0,(1-p)/2,p]
        end

        displacement = sample(possibilities,Weights(w))
        proposal = path[end] .+ displacement
    end
    return path
end

function SAW_brute_force(N::Integer)
    path = Tuple{Int, Int}[]
    push!(path, (0,0))

    possibilities = [(-1, 0),(+1, 0),(0, -1),(0, +1)]

    cur = path[end]
    proposal = cur .+ rand(possibilities)
    while !(proposal in path)
        cur = proposal
        push!(path,cur)
        proposal = cur .+ rand(possibilities)
    end
    return path
end

function SAW_stepback(N::Integer,maxfails)
    path = Tuple{Int, Int}[]
    push!(path, (0,0))

    max_stepsback = 15
    weights_backstep = [round(1/sqrt(i),digits=3) for i in 1:max_stepsback]
    weights_backstep = weights_backstep/sum(weights_backstep)

    i = 0 ; fails = 0
    while i < N && fails < maxfails
        x, y = path[end]
        valid = []

        for pos in [(x-1, y),(x+1, y),(x, y-1),(x, y+1)]
            if !(pos in path)
                push!(valid,pos)
            end
        end
        if isempty(valid)
            if fails < maxfails - 1
                n = min(sample(1:max_stepsback,Weights(weights_backstep)),length(path)-1)
                for j in 1:n pop!(path) end
                i -= n
            end
            fails += 1
        elseif length(valid) > 0
            push!(path, rand(valid))
            i += 1
        end
    end
    return path
end


function SAW_simple_sampling(N::Integer)
    path = Tuple{Int, Int}[]
    push!(path, (0,0))
    push!(path, (0,1))
    weights = [1.0,1.0]

    i = 0
    while i < N-2
        x, y = path[end]
        valid = []

        for pos in [(x-1, y),(x+1, y),(x, y-1),(x, y+1)]
            if !(pos in path)
                push!(valid,pos)
            end
        end
        possible = length(valid)
        if possible == 0
            break
        elseif possible > 0
            push!(path, rand(valid))
            push!(weights, weights[end]*possible/3)
            i += 1
        end
    end
    # path_vec = Vector{Tuple{Float64,Float64}}(undef,N)
    # path_vec[1:length(path)] = path
    # for i in 1+length(path):N path_vec[i] = (0,0) end
    #
    # weights_vec = Vector{Float64}(undef,N)
    # weights_vec[1:length(weights)] = weights
    # for i in 1+length(weights):N weights_vec[i] = 0 end

    return path,weights
end
