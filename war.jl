## Simple Sampling SAW (pas fini, ya des pb)
T = 2000
M = Int(1E5)
p = "dummy"
dx = NaN*zeros(Int,T,M)
dy = NaN*zeros(Int,T,M)
w  = NaN*zeros(T,M)
z  = @elapsed for m in 1:M
    walk,weights = SAW(T,p,"SS")
    duree = length(walk)
    dx[1:duree,m] = [walk[t][1] for t in 1:duree]
    dy[1:duree,m] = [walk[t][2] for t in 1:duree]
    w[1:duree,m] = weights
end
prinz(z)


m= 1
walk,weights = SAW(T,p,"SS")
duree = length(walk)
dx[1:duree] = [walk[t][1] for t in 1:duree]
dy[1:duree] = [walk[t][2] for t in 1:duree]
w[1:duree] = weights


times = collect(10:20:500)
xmax = maximum(abs,filter(!isnan,vcat(dx,dy)))
dx_abs = abs.(dx)
dy_abs = abs.(dy)
Gxt = NaN*zeros(length(1:1:xmax),length(times))
z = @elapsed for i in 1:1:Int(xmax), t in each(times)
    tmp = 0
    tt = times[t]
    for m in 1:M
        if dx_abs[tt,m] == i tmp += w[tt,m] end
        if dy_abs[tt,m] == i tmp += w[tt,m] end
    end
    tmp ≠ 0 ? Gxt[i,t] = tmp : Gxt[i,t] = NaN
end
prinz(z)
Gxt = Gxt/2M
plot(Gxt[:,2],axis=:log)
p=plot()
    for t in each(times)
        # plot!(Gxt[:,t])
        plot!(Gxt[:,t],axis=:log)
    end
    p

## Regarder MSD as a function of p
T = 1000
M = Int(1E6)
ps = collect(1/3:0.1:1)
SD = [Float64[] for m in 1:T,j in each(ps) ]

z = @elapsed for m in 1:M, j in each(ps)
    walk = SAW(T,ps[j],"BFPersistence")
    for t in 1:min(T,length(walk))
        push!(SD[t,j],walk[t]...)
    end
end
prinz(z)

MSD = [mean(abs2,SD[t,j]) for t in 1:T , j in each(ps)]
p=plot(legend=(0.45,0.01))
    plot!(collect(1:T),2remove_negative(MSD[:,1])*sqrt(log(1/ps[1])),axis=:log,label="p=1/3",rib=0)
    for j in 2:length(ps)
        # plot!(collect(1:T)*(log(1/ps[j])),remove_negative(MSD[:,j])*(log(1/ps[j]))^2,axis=:log,label="p=$(round(ps[j],digits=2))",rib=0)
        plot!(collect(1:T),2remove_negative(MSD[:,j])*sqrt(log(1/ps[j])),axis=:log,label="p=$(round(ps[j],digits=1))",rib=0)
    end
    xx = 1:500
    plot!(xx,1.3xx.^1.5,c=:black,line=:dash,lw=2)
    annotate!(6,20000,text("MSD/λ²"))
    annotate!(7,7000,text("λ ∼ 1/(log(1/p))",11,:center))
    annotate!(8,70,text(L"∼\,(t/λ)^{3/2}",11,:center,35.0))
    annotate!(650,0.03,text(L"t/λ"))
    annotate!(1,3E-2,text("(a)",14,:left))
    ylims!(0.01,60000)
    yticks!([0.1,1,10,100,1000,1E4])
# savefig("figures\\MSD_SAW_persistence.png")
# savefig("figures\\MSD_SAW_persistence.pdf")

using SpecialFunctions
p2=plot(yticks=[1E-5,1E-4,1E-3,1E-2,1E-1,1])
# p2=plot(xlabel=L"x^4/(λ^{1/4}\,t^3)",ylabel=L"G(x,t)⋅λ^{1/16}\,t^{3/4}",yticks=[1E-5,1E-4,1E-3,1E-2,1E-1,1])
    i = 30
    expo = 4
    for ind_p in each(ps)
        data = abs.(SD[i,ind_p])
        N = length(data)
        h = fit(Histogram,data,0:floor(Int,maximum(data)))
        r = Array(h.edges[1])
        plot!(r[2:end-1].^expo*log(1/ps[ind_p])/i^3,i^0.75/log(1/ps[ind_p])^(1/4)*smooth(remove_negative(h.weights[2:end]/N),1),rib=0,yaxis=:log)
    end
    xx= 0:0.01:8; plot!(xx,0.7exp.(-xx*0.75),c=:black,line=:dash,lw=2)
    xlims!(-0.2,10.5)
    ylims!(3E-4,3.5)
    annotate!(8.3,5.5E-4,text(L"x^4/(λ\,t^3)"))
    annotate!(0.1,1.8,text(L"G(x,t)⋅(λ\,t^{3})^{1/4}",:left))
    annotate!(0.05,5.5E-4,text("(b)",14,:left))
    p2

plot(p,p2,size=(800,400))
# savefig("figures\\G_MSD_SAW_persistence.png")
# savefig("figures\\G_MSD_SAW_persistence.pdf")

p2=plot(xlabel=L"x^4/t^3",ylabel=L"t^{3/4}\,G",yticks=[1E-5,1E-4,1E-3,1E-2,1E-1,1])
    expo = 2.65
    for i in 15:3:25
        data = abs.(SD[i,ind_p])
        N = length(data)
        h = fit(Histogram,data,0:floor(Int,maximum(data)))
        r = Array(h.edges[1])
        plot!(r[2:end-1].^expo/(i)^(expo*0.75),i^(5/7)*smooth(remove_negative(h.weights[2:end]/N),1),rib=0,yaxis=:log)
    end
    xx = 0:0.1:15
    # plot!(xx,2exp.(-0.8xx),c=:black,line=:dash)
    # annotate!(8,0.005,text(L"\sim f(y) = e^{-0.8\,y}",10,:left))
    p2
&
## Films détaillés des vortex/défauts pour Demian
#= 1. Creer une bonne realisation et la sauvegarder
2. Tracker le vortex
3. Left Panel of the movie : global view
3. Right Panel of the movie : local 20x20 view with arrows (and light colors ?) =#
# locsave = copy(history_locations)
# thetasave = copy(history_thetas)
# omegas_save = copy(omegas)

# history_locations[1:1001] = locsave
# history_thetas[:,:,1:1001] = thetasave
# thetas = thetasave[:,:,end]

# 1. Create Realisation
L = 300
    T = 0.02
    Var = 0.1
    init = "isolated"
    BC = "free"
    q = +1 ; r0 = 100
    tmax = 5000 ; transients = 500
    every = 1 ; times = collect(transients:every:tmax)

    history_locations = Vector{Union{Tuple{Int16,Int16},Missing}}(missing,length(times))
    history_thetas    = zeros(L,L,length(times))

thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
    t = 0.0
    while t < transients
        t += dt ; thetas = update_FBC(thetas,omegas,L,T,dt)
    end

last_location,~ = spot_single_default_global(thetas)
    history_locations[1] = last_location
    history_thetas[:,:,1] = thetas  # O(1)

    token = 2 ;
    while t < tmax
        t += dt
        thetas = update_FBC(thetas,omegas,L,T,dt)
        if t > times[token]
            try
                lastknown = (Int16(-1),Int16(-1)) # dummy because no -1 defect to be dealt with here
                last_location,alone = spot_single_default_local(thetas,last_location,lastknown,q)
                if alone # there are no other vortices in the surroundings (the '2margin' square)
                    saved_location = last_location
                else # there are other vortices in the surroundings that could pollute the data
                    saved_location = missing
                end
                history_locations[token] = saved_location  # O(1)
                history_thetas[:,:,token] = thetas  # O(1)
            catch e  # we lost our vortex (out of lattice, annihilation)
                println(e)
                printstyled("Warning : Vortex lost, simulation stopped at t = $(round(t,digits=2)). \n"; color = :yellow)
                break # we leave all the remaining values as they are : missing
            end
            token = min(token+1,length(times))
            print() ; println("t = $(times[token]) , Defect position = $(history_locations[token-1])")
        end
    end

# jldsave("data/trajs_Demian.jld2";history_locations,history_thetas,L,init,r0,tmax,transients,times,T,Var,BC,every,locsave,thetasave,omegas_save,anim)

p1 = heatmap(mod.(history_thetas[:,:,i]',2π),c=cols,colorbar=false)
    highlight_defects(p1,history_thetas[:,:,i])
    xlims!(100,L-50)
    ylims!(95,L-55)

window = 13
    posx,posy = history_locations[i]
    zoom_defect = history_thetas[posx-window:posx+window,posy-window:posy+window,i]
    p2 = heatmap(mod.(zoom_defect,2π)',c=cols,colorbar=false,alpha=1,axis=false,ticks=false)
    for j in 1:2window+1
        quiver!(j*ones(2window+1),collect(1:2window+1),
            quiver=(cos.(zoom_defect[j,:]),sin.(zoom_defect[j,:])),
            c=:white,lw=0.8)
    end
    xlims!(1,2window+1)
    ylims!(1,2window+1)
    p2

window = 4
        posx,posy = history_locations[i]
        zoom_defect = history_thetas[posx-window:posx+window,posy-window:posy+window,i]
        p3 = heatmap(mod.(zoom_defect,2π)',c=cols,colorbar=false,alpha=1,axis=false,ticks=false)
        for j in 1:2window+1
            quiver!(j*ones(2window+1),collect(1:2window+1),
                quiver=(cos.(zoom_defect[j,:]),sin.(zoom_defect[j,:])),
                c=:white,lw=1.4)
        end
        xlims!(1,2window+1)
        ylims!(1,2window+1)
        p3

plot(p1,p2,p3,layout=(1,3),size=(1200,400))


z = @elapsed anim = @animate for i in 350:1200#:length(times)
    p1 = heatmap(mod.(history_thetas[:,:,i]',2π),c=cols,colorbar=false)
        highlight_defects(p1,history_thetas[:,:,i])
        xlims!(100,L-50)
        ylims!(95,L-55)

    window = 13
        posx,posy = history_locations[i]
        zoom_defect = history_thetas[posx-window:posx+window,posy-window:posy+window,i]
        p2 = heatmap(mod.(zoom_defect,2π)',c=cols,colorbar=false,alpha=1,axis=false,ticks=false)
        for j in 1:2window+1
            quiver!(j*ones(2window+1),collect(1:2window+1),
                quiver=(cos.(zoom_defect[j,:]),sin.(zoom_defect[j,:])),
                c=:white,lw=0.8)
        end
        xlims!(1,2window+1)
        ylims!(1,2window+1)

    window = 4
        posx,posy = history_locations[i]
        zoom_defect = history_thetas[posx-window:posx+window,posy-window:posy+window,i]
        p3 = heatmap(mod.(zoom_defect,2π)',c=cols,colorbar=false,alpha=1,axis=false,ticks=false)
        for j in 1:2window+1
            quiver!(j*ones(2window+1),collect(1:2window+1),
                quiver=(cos.(zoom_defect[j,:]),sin.(zoom_defect[j,:])),
                c=:white,lw=1.4)
        end
        xlims!(1,2window+1)
        ylims!(1,2window+1)
        p3

    plot(p1,p2,p3,layout=(1,3),size=(1200,400))

end
prinz(z)
mp4(anim,"figures/films/trajs_Demian1_v2.mp4",fps=25)
# mp4(anim,"figures/films/trajs_Demian_short.mp4",fps=15)

## MSD of SAW
Ns = []
msds = []

z = @elapsed for i in 1:Int(3E5)
    saw = SAW(100000)
    push!(Ns,length(saw))
    push!(msds,sum(abs2,saw[end]))
end
prinz(z)
sp = sortperm(Ns)
    Ns = Ns[sp]
    msds = msds[sp]
    Ns_reduced = []
    msds_reduced = []
for i in minimum(Ns):maximum(Ns)
    ind = findall(isequal(i),Ns)
    if length(ind) > 0
        push!(Ns_reduced,Ns[ind[1]])
        push!(msds_reduced,mean(msds[ind]))
    end
end
plot(xlims=(1,3maximum(Ns)),axis=:log,xlabel="t",ylabel="MSD",xticks=[1,10,100,1000])
    scatter!(Ns,msds,ms=0.5,c=:grey)
    # plot!(Ns_reduced,msds_reduced,c=:green)
    plot!(Ns_reduced,smooth(msds_reduced,3),c=:red)
    # plot!(Ns_reduced,Ns_reduced.*sqrt.(log.(Ns_reduced)))
    plot!(Ns_reduced,0.22Ns_reduced.^1.5,c=:black)
# savefig("figures\\MSD_SAW.png")
&

## SAW with persistence
function SAW(N::Integer,maxfails=1000)
    @assert N >= 0
    # this will store the path we have traversed
    path = Tuple{Int, Int}[]
    cur = (0,0)
    push!(path, cur)
    i = 0 ; fails = 0
    while i < N && fails < maxfails
        x, y = path[end]
        # now check the four possible movement directions
        # store moves that are allowed in valid vector
        valid = Tuple{Int, Int}[]
        if !((x-1, y) in path)
            push!(valid, (x-1, y))
        end
        if !((x+1, y) in path)
            push!(valid, (x+1, y))
        end
        if !((x, y-1) in path)
            push!(valid, (x, y-1))
        end
        if !((x, y+1) in path)
            push!(valid, (x, y+1))
        end
        if isempty(valid)
            n = sample(1:3,Weights([0.8,0.1,0.1]))
            for j in 1:n pop!(path) end
            i -= n
            fails += 1
        else
            cur = rand(valid)
            push!(path, cur)
            i += 1
        end
    end
    if fails ≥ maxfails
        # println("Too many fails !")
    end
    return path
end

function SAW_persist(N::Integer,p = 2, maxfails=1000)
    @assert N >= 0
    # this will store the path we have traversed
    path = Tuple{Int, Int}[(0,0)]
    lastpos = (0,0)
    cur = rand([(0,1),(0,-1),(1,0),(-1,0)])
    push!(path, cur)
    i = 0 ; fails = 0
    while i < N-2 && fails < maxfails
        x, y = path[end]
        a,b = (x,y) .- path[end-1] # (a,b) = vector from lastpos to cur

        valid = [(x,y) .+ (a,b), (x,y) .+ (-b,a), (x,y) .+ (b,-a)]
        wei = [p,1,1]
        for i in 1:3
            if valid[i] in path wei[i] = 0 end
        end
        if sum(wei) == 0
            n = sample(1:3,Weights([0.8,0.1,0.1]))
            for j in 1:n pop!(path) end
            i -= n
            fails += 1
        else
            cur = sample(valid,Weights(wei / sum(wei)))
            push!(path, cur)
            i += 1
        end
    end
    if fails ≥ maxfails println("Too many fails !") end
    return path
end

tra = SAW(1000)
tra = SAW_persist(300,1.5)
    plot(tra,m=false,c=:black)

N = 1000
    aim = 3
    trajs = []
    for n in 1:aim
        succes = false
        traj_SAW = Vector{Tuple{Int64, Int64}}(undef,N+1)
        while succes == false
            traj_SAW = SAW_persist(N,2)
            if length(traj_SAW) ≥ N succes = true end
        end
        push!(trajs,traj_SAW)
        println("Success!")
    end
&
plot(trajs[2][1:3:end],m=false,c=:black)


## SM : irrelevant Parameters
Lbase = 200
    dilation = 2
    L = Lbase*dilation
    nbs_nucleation = [3,4,5,6]
    metric = ["manhattan","euclidian"]
    topology = "square"
    R = 60
    tmax = Int(3E3) ; every = 10 ; times = vcat(0:10,20:every:tmax)
    pourcentage = 0.6
    N = 5

history_locations = Array{Tuple{Int,Int}}(undef,2,N,length(times),length(nbs_nucleation),length(pourcentage),R)
z1 = @elapsed history_locations[1,:,:,:,:,:] = run_simu(Lbase,dilation,nbs_nucleation,metric[1],topology,R,times,N, pourcentage;remove=true,erasure=0.0)
z2 = @elapsed history_locations[2,:,:,:,:,:] = run_simu(Lbase,dilation,nbs_nucleation,metric[2],topology,R,times,N, pourcentage;remove=true,erasure=0.0)
prinz(z1+z2)

SD = zeros(2,N,length(times)-1,length(nbs_nucleation),length(pourcentage),R)
    for r in 1:R , i in 1:length(times)-1 , a in 1:length(nbs_nucleation), b in 1:length(pourcentage) , n in 1:N
        SD[1,n,i,a,b,r] = dist2(history_locations[1,n,i+1,a,b,r] , history_locations[1,n,1,a,b,r],L)
        SD[2,n,i,a,b,r] = dist2(history_locations[2,n,i+1,a,b,r] , history_locations[2,n,1,a,b,r],L)
    end
    MSD = mean(SD,dims=(2,6))

comp_metric = plot(legend=:bottomright)
    plot!(times[2:end],MSD[1,1,:,1,1,1],axis=:log,lw=1,line=:solid,rib=0,label="Metric = manhattan")
    plot!(times[2:end],MSD[2,1,:,1,1,1],axis=:log,lw=1,line=:solid,rib=0,label="Metric = euclidian")
    plot!(times[8:40],.42times[8:40].^1.5,axis=:log,c=:black)
    # plot!(times[8:40],.242times[8:40].^1.8,axis=:log,c=:black)

pmsd_manhattan = plot(legend=:topleft,xticks=[1,10,100,1000],yticks=[1,10,100,1000,1E4,1E5],xlabel="t",ylabel="MSD")
    for b in eachindex(nbs_nucleation)
        plot!(times[2:end],MSD[1,1,:,b,1,1],c=b,axis=:log,lw=1,line=:solid,label="n = $(nbs_nucleation[b])",rib=0)
    end
    # plot!(times[35:end],5.42times[35:end],axis=:log,c=:black,line=:dot)
    plot!(times[8:40],.242times[8:40].^1.5,axis=:log,c=:black)
    # plot!(times[5:150],0.15times[5:150].^2,c=:black,line=:dash)
    pmsd_manhattan

pmsd_euclidian = plot(legend=:topleft,xticks=[1,10,100,1000],yticks=[1,10,100,1000,1E4,1E5])
    for b in eachindex(nbs_nucleation)
        plot!(times[2:end],MSD[2,1,:,b,1,1],c=b,axis=:log,lw=1,line=:solid,label="n = $(nbs_nucleation[b])",rib=0)
    end
    # plot!(times[35:end],5.42times[35:end],axis=:log,c=:black,line=:dot)
    plot!(times[2:40],.542times[2:40].^1.5,axis=:log,c=:black)
    # plot!(times[5:150],0.15times[5:150].^2,c=:black,line=:dash)
    pmsd_euclidian

# now examples of boundary network
Lbase = 200
    dilation = 1
    L = Lbase*dilation
    nb_nucleation = 5
    topology = "square"

space_manhattan = create_space(Lbase,k=1,metric="manhattan",nb_nucleation=10,thicken=false)
net_manhattan   = heatmap(space_manhattan,c=cgrad([:black,:white]),colorbar=false,axis=false)
space_euclidian = create_space(Lbase,k=1,metric="euclidian",nb_nucleation=10,thicken=false)
net_euclidian   = heatmap(space_euclidian,c=cgrad([:black,:white]),colorbar=false,axis=false)

plot()
# savefig("figures/SM_irrelevance_nb_metric.p")


## Impact of the -1/2 on MSD
Lbase = 200
    dilation = 5
    L = Lbase*dilation
    nb_nucleation = [5]
    topology = "square"
    R = 200
    tmax = Int(2E3) ; every = 10 ; times = vcat(0:10,20:every:tmax)
    pourcentage = 0.4:0.1:0.8
    N = 1

z1 = @elapsed history_locations_no_minus12_scan_r = run_simu(Lbase,dilation, nb_nucleation,metric,topology,R,times,N,pourcentage,remove=true,erasure=0)
z2 = @elapsed history_locations_minus12_scan_r    = run_simu(Lbase,dilation, nb_nucleation,metric,topology,R,times,N,pourcentage,remove=true,erasure=-1/2)
prinz(z1+z2)

SD_no_minus12_scan_r = zeros(N,length(times)-1,length(nb_nucleation),length(pourcentage),R)
    SD_minus12_scan_r = zeros(N,length(times)-1,length(nb_nucleation),length(pourcentage),R)
    for r in 1:R , i in 1:length(times)-1 , a in 1:length(nb_nucleation), b in 1:length(pourcentage) , n in 1:N
        SD_no_minus12_scan_r[n,i,a,b,r] = dist2(history_locations_no_minus12_scan_r[n,i+1,a,b,r] , history_locations_no_minus12_scan_r[n,1,a,b,r],L)
        SD_minus12_scan_r[n,i,a,b,r] = dist2(history_locations_minus12_scan_r[n,i+1,a,b,r] , history_locations_minus12_scan_r[n,1,a,b,r],L)
    end
    MSD_no_minus12_scan_r = mean(SD_no_minus12_scan_r,dims=(1,5))
    MSD_minus12_scan_r    = mean(SD_minus12_scan_r,dims=(1,5))
p = plot(legend=:topleft,xticks=[1,10,100,1000],yticks=[1,10,100,1000,1E4,1E5])
    for b in [1,2,3,4]
        # plot!(times[2:200],MSD_no_minus12_scan_r[1,1:199,1,b,1],c=b,axis=:log,lw=1,line=:solid)
        plot!(times[2:end],MSD_minus12_scan_r[1,:,1,b,1],c=b,axis=:log,lw=1,line=:dash)
    end
    plot!(times[35:end],5.42times[35:end],axis=:log,c=:black,line=:dot)
    plot!(times[2:40],.542times[2:40].^1.5,axis=:log,c=:black)
    plot!(times[5:150],0.15times[5:150].^2,c=:black,line=:dash)
    p


## Figure to show that discontinuity and -1/2 are necessary
Lbase = 200
    dilation = 2
    L = Lbase*dilation
    nb_nucleation = [5]
    topology = "square"
    R = 200
    tmax = Int(2E3) ; every = 10 ; times = vcat(0:10,20:every:tmax)
    pourcentage = 0.6
    N = 1
    L^2/4

z = @elapsed history_locations_no_minus12 = run_simu(Lbase,dilation, nb_nucleation,metric,topology,R,times,N,pourcentage,remove=true)
prinz(z)

SD_no_minus12 = zeros(N,length(times)-1,length(nb_nucleation),length(pourcentage),R)
    for r in 1:R , i in 1:length(times)-1 , a in 1:length(nb_nucleation), b in 1:length(pourcentage) , n in 1:N
        SD_no_minus12[n,i,a,b,r] = dist2(history_locations_no_minus12[n,i+1,a,b,r] , history_locations_no_minus12[n,1,a,b,r],L)
    end
    MSD_no_minus12 = mean(SD_no_minus12,dims=(1,5))
p = plot(legend=:topleft,xticks=[1,10,100,1000],yticks=[1,10,100,1000,1E4,1E5])
    # plot!(times[2:end],MSD[1,:,1,1,1],axis=:log,lw=2)
    plot!(times[2:200],MSD_no_minus12[1,1:199,1,1,1],axis=:log,lw=2.5,c=3)
    plot!(times[2:200],MSD_not_discontinuous[1,1:199,1,1,1],axis=:log,lw=2)
    # plot!(times[15:end],10times[15:end],c=:black,line=:dot)
    plot!(times[3:55],.542times[3:55].^1.5,axis=:log,c=:black,lw=1)
    plot!(times[12:150],0.15times[12:150].^2,c=:black,line=:dash,lw=1)
    annotate!(2.5,1.2E5,text("MSD",15))
    annotate!(1E3,1.8,text("t",15))
    annotate!(200,29000,text(L"\sim t^2",12))
    annotate!(200,300,text(L"\sim t^{3/2}",12))
    xlims!(0.75,2300)
    p
# savefig("figures\\clean32_R200.pdf")

@unpack Lbase,dilation,L,R,times,tmax,topology,pourcentage,N,nb_nucleation,history_locations,MSD,SD,history_locations_not_discontinuous,MSD_not_discontinuous,SD_not_discontinuous,history_locations_no_minus12,MSD_no_minus12,SD_no_minus12 = load("data/hunt_for32/necessary_ingredients.jld2")
# jldsave("data/hunt_for32/necessary_ingredients.jld2";Lbase,dilation,L,R,times,tmax,topology,pourcentage,N,nb_nucleation,history_locations,MSD,SD,history_locations_not_discontinuous,MSD_not_discontinuous,SD_not_discontinuous,history_locations_no_minus12,MSD_no_minus12,SD_no_minus12)


## Figure to show boundary networks
# Toy model
Lbase = 200
    dilation = 1
    L = Lbase*dilation
    nb_nucleation = 5
    topology = "square"
    radius_vision = 5 # ::Int

space = create_space(Lbase,k=1,topology=topology,nb_nucleation=10,thicken=false)
p0 = heatmap(space,yflip = false,c=cgrad([:black,:white]),colorbar=false,axis=false)
space_removed = remove_parts(space,c=30,r=0.6)
    p3 = heatmap(space_removed,yflip = false,c=cgrad([:black,:white]),colorbar=false,axis=false)
    plot(p1,p2,p3,plot(ticks=false),layout=(2,2),size=(1600,1600))
plot(p0,p3,layout=(1,2),size=(800,400))
# savefig("figures\\network_2.png")


# Spin model
L = 200
    T = 0.05
    Var = 0.2
    init = "pair"
    BC = "periodic"
    q = +1 ; r0 = Int(L/2) # dummy
    tmax = 200

t = 0.0 ; thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
while t < 600 t += dt ; update(thetas,omegas,L,T,dt) end
p_spin_model = heatmap(mod.(thetas',2π),c=cols,colorbar=false,ticks=false,size=(1000,1000),alpha=0.2)
# savefig("figures\\spin_model_network.png")

## Question Demian : can continous network & -1/2 & many defects help retrieve 3/2 ?
Lbase = 200
    dilation = 2
    L = Lbase*dilation
    nb_nucleation = [5]
    topology = "square"
    R = 10
    tmax = Int(2E3) ; every = 10 ; times = vcat(0:10,20:every:tmax)
    pourcentage = 0.6
    N = 20
    L^2/4

z = @elapsed history_locations = run_simu(Lbase,dilation, nb_nucleation,topology,R,times,N,pourcentage,remove=false,erasure=-0.5)
prinz(z)

SD = zeros(N,length(times)-1,length(nb_nucleation),length(pourcentage),R)
    for r in 1:R , i in 1:length(times)-1 , a in 1:length(nb_nucleation), b in 1:length(pourcentage) , n in 1:N
        SD[n,i,a,b,r] = dist2(history_locations[n,i+1,a,b,r] , history_locations[n,1,a,b,r],L)
    end
    MSD = mean(SD,dims=(1,5))
p = plot(legend=:topleft,xticks=[1,10,100,1000],ylabel="MSD",xlabel="t")
    for b in eachindex(pourcentage)
        plot!(times[2:end],MSD[1,:,1,b,1],axis=:log,lw=3)
    end
    # plot!(times[15:end],10times[15:end],c=:black,line=:dot)
    plot!(times[2:40],.542times[2:40].^1.5,axis=:log,c=:black)
    plot!(times[2:80],0.15times[2:80].^2,c=:black,line=:dash)
    # annotate!(2.5,1.2E4,text("MSD",15))
    # annotate!(1E3,0.8,text("t",15))
    annotate!(200,35000,text(L"\sim t^2",12))
    annotate!(50,30,text(L"\sim t^{3/2}",12))
    # xlims!(0.75,2300)
    p


## Sample trajectories
# Levy Walks
using Distributions

function LW(T,D = Pareto())
    t = 0.0 ;  v = 1/2
    traj = [(0.0,0.0)]
    times = [0.0]
    while t < T
        τ = min(T,rand(D)) # time of the next walk
        θ = rand()*2π # direction of the next walk
        next_loc = traj[end] .+ (v*τ) .*(cos(θ),sin(θ))
        push!(traj,next_loc)
        t += τ
        push!(times,t)
    end
    return traj,times
end


function LWNonBack(T,D,epsilon)
    t = 0.0 ;  v = 1
    traj = [(0.0,0.0)] ; θprevious = 0.0
    times = [0.0]

    while t < T
        τ = min(T,rand(D)) # duration of the next walk
        θ = θprevious + epsilon - π + rand()*(2π-2epsilon) # direction of the next walk
        next_loc = traj[end] .+ (v*τ) .*(cos(θ),sin(θ)) ; push!(traj,next_loc)
        t += τ
        θprevious = θ
        push!(times,t)

    end
    return traj,times
end

# traj_levy,timesLW = LW(1000, Pareto(1.5,2))
traj_levy2,timesLW2 = LW(1000, Pareto(1.5,2))
# traj_levyNB,timesLWNB = LWNonBack(1000, Pareto(1.5,2),π/3)

p3 = plot(legend=:bottomleft,xlims=(0,200),ylims=(0,200))
    plot!([traj_levy[i] .* 0.5 .+ (100,140) for i in eachindex(traj_levy)],line_z=timesLW,c=cgrad([:blue,:green,:orange,:red]),colorbar=false)
    plot!([traj_levy2[i] .* 0.5 .+ (80,30) for i in eachindex(traj_levy2)],line_z=timesLW2,c=cgrad([:blue,:green,:orange,:red]),colorbar=false)
    plot!([traj_levyNB[i] .* 0.5  .+ (130,80) for i in eachindex(traj_levyNB)],line_z=timesLWNB,c=cgrad([:blue,:green,:orange,:red]),colorbar=false)
    scatter!((100,140),c=:black,ms=5)
    scatter!((130,80),c=:black,ms=5)
    scatter!((80,30),c=:black,ms=5)
    annotate!(121,79,text("5",6))
    annotate!(100.5,148,text("6",6))
    annotate!(80,40,text("7",6))

# jldsave("data/trajectories_Levy_figure.jld2";traj_levy,traj_levy2,traj_levyNB,timesLW,timesLW2,timesLWNB,p3)
@unpack traj_levy,traj_levy2,traj_levyNB,timesLW,timesLW2,timesLWNB,p3 = load("data/trajectories_Levy_figure.jld2")

# Trajectories of the mini model
Lbase = 200
    dilation = 1
    L = Lbase*dilation
    nb_nucleation = [6]
    topology = "square"
    tmax = Int(1E3) ; every = 5 ; times = vcat(0:every:tmax)
    pourcentage = 0.6
    N = 1
    R = 10

history_loc = run_simu(Lbase,dilation, nb_nucleation,topology,R,times,N,pourcentage,remove=true,erasure=0)
history_locminus12 = run_simu(Lbase,dilation, nb_nucleation,topology,R,times,N,pourcentage,remove=true,erasure=-0.5)


# rs_cool_minus12 = [2,3,5,6,10]
# shifts_minus12 = [(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),]


rs_cool = [1,4,8]
shifts = [(20,-40),(30,-60),(90,0)]
    p4 = plot(xlims=(0,L),ylims=(0,L))
    for r in 1:3
        plot!(shift_trajectory(history_loc[1,:,1,1,rs_cool[r]],shifts[r]),line_z=times,c=cgrad([:blue,:green,:orange,:red]),colorbar=false) ; scatter!(history_loc[1,1,1,1,rs_cool[r]] .+ shifts[r],c=:black,m=:circle,ms=5)
        # plot!(shift_trajectory(history_locminus12[1,:,1,1,r],shifts_minus12[r]),line_z=times,c=cgrad([:blue,:green,:orange,:red])) ; scatter!(history_locminus12[1,1,1,1,r] .+ shifts_minus12[r],c=:black,m=:utriangle,ms=8)
    end
    annotate!(34,158,text("13",6))
    annotate!(84,80,text("12",6))
    annotate!(126,135,text("11",6))
    p4
# savefig(pathfig*"sample_traj_model.svg")

# jldsave("data/trajectories_minimodel_figure.jld2";history_loc,rs_cool,shifts,p4,L,nb_nucleation,pourcentage,N,R,tmax,every,times,comment="Erasure = 0 and not -1/2")
@unpack history_loc,rs_cool,shifts,p4,L,nb_nucleation,pourcentage,N,R,tmax,every,times = load("data/trajectories_minimodel_figure.jld2")

# Trajectories of the actual model
L = 200
    T = 0.1
    Var = 0.1
    init = "isolated"
    BC = "free"
    q = +1 ; r0 = 30
    transients = 200 ; tmax = transients + 1000
    every = 10 ; times = collect(transients:every:tmax)
    R = 30

history_locations = Matrix{Union{Tuple{Int16,Int16},Missing}}(missing,length(times),R)
# history_locations = Vector{Union{Tuple{Int16,Int16},Missing}}(missing,length(times))
lastknown = (Int16(-1),Int16(-1)) # dummy because no -1 defect to be dealt with here

z = @elapsed Threads.@threads for r in 1:R
    thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
    t = 0.0
    while t < transients
        t += dt ; thetas = update(thetas,omegas,L,T,dt)
    end

    try
        last_location,~ = spot_single_default_local(thetas,(Int16(L/2),Int16(L/2)),lastknown,q,round(Int,L/4))
        history_locations[1,r] = last_location
    catch
        println("Vortex not detected au tout debut.")
        break
    end
    token = 2 ;
    while t < tmax
        t += dt ; thetas = update(thetas,omegas,L,T,dt)
        if t > times[token]
            if r == 1 print() ; println("t = $(times[token])") end
            try
                last_location,alone = spot_single_default_local(thetas,history_locations[token-1,r],lastknown,q)
                if alone # there are no other vortices in the surroundings (the '2margin' square)
                    saved_location = last_location
                else # there are other vortices in the surroundings that could pollute the data
                    saved_location = missing
                end
                history_locations[token,r] = saved_location  # O(1)
            catch e  # we lost our vortex (out of lattice, annihilation)
                println(e)
                printstyled("Warning : Vortex lost, simulation r$r stopped at t = $(round(t,digits=2)). \n"; color = :yellow)
                break # we leave all the remaining values as they are : missing
            end
            token = min(token+1,length(times))
            # print() ; println("Defect position = $(history_locations[token-1,r])")
        end
    end
end
prinz(z)

history_locations_plot = Matrix{Tuple{Float64,Float64}}(undef,length(times),R)
for i in eachindex(history_locations_plot)
    if history_locations[i] === missing
        history_locations_plot[i] = (NaN,NaN)
    else
        history_locations_plot[i] = history_locations[i]
    end
end

# jldsave("data/trajectories_actual_model_T$(T)_Var$(Var).jld2";L,T,Var,init,times,tmax,transients,R,history_locations_plot,history_locations,BC)
@unpack L,Var,init,times,tmax,transients,R,history_locations_plot,history_locations = load("data/trajectories_actual_model_T0.1_Var0.1.jld2")

# for r in [4,5,10,18,20,26]
#     p = plot(size=(1000,800),show=true,title="r = $r")
#         plot!(history_locations_plot[:,r],xlims=(50,125),ylims=(50,150),line_z=times.-transients,c=cgrad([:blue,:green,:orange,:red]))
#         display(p)
#         sleep(0.1)
#     end

rs_cool = [4,5,10,18,20,26]
shifts = [(-20,+25),(+35,0),(0,-20),(+25,-10),(20,+40),(-10,20)]
p1 = plot(xlims=(50,150),ylims=(50,150))
    for r in [1,2,3,6]
        plot!(shift_trajectory(history_locations_plot[:,rs_cool[r]],shifts[r]),line_z=times.-transients,c=cgrad([:blue,:green,:orange,:red]),colorbar=false)
        scatter!(history_locations_plot[1,rs_cool[r]] .+ shifts[r],c=:black,m=:circle)
    end
    annotate!(82.2,113.5,text(L"1",:black,6))
    annotate!(134.5,97.5,text(L"4",:black,6))
    annotate!(105,87,text(L"3",:black,6))
    annotate!(96.5,100,text(L"2",:black,6))

# p2 = plot(xlims=(50,150),ylims=(50,150))
#     for r in 4:6
#         plot!(shift_trajectory(history_locations_plot[:,rs_cool[r]],shifts[r]),line_z=times.-transients,c=cgrad([:blue,:green,:orange,:red]),colorbar=false)
#         scatter!(history_locations_plot[1,rs_cool[r]] .+ shifts[r],c=:black,m=:circle)
#     end
#     annotate!(130,82,text(L"4",:black,6))
#     annotate!(122.5,140.7,text(L"5",:black,6))

# jldsave("data/trajectories_actualdefects_figure.jld2";history_locations_plot,history_locations,rs_cool,shifts,p1,p2,nb_nucleation,pourcentage,N,every,L,T,Var,init,times,tmax,transients,R,BC)

# SAW
function SAW(N::Integer,maxfails=1000)
    @assert N >= 0
    # this will store the path we have traversed
    path = Tuple{Int, Int}[]
    cur = (0,0)
    push!(path, cur)
    i = 0 ; fails = 0
    while i < N && fails < maxfails
        x, y = path[end]
        # now check the four possible movement directions
        # store moves that are allowed in valid vector
        valid = Tuple{Int, Int}[]
        if !((x-1, y) in path)
            push!(valid, (x-1, y))
        end
        if !((x+1, y) in path)
            push!(valid, (x+1, y))
        end
        if !((x, y-1) in path)
            push!(valid, (x, y-1))
        end
        if !((x, y+1) in path)
            push!(valid, (x, y+1))
        end
        if isempty(valid)
            n = sample(1:3,Weights([0.8,0.1,0.1]))
            for j in 1:n pop!(path) end
            i -= n
            fails += 1
        else
            cur = rand(valid)
            push!(path, cur)
            i += 1
        end
    end
    if fails ≥ maxfails println("Too many fails !") end
    return path
end

function SAW_persist(N::Integer,maxfails=1000)
    @assert N >= 0
    # this will store the path we have traversed
    path = Tuple{Int, Int}[(0,0)]
    lastpos = (0,0)
    cur = rand([(0,1),(0,-1),(1,0),(-1,0)])
    push!(path, cur)
    i = 0 ; fails = 0
    while i < N && fails < maxfails
        x, y = path[end]
        a,b = cur .- path[end-1] # (a,b) = vector from lastpos to cur

        valid = [cur .+ (a,b), cur .+ (-b,a), cur .+ (b,-a)]
        wei = [2,1,1]
        for i in 1:3
            if valid[i] in path wei[i] = 0 end
        end
        wei = wei / sum(wei)
        if isempty(valid)
            n = sample(1:3,Weights([0.8,0.1,0.1]))
            for j in 1:n pop!(path) end
            i -= n
            fails += 1
        else
            cur = sample(valid,Weights(wei))
            push!(path, cur)
            i += 1
        end
    end
    if fails ≥ maxfails println("Too many fails !") end
    return path
end

N = 1000
    aim = 3
    trajs = []
    for n in 1:aim
        succes = false
        traj_SAW = Vector{Tuple{Int64, Int64}}(undef,N+1)
        while succes == false
            traj_SAW = SAW(N)
            if length(traj_SAW) ≥ N succes = true end
        end
        push!(trajs,traj_SAW)
    end

trajs_modified = copy(trajs)
    trajs_modified[1] = [trajs[1][j] .+ (-50,140) for j in 1:1001]
    trajs_modified[2] = [trajs[2][j] .+ (80,40) for j in 1:1001]
    trajs_modified[3] = [trajs[3][j] .+ (40,25) for j in 1:1001]
    p2 = plot(xlims=(-120,150),ylims=(-60,160))
    for i in 1:3
        plot!(trajs_modified[i][1:10:end-1],line_z=collect(0:10:1000),c=cgrad([:blue,:green,:orange,:red]),colorbar=false)
        scatter!(trajs_modified[i][1],m=:circle,ms=5,c=:black)
    end
    annotate!(51,25,text(L"8",:black,6))
    annotate!(68,39,text(L"9",:black,6))
    annotate!(-50,150,text(L"10",:black,6))
    p2

# Final Figure
plot(p1,p3,p2,p4,layout=(2,2),size=(800,800))
# savefig(pathfig*"trajectories.svg")

# Time legend
pl = plot(xlims=(50,150),ylims=(50,150))
    for r in 1
        plot!(shift_trajectory(history_locations_plot[:,rs_cool[r]],shifts[r]),line_z=times,c=cgrad([:blue,:green,:orange,:red]))
    end
    pl
# savefig(pathfig*"trajectories_colorbar.svg")

## Clean Figure 3/2
Lbase = 200
    dilation = 5
    L = Lbase*dilation
    nb_nucleation = [5]
    topology = "square"
    R = 200
    tmax = Int(2E3) ; every = 10 ; times = vcat(0:10,20:every:tmax)
    pourcentage = 0.6
    N = 1
    L^2/4

z = @elapsed history_locations = run_simu(Lbase,dilation, nb_nucleation,topology,R,times,N,pourcentage,remove=true)
prinz(z)

SD = zeros(N,length(times)-1,length(nb_nucleation),length(pourcentage),R)
    for r in 1:R , i in 1:length(times)-1 , a in 1:length(nb_nucleation), b in 1:length(pourcentage) , n in 1:N
        SD[n,i,a,b,r] = dist2(history_locations[n,i+1,a,b,r] , history_locations[n,1,a,b,r],L)
    end
    MSD = mean(SD,dims=(1,5))
p = plot(legend=:topleft,xticks=[1,10,100,1000])
    for b in eachindex(pourcentage)
        plot!(times[2:end],MSD[1,:,1,b,1],axis=:log,lw=3)
        # plot!(times[2:end],MSD_not_discontinuous[1,:,1,b,1],axis=:log,lw=1.3)
    end
    plot!(times[15:end],10times[15:end],c=:black,line=:dot)
    plot!(times[2:40],.542times[2:40].^1.5,axis=:log,c=:black)
    # plot!(times[2:end],0.2times[2:end].^2,c=:black,line=:dot)
    annotate!(2.5,1.2E4,text("MSD",15))
    annotate!(1E3,0.8,text("t",15))
    annotate!(200,10000,text(L"\sim t",12))
    annotate!(50,30,text(L"\sim t^{3/2}",12))
    xlims!(0.75,2300)
    p
# savefig("figures\\clean32_R200.svg")

## Timing creation of space
using BenchmarkTools
function create_space_before(L;k=1,topology="square",nb_nucleation=Int(round(log(L),RoundUp)),thicken=false)
    grid,~ = init_grid(L,nb_nucleation)
    while minimum(grid) == 0
        grid = step(grid,L,topology)
    end
    return DB(grid,L,k,thicken)
end

function create_space_after(L;k=1,topology="square",nb_nucleation=Int(round(log(L),RoundUp)),thicken=false)
    grid,nucleation_sites = init_grid(L,nb_nucleation)
    for j in 1:L , i in 1:L
        if grid[i,j] == 0
            distances_to_nucl = [manhattan_dist((i,j),nuc,L) for nuc in nucleation_sites]
            grid[i,j] = argmin(distances_to_nucl)
        end
    end
    return DB(grid,L,k,thicken)
end

space = create_space_after(200,k=1,topology="square",nb_nucleation=6,thicken=false)
    heatmap(space,yflip = true,size=(450,400))

@btime create_space_after(200,k=1,topology="square",nb_nucleation=6,thicken=false)
@btime create_space_before(200,k=1,topology="square",nb_nucleation=6,thicken=false)
@btime begin
    space = create_space_after(200,k=1,topology="square",nb_nucleation=6,thicken=false)
    loc = [Tuple(rand(1:L,2)) for n in 1:N]
    for i in 1:1000
        loc,space,~ = update(loc,2,space,200 ; alpha=0.95,topology="square")
    end
end
## Total length as a function of nb_nucleation
nb = 2:2:20
Lbase = 200
R = 50
total_length = zeros(length(nb),R)
Threads.@threads for r in 1:R
    for i in eachindex(nb)
        total_length[i,r] = sum( create_space(Lbase,k=1,topology="square",nb_nucleation=nb[i],thicken=false) .== 1 )
    end
end
total_length_avg = mean(total_length,dims=2)
total_length_std = std(total_length,dims=2)
plot(nb,total_length_avg,m=true,rib = total_length_std)
    plot!(nb,380(nb).^0.5,c=:black)


radii_vision = round.(Int,10 ./sqrt.(nb))
plot(nb,radii_vision)

## Emulation of domains and domain boundaries
# Get MSD
Lbase = 200
    dilation = 3
    L = Lbase*dilation
    nb_nucleation = [5]
    topology = "square"
    R = 50
    tmax = Int(2E3) ; every = 10 ; times = vcat(0:10,20:every:tmax)
    pourcentage = 0.4:0.1:1
    N = 4

    L^2/4
z = @elapsed history_locations = run_simu(Lbase,dilation, nb_nucleation,topology,R,times,N,pourcentage,remove=true)
prinz(z)

SD = zeros(N,length(times)-1,length(nb_nucleation),length(pourcentage),R)
    for r in 1:R , i in 1:length(times)-1 , a in 1:length(nb_nucleation), b in 1:length(pourcentage) , n in 1:N
        SD[n,i,a,b,r] = dist2(history_locations[n,i+1,a,b,r] , history_locations[n,1,a,b,r],L)
    end
    MSD = mean(SD,dims=(1,5))
p = plot(legend=:topleft,xlabel="t",ylabel="MSD")
    for b in eachindex(pourcentage)
        plot!(times[2:end],pourcentage[b]*MSD[1,:,1,b,1],axis=:log,label="r = $(pourcentage[b])",rib=0)
    end
    # plot!(times[2:end],10times[2:end],c=:black,line=:dash)
    # plot!(times[2:end],.52times[2:end].^1.5,axis=:log,c=:black)
    # plot!(times[2:end],0.1times[2:end].^2,c=:black,line=:dot)
    p
# savefig("figures\\explanation_MSD_different_alphas_tmax1E5.pdf")

## Movie time !
gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5)
Lbase = 200
    dilation = 3
    L = Lbase*dilation
    nb_nucleation = 5
    topology = "square"
    radius_vision = 5 # ::Int
    tmax = Int(1E3) ; every = 10
    alpha = 0.95
    N = 4
    space = create_space(Lbase,k=dilation,topology=topology,nb_nucleation=nb_nucleation,thicken=false)
    space = remove_parts(space)
    global loc = [Tuple(rand(1:L,2)) for n in 1:N]
    p = heatmap(space,yflip = true)
    for n in 1:N
        scatter!(reverse(loc[n]),m=:circle,c=:red,xlims=(1,L),ylims=(1,L),size=3 .*(450,400),ms=16)
    end
    p


z = @elapsed anim = @animate for i in 1:every:tmax
    for n in 1:every
        global loc,space,col = update(loc,radius_vision,space,L ; alpha=alpha,topology=topology)
    end
    p = heatmap(space,yflip = true)
    for n in 1:N
        scatter!(reverse(loc[n]),m=:circle,c=col[n],xlims=(1,L),ylims=(1,L),size=3 .*(450,400),ms=16)
    end
    p
end
prinz(z)
mp4(anim,"figures/films/hunt_for_32/model_$(N)defects_tmax$(tmax)_alpha$(alpha)_L$(L)_removeparts_A.mp4",fps = 10)

# Enlever des bouts de réseau dans l'espoir de supprimer ce parametre arbitraire alpha ~ 0.6
space = create_space(400,k=1,topology=topology,nb_nucleation=20,thicken=false)
    heatmap(space,yflip = true,size=(450,400))
heatmap(remove_parts(space),yflip = true,size=(450,400))

## influence on MSD (pas ouf)
Lbase = 200
    dilation = 3
    L = Lbase*dilation
    nb_nucleation = Int(round(log(Lbase),RoundUp))
    nb_nucleation = 5
    topology = "square"
    radius_vision = 2 # ::Int
    R = 100
    tmax = Int(1E4) ; every = 25 ; times = 0:every:tmax
    alphas = [0.9]
    history_locations = Array{Tuple{Int,Int}}(undef,2,length(times),length(alphas),R)

z = @elapsed history_locations[1,:,:,:] = run_simu(Lbase,dilation, nb_nucleation,topology,R,times,radius_vision,alphas)
z = @elapsed history_locations[2,:,:,:] = run_simu(Lbase,dilation, nb_nucleation,topology,R,times,radius_vision,alphas,true)
prinz(2z)

SD = zeros(2,length(times)-1,length(alphas),R)
    for r in 1:R , i in 1:length(times)-1 , a in 1:length(alphas)
        SD[1,i,a,r] = dist2(history_locations[1,i+1,a,r] , (Int(L/2),Int(L/2)),L)
        SD[2,i,a,r] = dist2(history_locations[2,i+1,a,r] , (Int(L/2),Int(L/2)),L)
    end
    MSD = mean(SD,dims=4)

#p = plot(legend=:topleft,xlabel="t",ylabel="MSD")
    for a in eachindex(alphas)
        plot!(times[2:end],MSD[1,:,a,1],c=a,axis=:log,label="α = $(alphas[a])",rib=0)
        plot!(times[2:end],MSD[2,:,a,1],c=a,axis=:log,label="α = $(alphas[a])",rib=0,ls=:dash)
        # plot!(times[2:end],MSD[:,a],axis=:log,label="α = $(alphas[a])",rib=0)
    end
    plot!(times[2:end],times[2:end],c=:black,line=:dash)
    plot!(times[2:end],0.4times[2:end].^1.5,axis=:log,c=:black)
    plot!(times[2:end],0.1times[2:end].^2,c=:black,line=:dot)
&

## Test remove parts of the network
Lbase = 200
    dilation = 3
    L = Lbase*dilation
    nb_nucleation = Int(round(log(L),RoundDown))
    topology = "square"
    radius_vision = 5 # ::Int
    R = 10
    tmax = Int(1E3) ; every = 10 ; times = 0:every:tmax
    alphas = 0.99

space = create_space(Lbase,k=2,topology=topology,nb_nucleation=5,thicken=false)
    heatmap(space,yflip = true,size=(450,400))
heatmap(remove_parts(space,c=25,r=0.6),yflip = true,size=(450,400))

# Get trajectories from my model
plot(xlims=(0,600),ylims=(0,600))
    scatter!(history_locations[2,:,1,rand(1:100)],line=false,ms=3)

## Movies of actual trajectories
L = 200
    T = 0.1
    Var = 0.1
    init = "pair"
    BC = "periodic"
    q = +1 ; r0 = Int(L/2) # dummy
    transients = 200 ; tmax = 10000
    times = transients:5:tmax
    R = 15
    save_thetas = zeros(L,L,length(times),R)

z = @elapsed for r in 1:R
    t = 0.0 ; thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
        p = display_thetas(thetas) ; highlight_defects(p,thetas)
    while t < transients
            t += dt ; update(thetas,omegas,L,T,dt)
        end
    # p = display_thetas(thetas) ; highlight_defects(p,thetas)

    token = 1
    while t < tmax
        t += dt ; update(thetas,omegas,L,T,dt)
        if t > times[token]
            save_thetas[:,:,token,r] = thetas
            token = min(length(times),1+token)
        end
    end
    # p = display_thetas(thetas) ; highlight_defects(p,thetas)

    anim = @animate for i in 1:length(times)
        p=display_thetas(save_thetas[:,:,i,r]) ; ; highlight_defects(p,save_thetas[:,:,i,r])
    end
    mp4(anim,"figures/films/hunt_for_32/trajectories_defect_T$(T)_Var$(Var)_r$(r+10).mp4",fps=10);
end

prinz(z)
save_thetas = Float16.(save_thetas)
@summary save_thetas
jldsave("data/")


## Defect velocity
@unpack L,TV,init,BC,tmax,times,history_locations,R = load("data\\pushMSD.jld")
history_locations = history_locations[1,:,:]
function get_v(x,c)
    T,R = size(x)
    v = zeros(2,T-c,R)
    for r in 1:R , t in 1:T-c
        if ismissing(x[t+c,r]) || ismissing(x[t,r]) v[1,t,r] = v[2,t,r] = NaN
        else
            v[1,t,r] = (x[t+c,r] .- x[t,r])[1]
            v[2,t,r] = (x[t+c,r] .- x[t,r])[2]
        end
    end
    return v
end
v = get_v(history_locations,200)
    histogram(vec(v[1,:,:]),yaxis=:log)
    histogram!(vec(v[2,:,:]),yaxis=:log)

## Sample trajectories and comparaison to Levy Flights and SAW
# NKM
L = 200
T = 0.1
Var = 0.1
init = "isolated"
BC = "free"
q = +1 ; r0 = 30
tmax = 500 ; transients = 50
every = 10 ; times = collect(transients:every:tmax)
R = 10

history_locations = Matrix{Union{Tuple{Int16,Int16},Missing}}(missing,length(times),R)
# history_locations = Vector{Union{Tuple{Int16,Int16},Missing}}(missing,length(times))
lastknown = (Int16(-1),Int16(-1)) # dummy because no -1 defect to be dealt with here

Threads.@threads for r in 1:R
    thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
    t = 0.0
    while t < transients
        t += dt ; thetas = update(thetas,omegas,L,T,dt)
    end

    last_location,~ = spot_single_default_local(thetas,(Int16(L/2),Int16(L/2)),lastknown,q)
    history_locations[1,r] = last_location
    token = 2 ;
    while t < tmax
        t += dt ; thetas = update(thetas,omegas,L,T,dt)
        if t > times[token]
            if r == 1 print() ; println("t = $(times[token])") end
            try
                last_location,alone = spot_single_default_local(thetas,last_location,lastknown,q)
                if alone # there are no other vortices in the surroundings (the '2margin' square)
                    saved_location = last_location
                else # there are other vortices in the surroundings that could pollute the data
                    saved_location = missing
                end
                history_locations[token,r] = saved_location  # O(1)
            catch e  # we lost our vortex (out of lattice, annihilation)
                println(e)
                printstyled("Warning : Vortex lost, simulation r$r stopped at t = $(round(t,digits=2)). \n"; color = :yellow)
                break # we leave all the remaining values as they are : missing
            end
            token = min(token+1,length(times))
            # print() ; println("Defect position = $(history_locations[token-1,r])")
        end
    end
end

history_locations_plot = Matrix{Tuple{Float64,Float64}}(undef,length(times),R)

for i in eachindex(history_locations_plot)
    if history_locations[i] === missing
        history_locations_plot[i] = (NaN,NaN)
    else
        history_locations_plot[i] = history_locations[i]
    end
end

plot(history_locations_plot[:,1],xlims=(0,L),ylims=(0,L))


function pplus(a::Vector{Tuple{Float64, Float64}}, b::Tuple{Int64, Int64}) return [a[i] .+ b for i in eachindex(a)] end
function pplus(a::Vector{Tuple{Int64, Int64}}, b::Tuple{Int64, Int64}) return [a[i] .+ b for i in eachindex(a)] end
good_real = [1,2,5,6,7,9]
deltas = [(2,-2),(13,10),(-5,-5),(0,+17),(-10,0),(-8,0)]
delim = (75,125)
tstart = 1
p=plot(xlims=delim,ylims=delim)
for i in eachindex(deltas)
    plot!(pplus(history_locations_plot[tstart:end,good_real[i]] , deltas[i]),c=i) ; scatter!(history_locations_plot[tstart,good_real[i]] .+ deltas[i],c=i)
    # plot!(pplus(history_locations_plot[:,good_real[i]] , (0,0)),c=i) ; scatter!(history_locations_plot[1,good_real[i]] .+ (0,0),c=i)
end
p

# savefig(pathfig*"traj_NKM.svg")
# savefig(pathfig*"traj_NKM.pdf")
# jldsave("data\\traj_NKM.jld2";L,T,Var,transients,tmax,every,R,times,history_locations,history_locations_plot)
@unpack L,T,Var,transients,tmax,every,R,times,history_locations,history_locations_plot = JLD2.load("data\\traj_NKM.jld2")

plot(history_locations_plot)

# Sample trajectories for SAW
function SAW(L::Integer, N::Integer)
    @assert L >= 1
    @assert N >= 0
    # this will store the path we have traversed
    path = Tuple{Int, Int}[]
    cur = ((L+1) ÷ 2, (L+1) ÷ 2)
    push!(path, cur)
    for i in 1:N
        x, y = cur
        # now check the four possible movement directions
        # store moves that are allowed in valid vector
        valid = Tuple{Int, Int}[]
        if x > 1 && !((x-1, y) in path)
            push!(valid, (x-1, y))
        end
        if x < L && !((x+1, y) in path)
            push!(valid, (x+1, y))
        end
        if y > 1 && !((x, y-1) in path)
            push!(valid, (x, y-1))
        end
        if y < L && !((x, y+1) in path)
            push!(valid, (x, y+1))
        end
        if isempty(valid)
            # we are stuck and failed to generate a walk of length N
            # no moves are possible and we still have not made N moves
            return path
        end
        # randomly pick a move from valid moves
        cur = rand(valid)
        push!(path, cur)
    end
    # we have successfully generated the walk - return it
    return path
end
N = 200

traj_SAW = SAW(10000,N)
plot!(pplus(traj_SAW,(-5000,-5000)))

# Levy Walks
using Distributions
D = Pareto()
rand(D)

function LW(T,D = Pareto())
    t = 0.0 ;  v = 1
    traj = [(0.0,0.0)]
    while t < T
        τ = min(T,rand(D)) # time of the next walk
        θ = rand()*2π # direction of the next walk
        next_loc = traj[end] .+ (v*τ) .*(cos(θ),sin(θ))
        push!(traj,next_loc)
        t += τ
    end
    return traj
end

function LWNonBack(T,D,epsilon)
    t = 0.0 ;  v = 1
    traj = [(0.0,0.0)] ; θprevious = 0.0
    while t < T
        τ = min(T,rand(D)) # duration of the next walk
        θ = θprevious + epsilon - π + rand()*(2π-2epsilon) # direction of the next walk
        next_loc = traj[end] .+ (v*τ) .*(cos(θ),sin(θ)) ; push!(traj,next_loc)
        t += τ
        θprevious = θ
    end
    return traj
end

traj = LWNonBack(1000, Pareto(1.5,2),π*0.)
p4 = plot(traj)
# plot(sd(traj))
plot(p1,p2,p3,p4,layout=(1,4),size=(1600,400))
savefig(pathfig*"sample_traj_levy.pdf")

# Velocity field
function V(x,t,p)
    k,ω,α,v0 = p
    return v0*(1+α*(cos(sum(k.*x) + ω*t)))
end
function traj_velocity_field(p,τ,dt,tmax)
    X = [(0.,0.)]
    η = (0.,0.)
    t = 0.0
    while t < tmax
        η = η .*(1-dt/τ) .+ √(2dt/τ).*Tuple(randn(2))
        x =  X[end]
        nextX = x .+ dt*V(x,t,p) .*η ; push!(X,nextX)
        t += dt
    end
    return X
end

# One trajectory
k = 3*[1,1] ; ω = 0 ; α = 1 ; v0 = 1 ; p = [k,ω,α,v0] # parameters
    τ = 2 # persistent time of the trajectory
    dt = 1E-3 ; tmax = 100 ; times = 0:dt:tmax

traj = traj_velocity_field(p,τ,dt,tmax)
    plot(traj)




function sd(x)
    msd = zeros(length(x))
    for i in eachindex(x)
        msd[i] = sum((x[i][1] - x[1][1],x[i][2] - x[1][2]).^2)
    end
    return msd
end

# Var = 0.1 ; σ
k = [1,1] ; ω = 0 ; α = 0.9 ; v0 = 1 ; p = [k,ω,α,v0] # parameters
τ = 1/2 # persistent time of the trajectory
dt = 1E-3 ; tmax = 100 ; times = 0:dt:tmax
R = 100
MSD = zeros(length(times),R)
for r in 1:R
    X = traj_velocity_field(p,τ,dt,tmax)
    # X = LWNB(1000, Pareto(1.5),π*0.25)
    MSD[:,r] = sd(X)
    # plot(X)
    # plot(MSD)
end
pmsd = plot(xlabel="t",ylabel="MSD",legend=:topleft)
plot!(times[2:end],mean(MSD,dims=2)[2:end,1],axis=:log,label="MSD")
plot!(times[2:100],10times[2:100].^2,c=:black,ls=:dash,label="∼t²")
plot!(times[200:end],times[200:end],c=:black,ls=:dot,label="∼t")

pt1 = plot(traj_velocity_field(p,τ,dt,tmax))
plot(pt,pt1,pmsd,size=(1200,400),layout=(1,3))
savefig(pathfig*"sample_traj_velfield.pdf")

# Soft SAW et MSD
function update_grid(grid,t,memory)
    if memory == "exp"
        halflife = 20
        return grid*(1-1/halflife)
    elseif memory == "powerlaw"
        exponent = 2
        return grid*(1-exponent/(t+10))
    elseif memory == "SAW"
        return grid
    else error("Choose 'memory' among 'SAW', 'exp' and 'powerlaw'. ")
    end
end

function softSAW(L,T,g,memory)
    X = fill((Int(L/2),Int(L/2)),T)
    grid = zeros(L,L) ; grid[Int(L/2),Int(L/2)] = 1
    for t in 2:T
        i,j = X[t-1]
        neighbours = [(i-1,j),(i+1,j),(i,j-1),(i,j+1)] # up, down, left, right
        n_passage = [grid[ind...] for ind in neighbours]
        probas = exp.(-g*n_passage) ; probas = probas / sum(probas)
        nextloc = sample(neighbours, Weights(probas))
        X[t] = nextloc
        grid[nextloc...] += 1
        grid = update_grid(grid,t,memory)
    end
    return X,grid
end
L = 200 ; T = 200 ; memorytype = "SAW" ; gs = [0.1,1,10,100] # [0.1,0.5,1,5,10]
X,grille = softSAW(L,T,1000,memorytype)
plot(X)
heatmap(grille)
for t in 1:2:T
    display(plot(X[1:t],xlims=(70,130),ylims=(70,130)))
    sleep(0.05)
end



# Etudier le MSD
function sd(x)
    msd = zeros(length(x))
    for i in eachindex(x)
        msd[i] = sum((x[i][1] - x[1][1],x[i][2] - x[1][2]).^2)
    end
    return msd
end

L = 300 ; T = 1000 ; memorytype = ["SAW"] ; gs = [0.1,1,10,100]
R = 100 ; MSD = zeros(T,length(gs),length(memorytype),R)
z = @elapsed Threads.@threads for r in 1:R
    for j in eachindex(gs) , k in eachindex(memorytype)
        g = gs[j] ; memory = memorytype[k]
        X,grid = softSAW(L,T,g,memory)
        MSD[:,j,k,r] = sd(X)
    end
end
println("Runtime : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")

styles = [:solid,:dash,:dot]
p = plot(legend=:topleft)
    plot!([NaN,NaN],ls=styles[1],c=:grey,label="SAW")
    plot!([NaN,NaN],ls=styles[2],c=:grey,label="Exp")
    plot!([NaN,NaN],ls=styles[3],c=:grey,label="Power")
    plot!([NaN,NaN],c=1,rib=0,label="g=0.1")
    plot!([NaN,NaN],c=2,rib=0,label="g=1")
    plot!([NaN,NaN],c=3,rib=0,label="g=10")
    plot!([NaN,NaN],c=4,rib=0,label="g=100")
    for k in 1 #eachindex(memorytype)
        for j in eachindex(gs)
            plot!(mean(MSD,dims=4)[2:end,j,k,1],axis=:log,c=j,ls=styles[k])
        end
    end
    plot!(1:L,collect(1:L).^(1.5),c=:black)
    plot!(1:L,collect(1:L).^1.,c=:black)


## Plot trajectories
@unpack L,TV,init,BC,tmax,times,history_locations,R = load("data\\pushMSD.jld")
history_locations_plot = Matrix{Tuple{Float64,Float64}}(undef,length(times),R)
for i in eachindex(history_locations_plot)
    if history_locations[i] === missing
        history_locations_plot[i] = (NaN,NaN)
    else
        history_locations_plot[i] = history_locations[i]
    end
end

p=plot(xlims=(0,L),ylims=(0,L))
    for r in 12
        plot!(history_locations_plot[:,r])
    end
    p

anim = @animate for r in 1:39
    plot(history_locations_plot[:,r],c=:black,size=(512,512),xlims=(0,L),ylims=(0,L),title="Trajectory r = $r")
    end
mp4(anim,"figures\\films\\trajectories.mp4",fps=0.5)

## Track vortices to get MSD
L = 300
    T = 0.1
    Var = 0.1
    init = "isolated"
    BC = "free"
    q = +1 ; r0 = 30
    tmax = 1000 ; transients = 300
    every = 10 ; times = collect(transients:every:tmax)

    history_locations = Vector{Union{Tuple{Int16,Int16},Missing}}(missing,length(times))

thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
    t = 0.0
    while t < transients
        t += dt
        thetas = update_FBC(thetas,omegas,L,T,dt)
    end

last_location,~ = spot_single_default_global(thetas)
    history_locations[1] = last_location
    token = 2 ;
    while t < tmax
        t += dt
        thetas = update_FBC(thetas,omegas,L,T,dt)
        if t > times[token]
            print() ; println("t = $(times[token])")
            try
                lastknown = (Int16(-1),Int16(-1)) # dummy because no -1 defect to be dealt with here
                last_location,alone = spot_single_default_local(thetas,last_location,lastknown,q)
                if alone # there are no other vortices in the surroundings (the '2margin' square)
                    saved_location = last_location
                else # there are other vortices in the surroundings that could pollute the data
                    saved_location = missing
                end
                history_locations[token] = saved_location  # O(1)
            catch e  # we lost our vortex (out of lattice, annihilation)
                println(e)
                printstyled("Warning : Vortex lost, simulation stopped at t = $(round(t,digits=2)). \n"; color = :yellow)
                break # we leave all the remaining values as they are : missing
            end
            token = min(token+1,length(times))
            print() ; println("Defect position = $(history_locations[token-1])")
        end
    end
end

## Analyse MSD
@unpack L,TV,init,BC,tmax,times,history_locations,R = load("data\\pushMSD.jld")
R = size(history_locations,3)
SD  = Matrix{Union{Missing,Float64}}(missing,length(times),R)
# SD  = zeros(length(times),R)
for r in 1:R , i in 1:length(times)
    a,b = history_locations[1,i,r] , history_locations[1,1,r]
    if ismissing(a) || ismissing(b)
        SD[i,r] = missing
    else
        SD[i,r] = sum(abs2,(Float64.(a) .- Float64.(b)))
    end
end
MSD  = zeros(length(times))
MSD_hereatend = zeros(length(times))
MSD_nothereatend = zeros(length(times))
for i in 1:length(times)
    MSD[i] = mean(skipmissing(SD[i,:]))
    MSD_hereatend[i] = mean(skipmissing(SD[i,hereatend]))
    MSD_nothereatend[i] = mean(skipmissing(SD[i,leftlattice]))
end

p=plot()
    plot!(times[2:end],MSD[2:end],axis=:log)
    plot!(times[2:end],MSD_hereatend[2:end],axis=:log)
    # plot!(times[2:end],MSD_nothereatend[2:end],axis=:log)
    # plot!(times[2:end],times[2:end],axis=:log)
    plot!(times[2:end],0.02*(times[2:end].^1.5 ),c=:black)
    # plot!([(times[findfirst(x->x>0, hereatend)],1),(times[findfirst(x->x>0, hereatend)],1E5)])
    plot!(times[2:end],1 .+ 1E3 .* leaving[2:end],axis=:right)

leaving = [sum(ismissing.(SD[i,:])) for i in eachindex(times)]
leftlattice = ismissing.(SD[end,:])
hereatend = BitVector(1 .- leftlattice)

## Considérations energetiques
plot()
L = 100
Vars = 0:0.01:0.3
T = 0.02
tmax = 500
q,r0 = 1,20
R = 10

E_omegas = zeros(length(Vars),R)
E_thetas = zeros(length(Vars),R)

z = @elapsed Threads.@threads for r in 1:R
    for i in eachindex(Vars)
        Var = Vars[i]
        omegas = randn(L,L)*sqrt(Var)
        thetas,omegas,dt = initialize(L,T,Var,"lowtemp",q,r0)
        t = 0.0
        while t < tmax
            t += dt
            thetas = update(thetas,omegas,L,T,dt)
        end
        E_thetas[i,r] = energy(thetas)[2]
        tmp = 0
        for j in 1:L , i in 1:L
            Delta_omegas = [omegas[mod1(i+1,L),j],omegas[mod1(i-1,L),j],omegas[i,mod1(j-1,L)],omegas[i,mod1(j+1,L)]] .- omegas[i,j]
            for i in eachindex(Delta_omegas)
                if abs(Delta_omegas[i]) > 2
                    Delta_omegas[i] = 2*sign(Delta_omegas[i])
                end
            end
            tmp -= sum(sqrt.(1 .- (Delta_omegas/2).^2))
        end
        E_omegas[i,r] = tmp
    end
end
E_omegas_avg = mean(E_omegas,dims=2)[:,1]/L^2 .+ 4
E_thetas_avg = (mean(E_thetas,dims=2)[:,1] .- mean(E_thetas,dims=2)[1,1]) /L^2
plot()
plot!(Vars,E_omegas_avg)
plot!(Vars,E_thetas_avg)
plot!(Vars,Vars/4,c=:black)

## Visualisation elastic energy
L = 200
Var = 0.1
T = 0.02
tmax = 500
q,r0 = 1,20

omegas = randn(L,L)*sqrt(Var)
thetas,omegas,dt = initialize(L,T,Var,"lowtemp",q,r0)
t = 0.0
while t < tmax
    t += dt
    thetas = update(thetas,omegas,L,T,dt)
end

heatmap(energy(thetas)[1].+4)

## Compare runtime of ODE and of manual integration (vector and matrix forms)
    #= Conclusions for L = 400, tmax = 1000, T = 0 and Var = 0.1:
    * Runtime for DifferentialEquations ~ 3100s (both for dt=0.15 fixed and for dt free)
    * Runtime for good old code under matrix form : 250s
    * Runtime for good old code under vector form : 1800s with neighbours precomputing, 370s without (lol)
    =#
function display_thetas(u,cols=cgrad([:black,:blue,:green,:orange,:red,:black]))
    N = length(u)
    L = Int(sqrt(N))
    heatmap(reshape(mod.(u,2π),L,L),c=cols,size=(550,512))
end

L = 400
    T = 0.0
    Var = 0.05
    omegas = randn(L^2)*sqrt(Var)
    p = [L,T,omegas]
    u0 = rand(L^2)*2pi
    tmax = 1000.0
    dt = 0.15
    tspan = (0.0,tmax)

    # ODE
function f(u::Vector{Float64},p::Vector{Any},t)::Vector{Float64} # p = [L,T,omegas]
    L,T,omegas = p
    interactions = zeros(L^2)
    for n in 1:L^2
        i,j = linear_to_square_index(n,L)
        nup    = square_to_linear_index(mod1(i-1,L),j,L)
        ndown  = square_to_linear_index(mod1(i+1,L),j,L)
        nright = square_to_linear_index(i,mod1(j+1,L),L)
        nleft  = square_to_linear_index(i,mod1(j-1,L),L)
        indices_nn = [nup,ndown,nleft,nright]
        interactions[n] = sum(sin.( u[indices_nn,1] .- u[n,1]))
    end
    return omegas + interactions + sqrt(2T)*randn(L^2)
end
prob = ODEProblem(f,u0,tspan,p)
zODE = @elapsed sol = solve(prob,save_everystep=false,dt=dt)
display_thetas(sol.u[end])

    # Manual
thetas = rand(L,L)*2pi
omegas = randn(L,L)*sqrt(Var)
t = 0.0
zmanuel = @elapsed while t < tmax
    t += dt
    thetas = update(thetas,omegas,L,T,dt)
end
display_thetas(thetas)

    # manual vector
function update_vec(thetas::Vector{Float64},omegas::Vector{Float64},nn::Matrix{Int},L::Int,T::Number,dt::Float64)::Vector{Float64}
    interactions = zeros(length(thetas))
    for n in eachindex(thetas)
        θ = thetas[n]
        interactions[n] = sum(sin.( thetas[indices_nn[n,:]] .- θ))
    end
    return thetas + dt*(omegas + interactions) + sqrt(2T*dt)*randn(length(thetas))
end
thetas = rand(L*L)*2pi
omegas = randn(L*L)*sqrt(Var)
indices_nn = zeros(Int,L^2,4)
for n in 1:L^2
    i,j = linear_to_square_index(n,L)
    nup    = square_to_linear_index(mod1(i-1,L),j,L)
    ndown  = square_to_linear_index(mod1(i+1,L),j,L)
    nright = square_to_linear_index(i,mod1(j+1,L),L)
    nleft  = square_to_linear_index(i,mod1(j-1,L),L)
    indices_nn[n,:] = [nup,ndown,nleft,nright]
end
t = 0.0
zmanuel = @elapsed while t < tmax
    t += dt
    thetas = update_vec(thetas,omegas,indices_nn,L,T,dt)
end
display_thetas(thetas)
zmanuel

## Does a naive implementation give the same results as Velocity Verlet ?
function update_underdamped_naive(thetas::Matrix{Float64},thetas_old::Matrix{Float64},omegas::Matrix{Float64},L::Int,T::Number,dt::Float64)
    thetas_new = zeros(L,L)
    for j in 1:L
        for i in 1:L
            angle_neighbours = get_angle_neighbours(thetas,i,j,L) # O(1)
            sin_angle_neighbours = sin.(angle_neighbours .- thetas[i,j])
            # if 1st derivative is backwards
            thetas_new[i,j] =  dt^2 * ( omegas[i,j] + sum(sin_angle_neighbours) + sqrt(2T/dt)*randn()) + dt* (thetas_old[i,j] - thetas[i,j]) + 2thetas[i,j] - thetas_old[i,j]
            # if 1st derivative is centered
            # thetas_new[i,j] =  ( omegas[i,j] + sum(sin_angle_neighbours) + sqrt(2T/dt)*randn()  + thetas[i,j]*2/dt^2 + (thetas_old[i,j])*(1/2/dt - dt^2) ) / (1/dt^2 + 2/dt)
        end
    end
    return thetas_new,thetas
end

function update_underdamped_VV(Q::Matrix{Float64},P::Matrix{Float64},omegas::Matrix{Float64},L::Int,T::Number,dt::Float64)
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
    P_new = P_new*(1-dt) + sqrt(2T*dt)*randn(L,L)
    return Q_new,P_new
end

L = 200
    TV = [(0.3,0.1),(0.2,0.0)]
    init = "hightemp"
    BC = "periodic"
    q = +1 ; r0 = 30
    tmax = 3000
    times = logspace(0.25,tmax,29)
    R = 10

    ξ = zeros(length(TV),length(times),2,R)
    n = zeros(length(TV),length(times),2,R)
    rs = 1:Int(L/2)


z = @elapsed Threads.@threads for r in 1:R
    for i in eachindex(TV)
        T,Var = TV[i]
        thetas_old,omegas,dt = initialize(L,T,Var,init,q,r0) ; thetas = thetas_old + 0.01*omegas
        Q = copy(thetas) ; P = omegas
        t = 0.0 ; token = 1
        # dt = 0.1

        while t < tmax
            t += dt
            ## Naive
            # thetas = update(thetas,omegas,L,T,dt)
            thetas,thetas_old = update_underdamped_naive(thetas,thetas_old,omegas,L,T,dt)
            ## VV
            Q,P = update_underdamped_VV(Q,P,omegas,L,T,dt)

            if t > times[token]
                println("t=$t")
                ξ[i,token,1,1] = corr_length(rs,get_C(thetas))
                ξ[i,token,2,1] = corr_length(rs,get_C(Q))
                n[i,token,1,1] = length(spot_defects(thetas,BC))
                n[i,token,2,1] = length(spot_defects(Q,BC))
                token = min(token+1,length(times))
            end
        end # while
    end # TV
end # r


ξ_avg = mean(ξ,dims=4)
    n_avg = mean(n,dims=4)
    plot(times,n_avg[1,:,1,1],axis=:log)
    plot!(times,n_avg[1,:,2,1])
    plot!(times,n_avg[2,:,1,1],l=:dot)
    plot!(times,n_avg[2,:,2,1],l=:dot)

using BenchmarkTools
thetas_old,omegas,dt = initialize(200,T,Var,init,q,r0) ; thetas = thetas_old + 0.01*omegas
Q = copy(thetas) ; P = omegas

@btime update_underdamped_naive(thetas,thetas_old,omegas,L,T,dt)
@btime update_underdamped_VV(Q,P,omegas,L,T,dt)

## Trouver cause incohérence de la limite overdamped
#= réponse :
    * incohérence pour XY à haute température (0.5 par exemple)
    expliquée par un manque de stat, si on augmente R, on a des
    résultats plus cohérents
    * incohérence pour KM expliquée par une erreur dans le code de
    Velocity Verlet, ce que je prenais en compte c'était pas omega
    comme pour le code overdamped mais 2omega donc xi -> xi/2. =#
L = 100
    global T = 0.2
    Var = 0.2
    ms = [0.0,0.01,0.1,0.5,1.0,2.0]
    init = "hightemp"
    BC = "periodic"
    q = +1 ; r0 = 30
    tmax = 300
    times = logspace(1,tmax,20)
    R = 10

    C = zeros(Int(L/2),length(ms),length(times),R)
    ξ = zeros(length(ms),length(times),R)
    n = zeros(length(ms),length(times),R)

    rs = 1:Int(L/2)

z = @elapsed Threads.@threads for r in 1:R
    for j in eachindex(ms)
        m = ms[j]
        if m == 0
            thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
            t = 0.0 ; token = 1
            while t < tmax
                    t += dt
                    thetas = update(thetas,omegas,L,T,dt)
                if t > times[token]
                    if r==1 println("m = $m , t = $(round(times[token],digits=2))") end
                    tmp = get_C(thetas)
                    C[:,j,token,r] = tmp
                    ξ[j,token,r] = corr_length(rs,tmp)
                    n[j,token,r] = length(spot_defects(thetas,BC))
                    token = min(token+1,length(times))
                end
            end # while
        else # m ≠ 0
            TT = T*sqrt(m)
            Q,omegas,dt = initialize(L,TT,Var,init,q,r0,m)
            t = 0.0 ; token = 1 ; P = omegas
            while t < tmax
                t += dt
                Q,P = update_underdamped_leapfrog(Q,P,omegas,L,TT,dt,m)
                if t > times[token]
                    if r==1 println("m = $m , t = $(round(times[token],digits=2))") end
                    tmp = get_C(Q)
                    C[:,j,token,r] = tmp
                    ξ[j,token,r] = corr_length(rs,tmp)
                    n[j,token,r] = length(spot_defects(Q,BC))
                    token = min(token+1,length(times))
                end
            end
        end
    end # for masses
end # for r
println("Runtime : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes = $(round(z/3660,digits=2)) hours.")

ξ_avg = mean(ξ,dims=3)[:,:,1]
    n_avg = mean(n,dims=3)[:,:,1]
    ξ_std = std(ξ,dims=3)[:,:,1]
pxixy = plot(xlabel="t",ylabel="ξ(t)",title="T=$T,σ²=$Var",legend=:outerright)
    plot!(times,ξ_avg[1,:],c=1,uaxis=:log,label="m=0",rib=0)
    for i in eachindex(ms)
        plot!(times*sqrt(ms[i]),ξ_avg[i,:],c=i,uaxis=:log,label="m=$(ms[i])",rib=ξ_std[2,:])
    end
    plot!(times[10:end-5],2.2sqrt.(times[10:end-5] ./ log.(times[10:end-5])),line=:dash,label="t/log t",c=:black)
    pxixy

## Tc(sigma,m)
L = 200
    Ts = collect(0:0.05:0.45)
    Vars = collect(0:0.025:0.25)
    ms = [0.1,0.5,1.0,2.0]
    init = "lowtemp" ; BC = "periodic"
    q = +1 ; r0 = Int(L/2) # dummy
    tmax = 600 ; every = 10 ; times = collect(every:every:tmax)
    R = 10

Var_critique = NaN*ones(length(Ts),length(ms),R)

z = @elapsed Threads.@threads for r in 1:R
    for j in eachindex(ms)
        for i in eachindex(Ts)
            T = Ts[i] ; m = ms[j]
            for Var in reverse(Vars)
                TT = T*sqrt(m)
                Q,omegas,dt = initialize(L,TT,Var,init,q,r0,m)
                t = 0.0 ; token = 1 ; P = omegas
                println()
                println("Simulation started for T = $T , m = $m and σ² = $Var.")
                while t < tmax
                    t += dt
                    Q,P = update_underdamped_leapfrog(Q,P,omegas,L,TT,dt,m)
                    if t > times[token]
                        # if r==1 println("m = $m , t = $(round(times[token],digits=2))") end
                        if length(spot_defects(Q,BC)) > 0
                            println("Abort simu for σ² = $Var at t = $(round(t,digits=2))")
                            break # sort de l'évolution temporelle et recommence avec une autre Var
                        end
                        token = min(token+1,length(times))
                    end
                end
                if length(spot_defects(Q,BC)) == 0
                    println("Critical Var found for T = $T : Var_c = $Var")
                    Var_critique[i,j,r] = Var
                    break
                end # laisse var critique as is, ie. la dernière
            end # for var
        end # for T
    end # for mass
end # for r
println("Runtime : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes = $(round(z/3660,digits=2)) hours.")
Varavg = mean(Var_critique,dims=3)

# JLD2.jldsave("data/Var_critique_2712.jld2";L,Ts,Vars,Var_critique,ms,init,BC,tmax,every,times,R,q,r0,runtime=z)


p = plot()
    for i in eachindex(ms)
        plot!(Varavg[:,i,1],Ts,m=:circle,c=i+1)
    end
    p
# @unpack Ts,Vars,nb_vortices = load("data/phase_space_number_vortices_hightemp.jld")
# nb_vortices_complement1 = load("data/phase_space_number_vortices_complement_hightemp.jld","nb_vortices")
# nb_vortices_complement2 = load("data/phase_space_number_vortices_complement2_hightemp.jld","nb_vortices")
# nb_final_time = nb_vortices[:,:,end,:]
# nb_final_time[1:5,:,:] = nb_vortices_complement1[:,:,end,:]
# nb_final_time[6:25,:,:] = nb_vortices_complement2[:,:,end,:]
# nb_vortices_avg = mean(nb_final_time,dims=3)[:,:,1]
# Var_critique_overdamped = zeros(length(Vars))
# for i in eachindex(Vars)
#     tmp = findfirst(x->x<2,nb_vortices_avg[:,i])
#     if tmp === nothing
#         Var_critique_overdamped[i] = NaN
#     else
#         Var_critique_overdamped[i] = tmp
#     end
# end
#
# findfirst(x->x==3,0:10)

## Symplectic integration of the underdamped XY hamiltonian dynamics
function update_underdamped_symplectic_alpha1(Q::Matrix{Float64},P::Matrix{Float64},omegas::Matrix{Float64},L::Int,T::Number,dt::Float64,m)
    Q_new = zeros(L,L)
    P_new = zeros(L,L)
    for j in 1:L
        for i in 1:L
            angle_neighbours = get_angle_neighbours(Q,i,j,L) # O(1)
            sin_angle_neighbours = sin.(angle_neighbours .- Q[i,j])
            Ptmp = P[i,j] + dt*sum(sin_angle_neighbours)
            P_new[i,j] = Ptmp + sqrt(2T*dt)*randn()
            Q_new[i,j] = Q[i,j] + dt * Ptmp
        end
    end
    return Q_new,P_new
end

function update_underdamped_symplectic_alpha0(Q::Matrix{Float64},P::Matrix{Float64},omegas::Matrix{Float64},L::Int,T::Number,dt::Float64,m)
    Q_new = zeros(L,L)
    for j in 1:L , i in 1:L
        Q_new[i,j] = Q[i,j] + dt * P[i,j]
    end
    P_new = zeros(L,L)
    for j in 1:L , i in 1:L
        angle_neighbours = get_angle_neighbours(Q_new,i,j,L) # O(1)
        sin_angle_neighbours = sin.(angle_neighbours .- Q_new[i,j])
        P_new[i,j] = P[i,j] + dt*sum(sin_angle_neighbours) + sqrt(2T*dt)*randn()
    end
    return Q_new,P_new
end

function update_underdamped_leapfrog(Q::Matrix{Float64},P::Matrix{Float64},omegas::Matrix{Float64},L::Int,T::Number,dt::Float64,m)
    #= Demarche :
        * the equation is q' = p/m && p' + p/m = f + noise , with p = mv and q = theta
        * first, one integrates the Langevin equation with \gamma = 0, ie p'(t) = f(t)
        * then, one introduces the damping and the noise p(t+dt) = p(t) (1-dt/m) + sqrt(2T*dt/m)*randn()
    =#
    Q_new = zeros(L,L)
    P_new = zeros(L,L)
    for j in 1:L , i in 1:L
        angle_neighbours_old = get_angle_neighbours(Q,i,j,L) # O(1)
        sin_angle_neighbours_old = sin.(angle_neighbours_old .- Q[i,j])
        Q_new[i,j] = Q[i,j] + dt * P[i,j] +  dt^2 / 2 * (sum(sin_angle_neighbours_old))
        # Q_new[i,j] = Q[i,j] + dt * P[i,j] +  dt^2 /2/m * (sum(sin_angle_neighbours_old) - P[i,j] )
    end
    for j in 1:L , i in 1:L
        P_new[i,j] = P[i,j] + dt/2 * (sum(sin.(get_angle_neighbours(Q,i,j,L) .- Q[i,j])) + sum(sin.(get_angle_neighbours(Q_new,i,j,L) .- Q_new[i,j])) )
        # P_new[i,j] = ( P[i,j] + dt/2m * ( sum(sin.(get_angle_neighbours(Q,i,j,L) .- Q[i,j])) + sum(sin.(get_angle_neighbours(Q_new,i,j,L) .- Q_new[i,j])) - P[i,j]) ) / (1+dt/2m)
    end
    P_new = P_new*(1-dt/sqrt(m)) + sqrt(2T*dt/sqrt(m))*randn(L,L)
    return Q_new,P_new
end

L = 200
    T = 0.3
    Var = 0.05
    m = 10.
    init = "hightemp"
    BC = "periodic"
    q = +1 ; r0 = Int(L/2) # dummy
    transients = 0 ; tmax = 1000 ; times = collect(transients:1:tmax)

t = 0.0 ; token = 2 ; Q,omegas,dt = initialize(L,T,Var,init,q,r0) ; P = omegas ; dt = sqrt(m) / 10
    Q_hist = zeros(L,L,length(times))
    Q_hist[:,:,1] = Q ;

z = @elapsed while t < tmax
    t = t + dt ; Q,P = update_underdamped_leapfrog(Q,P,omegas,L,T,dt,m)
    # t = t + dt ; Q = update(Q,omegas,L,T,dt)
    if t > times[token]
        Q_hist[:,:,token] = Q
        token = min(token+1,length(times))
    end
end
println("Simulation Time : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")
plot(exp_moving_average(Q_hist[rand(1:L),rand(1:L),:],1))

# Movies
zm = @elapsed anim = @animate for i in 1:size(Q_hist,3)
    p = display_thetas(Q_hist[:,:,i],title="t=$(round(times[i],digits=2))",defects=false)
end
println("Time Movie : $(round(Int,zm)) seconds = $(round(zm/60,digits=2)) minutes.")
mp4(anim,"figures/films/underdamped_symplectic_m$(m)_Var$(Var).mp4",fps=10);
# mp4(anim,"figures/films/control_T$T.mp4",fps=10);


## Correlation |∇θ| et creation de vortex
L = 200
    T = 0.05
    Var = 0.25
    init = "lowtemp"
    BC = "periodic"
    q = +1 ; r0 = 100 # dummy
    transients = 300 ; tmax = 600 ; times = collect(transients:.5:tmax)
    R = 10
    c = 1 # leads to (2c+1)x(2c+1) squares around locations of creation to compute |Δθ|


    function nabla_theta(thetas)
        l,L = size(thetas)
        result = Float64[]
        for i in 2:l-1 , j in 2:L-1
            neighbours = get_angle_neighbours(thetas,i,j,L)
            push!(result,(neighbours .- thetas[i,j])...)
        end
        return mod.(result.-pi,2pi) .- pi
    end

Dtheta_creation = Float64[]
save_one_config = zeros(L,L)
Threads.@threads for r in 1:R
    t = 0.0 ; token = 1 ; thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
    if r == 1 print() ; println("Simulation started.") end
    while t < transients t = t + dt ; thetas = update(thetas,omegas,L,T,dt) end
    fiche_vortices,fiche_antivortices,vortices_old,antivortices_old = initialize_fiches_defects(thetas,BC,t)
    thetas_old = thetas
    if r == 1 print() ; println("Transients over.") end
    while t < tmax
        t = t + dt ; thetas = update(thetas,omegas,L,T,dt)
        if round(Int,t,RoundDown) ≥ times[token]
            if r == 1 print() ; println("t = $(times[token])/$tmax") end
            vortices_new,antivortices_new = spot_defects_separated(thetas,BC)
            m_old = length(fiche_vortices)
            fiche_vortices,fiche_antivortices = defect_tracker(thetas,BC,fiche_vortices,fiche_antivortices,vortices_old,antivortices_old,t)
            m_new = length(fiche_vortices)
            if m_new > m_old
                loc_new_creations = [fiche_vortices[end-i][4] for i in 0:m_new-m_old-1]
                for i in 1:length(loc_new_creations)
                    i_creation,j_creation = loc_new_creations[i]
                    push!(Dtheta_creation, nabla_theta(thetas_old[max(i_creation-c,1):min(i_creation+c,L),max(j_creation-c,1):min(j_creation+c,L)])...)
                end
            end
            vortices_old,antivortices_old = vortices_new,antivortices_new
            thetas_old = thetas
            token = min(length(times),token+1)
        end
        if r == 1 save_one_config[:,:] = thetas end
    end
end

histogram(nabla_theta(save_one_config),normalize=true,yaxis=:log)
    histogram!(Dtheta_creation,normalize=true,yaxis=:log)

length(Dtheta_creation)
histogram(Dtheta_creation,normalize=true,yaxis=:log)

# thetasKM = thetas
# thetasXY = thetas


# histogram(nabla_theta(save_one_config),normalize=true,yaxis=:log)
histogram(Dtheta_creation,normalize=true,yaxis=:log)
    xx = -2.5:0.01:2.5
    Var = 0.45
    plot!(xx,exp.(-xx.^2/(2*Var))/(2π*Var)^0.5,c=:red,lw=2)

## Reproducing superdiffusion
T_rotation = 5 ; T_translation = 0.0 ;
    dt = 0.01 ; tmax = 3
    x = [0.5,0.5] ; theta = 0 ; t = 0.0 ;
    x_history = zeros(1+round(Int,tmax/dt),2) ; theta_history = zeros(1+round(Int,tmax/dt))
    x_history[1,:] = x ; theta_history[1] = theta

    function U(x,theta,k=[1,1])
        return (1 .+cos.(k .* x)) .* [cos(theta),sin(theta)]
    end

    for i in 1:1+round(Int,tmax/dt)
        t += dt
        theta += sqrt(2T_rotation*dt)*randn()
        x = mod.(x .+ dt*U(x,theta) + sqrt(2T_translation*dt)*randn(2),1)
        x_history[i,:] = x ; theta_history[i] = theta
    end

    plot(xlims=(0,1),ylims=(0,1))
    scatter!([(a[1],a[2]) for a in eachrow(x_history)],ms=1.5)

anim = @animate for i in 1:2:1+round(Int,tmax/dt)
    plot(xlims=(0,1),ylims=(0,1),size=(512,512))
    scatter!([(x_history[j,1],x_history[j,2]) for j in 1:i],ms=1.5)
end
mp4(anim,"figures/films/mini_model_test.mp4",fps=10);

## Number of defects over time (9h max)
L = 50
    TV = [(0.2,0.03),(0.5,0.03),(0.2,0.15),(0.5,0.15)]
    ms = [0.01,0.1,0.5,1.0,5.0]
    init = "hightemp"
    BC = "periodic"
    q = +1 ; r0 = 30
    tmax = 500
    times = logspace(1,tmax,10)
    R = 10

    n  = zeros(length(TV),length(ms),length(times),R)

# z = @elapsed Threads.@threads for r in 1:R
#     for i in eachindex(TV)
#         T,Var = TV[i]
#         for j in eachindex(ms)
#             m = ms[j]
#             print() ; println("Simulation started for r = $r , T = $T , σ² = $Var and m = $m ") ;
#             thetas_new,omegas,dt = initialize(L,T,Var,init,q,r0)
#             t = 0.0 ; token = 1 ; thetas_old = copy(thetas_new)
#             while t < tmax
#                     t += dt
#                     thetas_tmp = update_underdamped(thetas_new,thetas_old,omegas,L,T,dt,m)
#                     thetas_old = thetas_new
#                     thetas_new = thetas_tmp
#                 if t > times[token]
#                     n[i,j,token,r] = length(spot_defects(thetas_new,BC))
#                     token = min(token+1,length(times))
#                 end
#             end # while
#         end # for j (masses)
#     end # for i (TV)
# end # for r (realisations)
z = @elapsed Threads.@threads for r in 1:R
    r = 1
    for i in eachindex(TV)
        T,Var = TV[i]
        for j in eachindex(ms)
            m = ms[j]
            print() ; println("Simulation started for r = $r , T = $T , σ² = $Var and m = $m ") ;
            Q,omegas,dt = initialize(L,T,Var,init,q,r0,m)
            P = omegas
            t = 0.0 ; token = 1
            while t < tmax
                    t += dt
                    Q,P = update_underdamped_leapfrog(Q,P,omegas,L,T,dt,m)
                if t > times[token]
                    n[i,j,token,r] = length(spot_defects(Q,BC))
                    token = min(token+1,length(times))
                end
            end # while
        end # for j (masses)
    end # for i (TV)
end # for r (realisations)
println("Simulation Finished : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes = $(round(z/3600,digits=2)) hours.")

n_avg = mean(n,dims=4)[:,:,:,1]
# JLD2.jldsave("data/underdamped_symplectic_nt_1612.jld2";L,TV,ms,init,tmax,times,R,BC,n,n_avg)

@unpack L,TV,ms,init,tmax,times,R,BC,n,n_avg = JLD2.load("data/underdamped_symplectic_nt_1612.jld2")

j=1
    p=plot(title="T=$(TV[j][1]) , σ²=$(TV[j][2])")
    for i in reverse(eachindex(ms))
        plot!(1/sqrt(ms[i])*times,log.(1 .+ n_avg[j,i,:]),c=i,line=:solid,xaxis=:log)
        # plot!(times,n_avg[2,i,:],c=i,line=:dot,axis=:log)
        # plot!(times,n_avg[3,i,:],c=i,line=:dash,axis=:log)
    end
    plot!(times[15:end],1E3 * log.(times[15:end])./ times[15:end],c=:black)
@unpack L,TV,ms,init,tmax,times,R,BC,n,n_avg = JLD2.load("data/underdamped_nt_1012.jld2")
    for i in reverse(eachindex(ms))
        plot!(times,log.(1 .+ n_avg[2,i,:]),c=i,line=:solid,xaxis=:log)
    end
    p
n_avg[2,1,:]

## Movies underdamped dynamics with symplectic leapfrog integrator
L = 200
    T = 0.1
    Var = 0.1
    m = .01
    init = "pair"
    BC = "periodic"
    q = +1 ; r0 = 30
    tmax = 500
    every = 1 ; times = collect(0:every:tmax)

    Q,omegas,dt = initialize(L,T,Var,init,q,r0,m)
    dt = determine_dt(T,Var,m,L)
    t = 0.0 ; token = 2 ; P = omegas
    thetas_hist = zeros(L,L,length(times))
    thetas_hist[:,:,1] = Q

z = @elapsed while t < tmax
        t += dt
        Q,P = update_underdamped_leapfrog(Q,P,omegas,L,T,dt,m)
    if t > times[token]
        thetas_hist[:,:,token] = Q
        token = min(token+1,length(times))
    end
end
println("Simulation Time : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")
# plot(exp_moving_average(thetas_hist[12,12,:],20))

# Movies
zm = @elapsed anim = @animate for i in 1:size(thetas_hist,3)
    p = display_thetas(thetas_hist[:,:,i],title="t=$(round(times[i],digits=2))",defects=false)
end
println("Time Movie : $(round(Int,zm)) seconds = $(round(zm/60,digits=2)) minutes.")
mp4(anim,"figures/films/underdamped_leapfrog_m$(m)_Var$(Var)_T$(T).mp4",fps=25);

## Superdiffusion 3/2 Generalized Langevin Equation with memory kernel
using SpecialFunctions
t_memory = 2
    dt = 0.01 ; tt = collect(0:dt:1)
    α = 12
    dw = 0.01 ;  ; ww = Complex(α) .+ im*collect(-1:dw:1)

    K(t) = 1/(.1+t)
    K(t) = exp(-t)
    K(t) = exp(-t^2/2*0.0001)/sqrt(2pi*0.0001)

    K_t = K.(tt)
    K_w = zeros(Complex,length(ww))
    for i in eachindex(ww)
        K_w[i] = dt*sum([exp.(-ww[i]*t)*K(t) for t in tt])
    end
    chi_x = zeros(Complex,length(tt))
    for i in eachindex(tt)
        chi_x[i] = 1/(2π*im) * dw * sum([exp(ww[j]*tt[i])/(ww[j]*K_w[j]) for j in eachindex(ww)])
    end
    chi = zeros(Complex,length(tt))
    for i in eachindex(tt)
        chi[i] = sum(chi_x[1:i])*dt
    end
    plot(real(chi_x))
    plot(imag(chi))


## Track vortices to get MSD
# L = 300
#     T = 0.2
#     Var = 0.03
#     ms = [0.0,0.1,0.5,1.0,5.0]
#     ms = [0.0,0.1]
#     init = "isolated"
#     BC = "free"
#     q = +1 ; r0 = 30
#     tmax = 1000 ; transients = 300
#     every = 10 ; times = collect(transients:every:tmax)
#
#     history_locations = Matrix{Union{Tuple{Int16,Int16},Missing}}(missing,length(times),length(ms))
#     # history_locations = Vector{Union{Tuple{Int16,Int16},Missing}}(missing,length(times))
#
# for k in eachindex(ms)
#     m = ms[k]
#     print() ; println("Simulation started for m = $m")
#
#     thetas_new,omegas,dt = initialize(L,T,Var,init,q,r0)
#     t = 0.0 ; thetas_old = copy(thetas_new)
#     while t < transients
#         t += dt
#         thetas_tmp = update_underdamped_FBC(thetas_new,thetas_old,omegas,L,T,dt,m)
#         thetas_old = thetas_new
#         thetas_new = thetas_tmp
#     end
#
#     last_location,~ = spot_single_default_global(thetas_new)
#     history_[1,k] = last_location
#     token = 2 ;
#     while t < tmax
#             t += dt
#             thetas_tmp = update_underdamped_FBC(thetas_new,thetas_old,omegas,L,T,dt,m)
#             thetas_old = thetas_new
#             thetas_new = thetas_tmp
#         if t > times[token]
#             print() ; println("t = $(times[token])")
#             try
#                 lastknown = (Int16(-1),Int16(-1)) # dummy because no -1 defect to be dealt with here
#                 last_location,alone = spot_single_default_local(thetas_new,last_location,lastknown,q)
#                 if alone # there are no other vortices in the surroundings (the '2margin' square)
#                     saved_location = last_location
#                 else # there are other vortices in the surroundings that could pollute the data
#                     saved_location = missing
#                 end
#                 history_[token,k] = saved_location  # O(1)
#             catch e  # we lost our vortex (out of lattice, annihilation)
#                 println(e)
#                 printstyled("Warning : Vortex lost, simulation stopped at t = $(round(t,digits=2)). \n"; color = :yellow)
#                 break # we leave all the remaining values as they are : missing
#             end
#             token = min(token+1,length(times))
#             print() ; println("Defect position = $(history_[token-1,k])")
#         end
#     end
# end

## Data generation for Coarsening dynamics
L = 200
    TV = [(0.2,0.0),(0.2,0.1)]
    ms = [0.0,0.01,0.1,0.5,1.0,5.0]
    init = "hightemp"
    BC = "periodic"
    q = +1 ; r0 = 30
    tmax = 3000
    times = logspace(1,tmax,30)
    R = 10

    rs = collect(1:Int(L/2))
    C  = zeros(length(TV),length(ms),Int(L/2),length(times),R)
    ξ  = zeros(length(TV),length(ms),length(times),R)

z = @elapsed Threads.@threads for r in 1:R
    for i in eachindex(TV)
        T,Var = TV[i]
        for j in eachindex(ms)
            m = ms[j]
            print() ; println("Simulation started for r = $r , T = $T , σ² = $Var and m = $m ") ;
            thetas_new,omegas,dt = initialize(L,T,Var,init,q,r0)
            t = 0.0 ; token = 1 ; thetas_old = copy(thetas_new)
            while t < tmax
                    t += dt
                    thetas_tmp = update_underdamped(thetas_new,thetas_old,omegas,L,T,dt,m)
                    thetas_old = thetas_new
                    thetas_new = thetas_tmp
                if t > times[token]
                    C_tmp = get_C(thetas_new)
                    C[i,j,:,token,r] = C_tmp
                    ξ[i,j,token,r]   = corr_length(rs,C_tmp)
                    token = min(token+1,length(times))
                end
            end # while
        end # for j (masses)
    end # for i (TV)
end # for r (realisations)
println("Simulation Finished : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")
ξ_avg = mean(ξ,dims=4)[:,:,:,1]
# JLD2.jldsave("data/underdamped_coarsening_0812.jld2";L,TV,ms,init,tmax,times,R,BC,C,rs,ξ_avg,ξ)

p=plot()
    for i in reverse(eachindex(ms))
        plot!(times,ξ_avg[1,i,:],c=i,line=:solid,axis=:log)
        plot!(times,ξ_avg[2,i,:],c=i,line=:dot,axis=:log)
    end
    p

## First movies underdamped dynamics
L = 200
    T = 0.2
    Var = 0.2
    m = 1.0
    init = "pair"
    BC = "periodic"
    q = +1 ; r0 = 30
    tmax = 1000
    every = 10 ; times = collect(0:every:tmax)

    thetas_new,omegas,dt = initialize(L,T,Var,init,q,r0)
    # dt = 0.05
    t = 0.0 ; token = 2 ; thetas_old = copy(thetas)
    thetas_hist = zeros(L,L,length(times))
    thetas_hist[:,:,1] = thetas_new

z = @elapsed while t < tmax
        t += dt
        thetas_tmp = update_underdamped(thetas_new,thetas_old,omegas,L,T,dt,m)
        thetas_old = thetas_new
        thetas_new = thetas_tmp
    if t > times[token]
        thetas_hist[:,:,token] = thetas_new
        token = min(token+1,length(times))
    end
end
println("Simulation Time : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")
# plot(exp_moving_average(thetas_hist[12,12,:],20))

# Movies
zm = @elapsed anim = @animate for i in 1:size(thetas_hist,3)
    p = display_thetas(thetas_hist[:,:,i],title="t=$(round(times[i],digits=2))",defects=false)
end
println("Time Movie : $(round(Int,zm)) seconds = $(round(zm/60,digits=2)) minutes.")
mp4(anim,"figures/films/underdamped_m$(m)_test.mp4",fps=10);



## C(r,t) for Var = 0 and Var = 0.05
L = 200
    Ts = [0.4]
    Vars = [0.1]
    init = "hightemp"
    BC = "periodic"
    q = +1 ; r0 = 100 # dummy
    tmax = 10000
    times = logspace(1,tmax,25)
    R = 10

C = zeros(length(Ts),length(Vars),Int(L/2),length(times),R)

z = @elapsed for i in eachindex(Ts) , j in eachindex(Vars)
    T = Ts[i]
    Var = Vars[j]
    Threads.@threads for r in 1:R
        t = 0.0 ; token = 1 ; thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
        while t < tmax
            t = t + dt
            thetas = update(thetas,omegas,L,T,dt)
            if round(Int,t,RoundDown) ≥ times[token]
                if r == 1 print() ; println("t=$(times[token]) / $tmax") end
                C[i,j,:,token,r] = get_C(thetas)
                token = min(length(times),token+1)
            end
        end
    end
end
println("Time = $(round(z/60,digits=1)) minutes.")

# JLD2.jldsave("data/Crt_paper_Var.jld2";L,Ts,Vars,tmax,init,C,BC,times,R)

## MSD of vortices
L = 100
    T = 0.2
    Vars = [0.0]
    init = "lowtemp"
    BC = "periodic"
    q = +1 ; r0 = 100 # dummy
    transients = 1000 ; tmax = 3000
    times = transients:5:tmax
    save_thetas = zeros(L,L,length(times))

## Old routine
L = 100
    T = 0.2
    Var = 0.0
    r0 = (round(Int16,L/2),round(Int16,L/2)) # initial position
    tmax = 50
    t0 = 1
    nsave = 50
    q = +1
    R = 1
    ts = Array(2:2:tmax)
    # ts = round.(Int,range(t0,tmax,length=nsave))

    Δx = NaN*zeros(nsave,R)
    Δy = NaN*zeros(nsave,R)

z=@elapsed Threads.@threads for real in 1:R
    t = 0.0 ; thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
    token_time_t = 1 ; last_loc = r0
    while t < tmax
        t = t + dt
        thetas = update(thetas,omegas,L,T,dt)

        try
            if round(Int,t,RoundDown) ≥ ts[token_time_t]
                last_loc = spot_single_default_local(thetas,last_loc,(Int16.(-1),Int16.(-1)),Q)[1]
                Δx[token_time_t,real],Δy[token_time_t,real] = last_loc .- r0
                token_time_t             = min(token_time_t + 1,nsave)  # not very pretty but here to handle border effects due to rounding issues
            end
        catch e
            if t<tmax println("$real , t = $t : ",e) end
            break
        end
    end
end
println("Runtime = $(round(Int,z)) seconds")


## Snaps to show the progression of the defects on the boundaries
L = 200
    T = 0.1
    Var = 0.15
    init = "lowtemp"
    BC = "periodic"
    q = +1 ; r0 = 100 # dummy
    transients = 1000 ; tmax = 3000
    times = transients:5:tmax
    save_thetas = zeros(L,L,length(times))

t = 0.0 ; thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
    p = display_thetas(thetas) ; highlight_defects(p,thetas)
ztr = @elapsed while t < transients
        t += dt ; update(thetas,omegas,L,T,dt)
    end
    println("Time = $(round(ztr/60,digits=1)) minutes.")
p = display_thetas(thetas) ; highlight_defects(p,thetas)

token = 1
z = @elapsed while t < tmax
    t += dt ; update(thetas,omegas,L,T,dt)
    if t > times[token]
        save_thetas[:,:,token] = thetas
        token = min(length(times),1+token)
    end
end
println("Time = $(round(z/60,digits=1)) minutes.")
p = display_thetas(thetas) ; highlight_defects(p,thetas)

# anim = @animate for i in 1:length(times)
#     p=display_thetas(save_thetas[:,:,i]) ; ; highlight_defects(p,save_thetas[:,:,i])
# end
# mp4(anim,"figures/films/film_for_snaps_displacement_defect_T$(T)_Var$(Var).mp4",fps=10);
plots[6] = heatmap(mod.(save_thetas[30:80,80:200,times_to_plot[6]]',2pi),c=cols,ticks=false,colorbar=false,size=(140,240))
for defect in spot_defects(save_thetas[30:80,80:200,times_to_plot[6]],"free")
    if defect[3] > 0  scatter!((defect[1:2]), m = (8, 12.0, :circle,:transparent, stroke(2, :white))) end
end
plots[6]
times_to_plot = [260,270,280,290,300,308]
    plots = Vector{Any}(undef,length(times_to_plot))
    for i in eachindex(times_to_plot)
        plots[i] = heatmap(mod.(save_thetas[30:80,80:200,times_to_plot[i]]',2pi),c=cols,ticks=false,colorbar=false,size=(140,240))
        for defect in spot_defects(save_thetas[30:80,80:200,times_to_plot[i]],"free")
            if defect[3] > 0  scatter!((defect[1:2]), m = (8, 12.0, :circle,:transparent, stroke(2, :white))) end
        end
        ylims!(0,120)
        xlims!(0,50)
    end
    p=plot(plots...,layout=(1,length(plots)),size=(length(plots)*140,240))
    p
savefig("figures\\figures_paper\\defects_surf.pdf")

# JLD2.jldsave("data\\defects_surf.jld2";L,T,Var,init,BC,transients,tmax,times,save_thetas)
@unpack L,T,Var,init,BC,transients,tmax,times,save_thetas = load("data\\defects_surf.jld2")


## Lieux de création et annihilation comparaison XY et NKM
L = 200
    T = 0.4
    Var = 0.
    init = "lowtemp"
    BC = "periodic"
    q = +1 ; r0 = 50 # dummy
    nb_creation_wanted = 100
    tmax = 1000

    #= Explanation of the data structure
        - each fiche_vortices[i] is the data corresponding to a vortex
        - this data is a vector composed of
            1. the vortex ID ,
            2. the last known location (i,j)
            3. the creation time ,
            4. the location of the creation
            5. the annihilation time (nothing until known) ,
            6. the location of the annihilation (nothing until known) ---> to code !!
            7. the ID of the antivortex with which it annihilated (nothing until known)
    =#

t = 0.0 ; thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
fiche_vortices,fiche_antivortices,vortices_old,antivortices_old = initialize_fiches_defects(thetas,BC,t)
while length(fiche_vortices) < nb_creation_wanted || t > tmax
    println("Already $(length(fiche_vortices)) vortex creations.")
    t += dt ; update(thetas,omegas,L,T,dt)
    vortices_new,antivortices_new = spot_defects_separated(thetas,BC)
    fiche_vortices,fiche_antivortices = defect_tracker(thetas,BC,fiche_vortices,fiche_antivortices,vortices_old,antivortices_old,t)
    vortices_old,antivortices_old = vortices_new,antivortices_new
end

# p = display_thetas(thetas) ; highlight_defects(p,thetas)
# for i in 1:100 t += dt ; update(thetas,omegas,L,T,dt) end
#     vortices_new,antivortices_new = spot_defects_separated(thetas,BC)
#     fiche_vortices,fiche_antivortices = defect_tracker(thetas,BC,fiche_vortices,fiche_antivortices,vortices_old,antivortices_old,t)
#     vortices_old,antivortices_old = vortices_new,antivortices_new
#     p = display_thetas(thetas) ; highlight_defects(p,thetas)

# fiche_vortices
# p = display_thetas(thetas) ; highlight_defects(p,thetas)

function skinnycross(w = 0.1)
    x = [w,w,1]; rx = reverse(x)
    y = [1,w,w]; ry = reverse(y)
    Shape(vcat(x,rx,-x,-rx), vcat(y,-ry,-y,ry))
end
skinnyx(w=0.1) = rotate(skinnycross(w), 0.25π)

plot(xlims=(0,L),ylims=(0,L))
    scatter!(location_creation(fiche_vortices),m=:circle,c=:limegreen,ms=6)
    scatter!(location_annihiliation(fiche_vortices),m=skinnyx(0.2),c=:red)
    # lens!([145, 160], [180,188], inset = (1, bbox(0.065, 0.05, 0.2, 0.2)))
# savefig("figures\\figures_paper\\_creations_NKM_T$(T)_Var$(Var).svg")
# savefig("figures\\figures_paper\\_creations_XY_T$(T).svg")


## Deuxieme tentative pour comprendre la normalisation
L = 200
    T = 0.2 ; Var = 0.
    init = "hightemp"
    BC = "periodic"
    q = +1 ; r0 = Int(L/2) # dummy
    tmax = 100 ; times = collect(50:50:tmax) .- 1
    R = 100
    save_thetas = zeros(L,L,length(times),R)
    save_omegas = zeros(L,L,R)

zt = @elapsed Threads.@threads for r in 1:R
    t = 0.0 ; thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
    save_omegas[:,:,r] = omegas
    token = 1
    z1 = @elapsed while t < tmax
        t += dt ; update(thetas,omegas,L,T,dt)
        if t>times[token]
            save_thetas[:,:,token,r] = thetas
            token = min(length(times),1+token)
        end
    end
    # thetas_ext = zeros(2L,2L)
    # thetas_ext[1:L,1:L] = thetas_ext[L+1:end,1:L] = thetas_ext[1:L,L+1:end] = thetas_ext[L+1:end,L+1:end] = thetas
    # omegas_ext = zeros(2L,2L)
    # omegas_ext[1:L,1:L] = omegas_ext[L+1:end,1:L] = omegas_ext[1:L,L+1:end] = omegas_ext[L+1:end,L+1:end] = omegas
    # z2 = @elapsed while t < tmax
    #     t += dt ; update(thetas_ext,omegas_ext,2L,T,dt)
    #     if t>times[token]
    #         save_thetas[:,:,token,r] = thetas_ext
    #         token = min(length(times),1+token)
    #     end
    # end
    println()
    println("Elapsed time = $(round((z1)/60,digits=1)) minutes for r = $r.")
    println("Expected time for enlargment = $(round((0.4z1)/60,digits=1)) minutes for r = $r.")
end

# Fake enlargment
save_thetas_ext_finaltime = zeros(2L,2L,R)
dt = determine_dt(T,Var)
Threads.@threads for r in 1:R
    t = 0.0
    thetas = save_thetas[:,:,end,r]
    omegas = zeros(L,L)
    thetas_ext = zeros(2L,2L)
    thetas_ext[1:L,1:L] = thetas_ext[L+1:end,1:L] = thetas_ext[1:L,L+1:end] = thetas_ext[L+1:end,L+1:end] = thetas
    omegas_ext = zeros(2L,2L)
    omegas_ext[1:L,1:L] = omegas_ext[L+1:end,1:L] = omegas_ext[1:L,L+1:end] = omegas_ext[L+1:end,L+1:end] = omegas
    while t < 0.2tmax
        t += dt ; update(thetas_ext,omegas_ext,2L,T,dt)
    end
    save_thetas_ext_finaltime[:,:,r] = thetas_ext
end

# JLD2.jldsave("data\\data_gr_0311.jld2";L,init,BC,tmax,save_thetas,save_thetas_ext_finaltime,R,times,T,Var,dt)
# @unpack L,init,BC,tmax,save_thetas,save_thetas_ext_finaltime,R,times,T,Var,dt = JLD2.load("data\\data_gr_0411.jld2")

g_same = zeros(Int(L/2),length(times),R)
    g_opposite = zeros(Int(L/2),length(times),R)
    for r in 1:R ,t in eachindex(times)
            g_opposite[:,t,r],g_same[:,t,r] = compute_gr(save_thetas[:,:,t,r])
    end
    g_opposite_avg = mean(g_opposite,dims=3)[:,:,1]
    g_same_avg = mean(g_same,dims=3)[:,:,1]

g_same_ext_finaltime = zeros(L,R)
    g_opposite_ext_finaltime = zeros(L,R)
    for r in 1:R
            g_opposite_ext_finaltime[:,r],g_same_ext_finaltime[:,r] = compute_gr(save_thetas_ext_finaltime[:,:,r])
    end
    g_opposite_ext_finaltime_avg = mean(g_opposite_ext_finaltime,dims=2)[:,1]
    g_same_ext_finaltime_avg = mean(g_same_ext_finaltime,dims=2)[:,1]


p=plot()
    for t in [1,2]
        plot!(g_opposite_avg[3:end,t],label="t=$(1+times[t])",rib=0)
    end
    p

p=plot()
    # plot!(g_same_ext_finaltime_avg[3:end],axis=:log)
    plot!(g_opposite_ext_finaltime_avg[3:end],yaxis=:log)
    xlims!(1,200)
    ylims!(.9,2)
    # plot!(g_same_ext_finaltime_avg)


&
plot(log.(4:200),log.(g_opposite_ext_finaltime_avg[4:end]))

t = 2
    p=heatmap(mod.(save_thetas[:,:,t,1],2pi)',c=cols)
    highlight_defects(p,save_thetas[:,:,t,1])
length(spot_defects_separated(save_thetas[:,:,t,1],BC)[1])

t = 5
    plot(g_opposite_avg[:,t],label="t=$(1+times[t])",rib=0)
plot((compute_gr(save_thetas[:,:,t,1])[2]))

## Nouvelle tentative pour g(r)
L = 200
    T = 0.6 ; Var = 0.
    init = "hightemp"
    BC = "periodic"
    q = +1 ; r0 = Int(L/2) # dummy
    tmax = 300 ; times = collect(50:50:tmax) .- 1
    R = 6

    grr = zeros(L,R)

zt = @elapsed Threads.@threads for r in 1:R
    t = 0.0 ; thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
    z = @elapsed while t < tmax
        t += dt ; update(thetas,omegas,L,T,dt)
    end
    println()
    println("Part 1 done for r = $r. ETA in $(round(4*0.3*z/60,digits=1)) minutes")
    thetas_ext = zeros(2L,2L)
    thetas_ext[1:L,1:L] = thetas_ext[L+1:end,1:L] = thetas_ext[1:L,L+1:end] = thetas_ext[L+1:end,L+1:end] = thetas
    omegas_ext = zeros(2L,2L)
    omegas_ext[1:L,1:L] = omegas_ext[L+1:end,1:L] = omegas_ext[1:L,L+1:end] = omegas_ext[L+1:end,L+1:end] = omegas
    while t < 1.3tmax
        t += dt ; update(thetas_ext,omegas_ext,2L,T,dt)
    end


    grr[:,r,token] = compute_gr(thetas_ext)
end

plot(grr[2:end,1],axis=:log,xlabel="r",ylabel="g(r)",)
gg = zeros(L-1)
for i in eachindex(gg)
    if mean(grr,dims=2)[i,1]-1 > 0
        gg[i] = mean(grr,dims=2)[i,1]-1
    else
        gg[i] = NaN
    end
end
plot(gg,axis=:log,xlabel="r",ylabel="g(r)-1",m=true,ms=1.5,line=false)
# savefig("figures\\gr_T0.6.pdf")
plot!(1 ./(log.(1:L)))



## Comprendre la dynamique de la g(r)
L = 200
    T = 0.4 ; Var = 0.
    init = "hightemp"
    BC = "periodic"
    q = +1 ; r0 = Int(L/2) # dummy
    R = 20
    tmax = 200 ; times = collect(50:50:tmax) .- 1
    normalisation = Matrix{Vector{Float64}}(undef,length(times),R)
    distances = Matrix{Vector{Float64}}(undef,length(times),R)

save_thetas = ones(L,L)
z = @elapsed Threads.@threads for r in 1:R
    t = 0.0 ; token_t = 1 ;
    thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
    if r == 1 save_thetas = thetas end

    while t < tmax
        t += dt ; update(thetas,omegas,L,T,dt)
        if round(Int,t) > times[token_t]
            vortices,antivortices = spot_defects_separated(thetas,BC)
            m = length(vortices)
            println() ; println("Number of vortices = $m at t = $(round(Int,t)) for r = $r")

            # Opposite charge
            tmp = []
            for i in 1:m , j in 1:m
                push!(tmp,dist(vortices[i],antivortices[j],L))
            end
            distances[token_t,r] = tmp
            token_t = min(token_t + 1,length(times))
        end
    end
end
println("Runtime : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")
dr = 5 ; rr = 0:dr:Int(L/2)
    aire_system = L^2  # aire d'un cercle et non L² car on a filtered out les distances > L/2 lors de la génération des données
    aire_couronne = π*(2 .* rr[1:end-1] * dr .+ dr^2) # on garde le deuxième ordre car dr est souvent 1, donc pas si petit
    g = ones(length(times),R,length(rr)-1)
    for i in 1:length(times) , r in 1:R
        normalisation[i,r] = 1 ./ aire_couronne / length(distances[i,r]) * aire_system
        g[i,r,:]  = fit(Histogram,distances[i,r],rr).weights .* normalisation[i,r]
    end
    gavg = mean(g,dims=2)[:,1,:]
    p=plot(xlabel="r",ylabel="g(r)")
    for i in 1:length(times)
        plot!(rr[1:end-1],gavg[i,:],rib=0,label="t=$(1+times[i])")
    end
    p
    plog=plot(xlabel="r",ylabel="g(r)")
        for i in 1:length(times)
            plot!(rr[2:end-1],gavg[i,2:end],rib=0,label="t=$(1+times[i])",axis=:log)
        end
        ylims!(0.5,10)
        plog


vortices,antivortices = spot_defects_separated(save_thetas,BC)
    plot(title="T = $(T) , σ² = $(Var)")
    scatter!(vortices,label="+")
    scatter!(antivortices,label="-")


## Analyse des données pour la g(r)
dico = JLD2.load("data/data_gr_0810.jld2")

@unpack L,init,init,BC,R,tmax,TV,times = dico["metadata"] # @unpack found in Parameters.jl
dmax = floor(Int,L/2) ; dr = 1 ; r = 0:dr:dmax

grt_same_charge     = zeros(length(TV),length(r)-1)
grt_opposite_charge = zeros(length(TV),length(r)-1)
grt_all_charge = zeros(length(TV),length(r)-1)
for i in eachindex(TV)
    T,Var = TV[i] ; group_name = "T$(T)_Var$Var"

    distances_opposite_charge = dico[group_name*"/opposite_charge"]
    distances_same_charge = dico[group_name*"/same_charge"]

    aire_system = π/4*L^2  # aire d'un cercle et non L² car on a filtered out les distances > L/2 lors de la génération des données
    aire_couronne = π*(2 .* r[1:end-1] * dr .+ dr^2) # on garde le deuxième ordre car dr est souvent 1, donc pas si petit
    normalisation = 1 ./ aire_couronne / length(distances_opposite_charge) * aire_system # so that g(r->\inf) = 1

    grt_opposite_charge[i,:]  = fit(Histogram,distances_opposite_charge,r).weights .* normalisation
    grt_same_charge[i,:]      = fit(Histogram,distances_same_charge,r).weights     .* normalisation
    grt_all_charge[i,:]       = fit(Histogram,vcat(distances_same_charge,distances_opposite_charge),r).weights     .* normalisation/2
end

debut = 1
    p1 = plot(r[debut:end-1],(grt_opposite_charge[1,debut:end]),title="Opposite Charge",xlabel="r")
    p2 = plot(r[debut:end-1],(grt_same_charge[1,debut:end]),title="Same Charge",xlabel="r")
    p3 = plot(r[debut:end-1],(grt_all_charge[1,debut:end]),title="All Charges",xlabel="r")
    debut = 2
    p4 = plot(r[debut:end-1],(grt_opposite_charge[1,debut:end]),axis=:log,m=:circle,xlabel="r")
    p5 = plot(r[debut:end-1],(grt_same_charge[1,debut:end]),axis=:log,m=:circle,xlabel="r")
    p6 = plot(r[debut:end-1],(grt_all_charge[1,debut:end]),axis=:log,m=:circle,xlabel="r")
    plot!(r[debut:10],3000r[debut:10].^-1.7)
    plot(p1,p2,p3,p4,p5,p6,size=(1200,800))
    # plot(p1,p2,p3,size=(1200,800))
# savefig("figures/gr_T0.7_Var0.pdf")

vortices,antivortices = spot_defects_separated(JLD2.load("data/data_gr_0810.jld2","T0.2_Var0.0/thetas_r1"),BC)
p1 = plot(title="T = $(TV[1][1]) , σ² = $(TV[1][2])")
    scatter!(vortices,label="+")
    p2 = plot(title="T = $(TV[1][1]) , σ² = $(TV[1][2])")
    scatter!(vortices,label="+")
    scatter!(antivortices,label="-")
    plot(p1,p2,size=(800,400))
# savefig("figures/visu_defects_L$(L)_T$((TV[1][1]))_Var$((TV[1][2])).pdf")

sum(grt_all_charge.-1)


## XY relation ξ ∼ 1/√n holds ?
L = 200
    T = 0.4
    Vars = [0.0,0.02,0.05,0.1,0.015,0.2]
    init = "hightemp"
    BC = "periodic"
    q = +1 ; r0 = Int(L/2) # dummy
    tmax = 200 ; times = logspace(1,tmax,60)
    R = 40

    rs = 1:Int(L/2)
    C = zeros(Int(L/2),length(times),length(Vars),R)
    ξ = zeros(length(times),length(Vars),R)
    nb_defects = zeros(length(times),length(Vars),R)

zt = @elapsed Threads.@threads for r in 1:R
    for i in eachindex(Vars)
        println()
        println("σ² = $(Vars[i]) and r=$r")
        t = 0.0 ; thetas,omegas,dt = initialize(L,T,Vars[i],init,q,r0)
        token = 1
        while t < tmax
            t += dt ; update(thetas,omegas,L,T,dt)
            if t>times[token]
                C[:,token,i,r] = get_C(thetas)
                ξ[token,i,r] = corr_length(rs,C[:,token,i,r] )
                nb_defects[token,i,r] = length(spot_defects(thetas,BC))
                token = min(length(times),1+token)
            end
        end
    end
end
# JLD2.jldsave("data\\xi_sqrtn_holds.jld2";L,init,BC,tmax,C,ξ,nb_defects,R,times,T,Vars)


p=plot()
    for i in eachindex(Vars)
        plot!(times,smooth(mean(ξ,dims=3)[:,i,1].* sqrt.(mean(nb_defects,dims=3)[:,i,1])/L,over=5),xaxis=:log)
        # plot!(times,ξ[:,i,2].* sqrt.(nb_defects[:,i,2])/L,xaxis=:log,m=true)
        # plot!(times,sqrt.(mean(nb_defects,dims=3)[:,i,1])/L,xaxis=:log,m=true)
        # plot!(times,mean(ξ,dims=3)[:,i,1]/L,xaxis=:log)
    end
    p


# # p = plot()
# #     for r in 1:R
# #         plot!(rs,C[:,50,1,r])
# #     end
# #     plot!(rs,(mean(C[:,50,1,:],dims=4)[:,1,1,1]),c=:black)
# #     p
#
# # à quoi ressemble xi ?
# plot(times,ξ[:,1,2],xaxis=:log,m=true)
#     over = 10
#     # plot!(times[1:end-over+1],ema(ξ[:,1,2],over)/1.2,xaxis=:log,m=true)
#
# # à quoi ressemble nb ?
# plot(times,nb_defects[:,1,2],xaxis=:log,m=true)
#     plot!(times,nb_defects[:,2,2],xaxis=:log,m=true)
#
# # à quoi ressemble nb ?
# plot(times[1:end-over+1],ema(ξ[:,1,2],over) .* sqrt.(nb_defects[1:end-over+1,1,2]),xaxis=:log,m=true)
#     plot!(times[1:end-over+1],ema(ξ[:,1,2],over) .* sqrt.(nb_defects[1:end-over+1,2,2]),xaxis=:log,m=true)

## Correlation matrix Omega and creation/annihilation locations
L = 200
T,Var = (0.2,0.2)
    # T,Var = (0.5,0.)
    Random.seed!(1234); omegas_fixed = rand(Normal(0,sqrt(Var)),L,L) # the seed is only active for this line as I understand it
    init = "hightemp"
    BC = "periodic"
    q = +1 ; r0 = Int(L/2) # dummy
    R = 3
    transients = 1500 ; tmax = 500
    track_every = 1 ; times_tracking = Float64.(collect(0:track_every:tmax))

# p = plot(xlims=(0,L),ylims=(0,L))
thetas,~,dt = initialize(L,T,Var,init,q,r0)
    omegas = omegas_fixed
    t = 0.0 ; while t < transients t += dt ; update(thetas,omegas,L,T,dt) end

    t = 0.0 ; token_t = 1 ;
    fiche_vortices,fiche_antivortices,vortices_old,antivortices_old = initialize_fiches_defects(thetas,BC,t)
    len = length(fiche_vortices)
    z = @elapsed while t < tmax
        t += dt
        update(thetas,omegas,L,T,dt)
        if round(t,digits=2) > times_tracking[token_t]
            vortices_new,antivortices_new = spot_defects_separated(thetas,BC)
            fiche_vortices,fiche_antivortices = defect_tracker(thetas,BC,fiche_vortices,fiche_antivortices,vortices_old,antivortices_old,times_tracking[token_t])

            # if length(fiche_vortices) > len
            #     scatter!(location_creation(fiche_vortices[len+1:end]),m=:cross,c=:black,title="t = $(round(times_tracking[token_t],digits=1))")
            #     display(p)
            #     len = length(fiche_vortices)
            # end

            vortices_old,antivortices_old = vortices_new,antivortices_new
            token_t += 1
        end
    end
println("Runtime : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")

p = plot(xlims=(0,L),ylims=(0,L),title="T = $T σ² = $Var")
    scatter!(location_creation(fiche_vortices),m=:cross,c=:black)
    scatter!(location_annihiliation(fiche_vortices),m=:cross,c=:red)
# savefig("figures\\correlation_creation_omegas_T$(T)_Var$Var.pdf")


# scatter!(spot_runaways(omegas_fixed,1.5,3),c=:green)

    scatter!(get_fastest_oscillators(omegas_fixed,threshold=1.5),c=:green)
&

## Screenshots pour estimer les impacts des différents niveaux de randomness sur le champ theta
L = 200
    T,Var = (0.2,0.2)
    Random.seed!(1234); omegas_fixed = rand(Normal(0,sqrt(Var)),L,L) # the seed is only active for this line as I understand it
    q = +1 ; r0 = Int(L/2) # dummy
    R = 4
    tmax = 1000
    BC = "periodic"

# Impact du bruit thermique (init lowtemp et omegas fixed)
    # -> Réponse : ca change pas grand chose à l'oeil nu
init = "lowtemp"
p = Vector{Any}(undef,R)
thetass = zeros(L,L,R)
z = @elapsed Threads.@threads for r in 1:4
    # println("r = $r")
    thetas,~,dt = initialize(L,T,Var,init,q,r0)
    t = 0.0
    while t<tmax
        t += dt ; update(thetas,omegas_fixed,L,T,dt)
    end
    thetass[:,:,r] = thetas
    p[r] = heatmap(mod.(thetas',2π),c=cgrad([:black,:blue,:green,:orange,:red,:black]),clims=(0,2π),colorbar=nothing)
end
plots=plot(p...,layout=(2,2),size=(800,800))
JLD2.jldsave("data/impact_noise_$(init)_fixed_omegas_fixed_T$(T)_Var$(Var)_tmax$(tmax).jld2";L,tmax,thetass,omegas_fixed,T,Var,plots,init,BC,p)
savefig("figures/figures_paper/impact_noise_$(init)_fixed_omegas_fixed_T$(T)_Var$(Var)_tmax$(tmax).png")

# Impact du bruit thermique (init hightemp fixed et omegas fixed)
    # -> Réponse : ca change pas grand chose à l'oeil nu
init = "hightemp"
Random.seed!(123); thetas_hightemp_fixed = 2pi*rand(L,L)
p = Vector{Any}(undef,R)
z = @elapsed Threads.@threads for r in 1:4
    # println("r = $r")
    dt = determine_dt(T,Var,L)
    thetas = thetas_hightemp_fixed
    t = 0.0
    while t<tmax
        t += dt ; update(thetas,omegas_fixed,L,T,dt)
    end
    p[r] = heatmap(mod.(thetas',2π),c=cgrad([:black,:blue,:green,:orange,:red,:black]),clims=(0,2π),colorbar=nothing)
    # highlight_defects(p[r],thetas)
end
plots = plot(p...,layout=(2,2),size=(800,800))
JLD2.jldsave("data/impact_noise_$(init)_fixed_omegas_fixed_T$(T)_Var$(Var)_tmax$(tmax).jld2";L,tmax,omegas_fixed,T,Var,plots,init,BC,p)
savefig("figures/figures_paper\\impact_noise_$(init)_fixed_omegas_fixed_T$(T)_Var$(Var)_tmax$(tmax).png")

# Impact conjoint de l'init et du bruit thermique (init hightemp variable et omegas fixed)
    # -> Réponse, les vortex créés à hightemp ne sont pas au meme endroit pour les memes réalisations donc normal qu'on ait un champ theta différent au final
init = "hightemp"
p = Vector{Any}(undef,R)
z = @elapsed Threads.@threads for r in 1:4 # jsais pas pk mais ca supporte pas le Threads.@threads , ca shutdown julia
    # println("r = $r")
    thetas,~,dt = initialize(L,T,Var,init,q,r0)
    t = 0.0
    while t<tmax
        t += dt ; update(thetas,omegas_fixed,L,T,dt)
    end
    p[r] = heatmap(mod.(thetas',2π),c=cgrad([:black,:blue,:green,:orange,:red,:black]),clims=(0,2π),colorbar=nothing)
end
plots = plot(p...,layout=(2,2),size=(800,800))
JLD2.jldsave("data/impact_noise_$(init)_variable_omegas_fixed_T$(T)_Var$(Var)_tmax$(tmax).jld2";L,tmax,omegas_fixed,T,Var,plots,init,BC,p)
savefig("figures/figures_paper\\impact_noise_$(init)_variable_omegas_fixed_T$(T)_Var$(Var)_tmax$(tmax).png")


## Dynamic Coarsening
@unpack init,C,tsave,tmax,Ts,Vars = JLD2.load("data/scalingXI_Var_L500_HighTemp.jld")

# load("data/XY_saturationXI_L200_HighTemp.jld") # violet
# load("data/scalingXI_Var_L100_HighTemp.jld")
# load("data/scalingXI_Var_L200_HighTemp.jld")

seuil = exp(-1)
variances = [0.0,0.02,0.06,0.1,0.2]
Ls = [200,500]
ξ = zeros(length(Ls),length(variances),39)
for i in eachindex(Ls)
    println(i)
    L = Ls[i] ; Lover2 = round(Int,L/2,RoundDown) ; r_vector = sort(vcat([n for n in 1:Lover2],[n*sqrt(2) for n in 1:Lover2]))
    C = load("data/scalingXI_Var_L$(L)_HighTemp.jld","C")
    C_avg = mean(C,dims=5)[1,[1,3,7,11,21],:,:,1]
    for j in 1:length(variances)
        for n in 1:39
            # try
                # Linear interpolation
                i_after = findfirst(x -> x < seuil,C_avg[j,:,n])
                i_befor = i_after - 1
                r_after = r_vector[i_after]
                r_befor = r_vector[i_befor]
                c_after = C_avg[j,i_after,n]
                c_befor = C_avg[j,i_befor,n]
                ξ[i,j,n] = (seuil*(r_after-r_befor) -  (c_befor*r_after - r_befor*c_after))/(c_after-c_befor)
            # catch e
            #     println(e)
            #     ξ[i,j,n] = NaN
            # end
        end
    end
end

colors = [:purple,ColorSchemes.tab10[1],ColorSchemes.tab10[3],:darkgoldenrod1,ColorSchemes.tab10[4]]
p = plot(xlabel="t",ylabel="ξ(t)",legend=:bottomright)
    plot!([NaN],[NaN],c=:grey,m=:star,label="L = 200",line=false)
    plot!([NaN],[NaN],c=:grey,m=:utriangle,label="L = 500",line=false)
    for j in eachindex(variances)
        # plot!([NaN],[NaN],c=colors[j],rib=0,label="σ² = $(variances[j])")
        plot!(tsave[1:39],ξ[1,j,:],m=:star,axis=:log,c=colors[j];line=false)
        plot!(tsave[1:39],ξ[2,j,:],m=:utriangle,axis=:log,c=colors[j];line=false)
    end
    annotate!(70,12,text(L"\sim \sqrt{t\,/\,\ln t}",40.0,10))
    ylims!(0.8,1.2maximum(ξ))
    xlims!(0.8,1500)
    yticks!(vcat(1:10,20),string.([1,2,3,4,5,"","","","",10,20]))
    xticks!([1,10,100,1000],["10^0","10^1","10^2","10^3"])
    plot!(tsave[10:37],2.3 .* sqrt.(tsave[10:37] ./ log.(tsave[10:37])),c=:black)
    p
# savefig("figures/dynamic_coarsening_rest_legend.svg")

## Time to reach SS versus sigma
using LambertW,QuadGK
@unpack L,init,C,tsave,tmax,Ts,Vars = JLD2.load("data/scalingXI_Var_L200_HighTemp.jld")
Lover2 = round(Int,L/2,RoundDown) ; r_vector = sort(vcat([n for n in 1:Lover2],[n*sqrt(2) for n in 1:Lover2]))
C_avg = mean(load("data/scalingXI_Var_L$(L)_HighTemp.jld","C"),dims=5)[1,:,:,:,1]
ξ = zeros(length(Vars),39) ; seuil = exp(-1)
for j in 1:length(Vars)
    for n in 1:39
        # try
            # Linear interpolation
            i_after = findfirst(x -> x < seuil,C_avg[j,:,n])
            i_befor = i_after - 1
            r_after = r_vector[i_after]
            r_befor = r_vector[i_befor]
            c_after = C_avg[j,i_after,n]
            c_befor = C_avg[j,i_befor,n]
            ξ[j,n] = (seuil*(r_after-r_befor) -  (c_befor*r_after - r_befor*c_after))/(c_after-c_befor)
        # catch e
        #     println(e)
        #     ξ[i,j,n] = NaN
        # end
    end
end

# Calcul de la longueur caractérstique du petit argument 1D et noiseless du PRL
Varss = 0.01:0.01:0.25
esp = zeros(length(Varss))
for i in eachindex(Varss)
    σ = sqrt(Varss[i])
    # f(x) = asin(σ*x)*exp(-x^2)
    # esp[i] = quadgk(f,0,sin(1)/σ)[1]
    g(x) = cos(x)*x*exp(-(sin(x)/σ)^2)*2/σ/sqrt(pi)
    esp[i] = quadgk(g,0,pi/2)[1]
end

# # calcul du temps caractérstique pour atteindre la longueur caractérstique finale (au moins 10% de diff avec le cas XY)
# xi_XY = ξ[1,:]
# diff_to_xy = zeros(size(ξ))
# for j in 2:length(Vars)
#     for t in 1:39
#         diff_to_xy[j,t] = abs.(ξ[j,t] - ξ[1,t]) ./ ξ[1,t]
#     end
# end
# tau = Vector(undef,length(Vars))
# for i in 1:length(Vars)
#     ind_tmp = findfirst(x->x>0.1,diff_to_xy[i,:])
#     if ind_tmp ≠ nothing tau[i] = tsave[ind_tmp] end
# end

# calcul du temps nécessaire pour que le cas XY atteigne la valeur \ell = pi/2/esp
tau = Vector(undef,length(Vars))
for i in 1:length(Vars)
    ind_tmp = findfirst(x->x≥π/2/esp[i],ξ[1,:])
    if ind_tmp ≠ nothing tau[i] = tsave[ind_tmp] end
end

    plot(size=(200,200))
    plot!(Vars[3:end],tau[3:end],axis=:log,m=:circle,line=false,c=:grey)
    plot!(Varss[2:end],0.13exp.(-lambertw.(-(pi ./2esp[2:end]).^-2,-1)),c=:red)
    xticks!(vcat(0.01:0.01:0.1,0.2:0.1:0.5),string.([L"10^{-2}","","","","","","","","",L"10^{-1}","","","",""]))
    yticks!(vcat(10:10:100,200:100:1000),string.([L"10^{1}","","","","","","","","",L"10^{2}","","","","","","","","",L"10^{3}"]))
    xlims!(0.01,0.6)
    ylims!(10,1000)
    annotate!(.39,15,text(L"σ²",12))
    annotate!(0.015,600,text(L"τ",12))

# savefig("figures\\figures_paper\\scaling_transients.svg")

## Analyse Data Lifetime
JLD2.load("data/phase_space_lifetime_pair.jld")
lifetime = mean(JLD2.load("data/phase_space_lifetime_pair.jld","lifetime"),dims=3)[:,:,1]
Vars = JLD2.load("data/phase_space_lifetime_pair.jld","Vars")
Ts = JLD2.load("data/phase_space_lifetime_pair.jld","Ts")
heatmap(Vars,Ts,mean(lifetime,dims=3)[:,:,1],c=cgrad([:red,:darkgreen,:blue,:white]),xlabel = "σ²",ylabel = "T")
    Tkt = 0.8 ; xx= Array(0:0.0001:7 /log(200)^2) ; plot!(xx,Tkt*(1 .- sqrt.(xx)*log(200)/sqrt(1/0.43^2)),c=:black,lw=2)
    xlims!(0,maximum(Vars))
    ylims!(0,maximum(Ts))

# savefig("figures\\phase_space_lifetime_pair.png")

## Nombre total de vortex : phase space and over time
@unpack L,init,BC,tmax,Ts,Vars,times,nb_vortices = load("data/phase_space_number_vortices_lowtemp.jld")
nb_vortices_avg = mean(nb_vortices,dims=4)[:,:,:,1]
# Phase Space
anim = @animate for i in 1:size(nb_vortices_avg,3)
    heatmap(Vars,Ts,log.(1 .+ nb_vortices_avg[:,:,i]),title="t=$(times[i])",c=cgrad([:blue,:green,:orange,:red]),clims=(0,log.(1 .+ maximum(nb_vortices_avg[:,:,end]))))
    Tkt = 0.8 ; xx= Array(0:0.0001:7 /log(200)^2) ; plot!(xx,Tkt*(1 .- sqrt.(xx)*log(200)/sqrt(1/0.43^2)),c=:black,lw=2)
    xlims!((0,maximum(Vars)))
    ylims!((0,maximum(Ts)))
end
mp4(anim,"figures/films/phase_space_nb_vortices_$(init).mp4",fps=5);

# Final time
@unpack Ts,Vars,nb_vortices = load("data/phase_space_number_vortices_hightemp.jld")
tmax = 3000 ; tmax_complement = 15000
nb_final_time = nb_vortices[:,:,end,:]
nb_vortices_complement = load("data/phase_space_number_vortices_complement_hightemp.jld","nb_vortices")
nb_final_time[1:5,:,:] = nb_vortices_complement[:,:,end,:]
nb_vortices_avg = mean(nb_final_time,dims=3)[:,:,1]
p1 = heatmap(Vars,Ts,log.(1 .+ nb_vortices_avg),title="t=$(tmax)",c=cgrad([:blue,:green,:orange,:red]),clims=(0,log.(1 .+ maximum(nb_vortices_avg[:,:,end]))))
p1 = heatmap(Vars,Ts,log.(1 .+ nb_vortices_avg),title="t=$(tmax_complement)",c=cgrad([:blue,:green,:orange,:red]),clims=(0,log.(1 .+ maximum(nb_vortices_avg[:,:,end]))))

# n(t)
p=plot(xlabel="t",ylabel=L"n_v(t)",legend=nothing)
    token = 1
    plot!(times[2:end],log.(1 .+ nb_vortices_avg[6,1,2:end]),xaxis=:log,linewidth=2,c=:purple,label="σ² = 0")
    for j in ([6,11,16,21,26])
        plot!(times[2:end],log.(1 .+ nb_vortices_avg[9,j,2:end]),xaxis=:log,line=true,m=:circle,ms=2.5,c=ColorSchemes.tab10.colors[token],label="σ² = $(Vars[j])")
        token += 1
    end
    p
    # xticks!([1,10,100,1000],["10^0","10^1","10^2","10^3"])
    # yticks!([1,10,100,1000,10000],["10^0","10^1","10^2","10^3","10^4"])
    # plot!(times[15:end],1E3 * log.(times[15:end]) ./times[15:end],c=:black)
    # annotate!(350,8,text(L"\sim t\,/\,\ln t",-39.0,10))
    # savefig("figures/decay_number_defects.pdf")

## Comparaison nombre vortices lowtemp et hightemp
@unpack L,init,BC,tmax,Ts,Vars,times,nb_vortices = load("data/phase_space_number_vortices_lowtemp.jld")
nb_vortices_low_avg = mean(nb_vortices,dims=4)[:,:,:,1]
@unpack L,init,BC,tmax,Ts,Vars,times,nb_vortices = load("data/phase_space_number_vortices_hightemp.jld")
nb_vortices_high_avg = mean(nb_vortices,dims=4)[:,:,:,1]
p=plot(xlabel="t",ylabel=L"\log_{10}(1+n_v)",legend=nothing)
    plot!(times[1:end],log10.(1 .+ nb_vortices_low_avg[30,1,1:end]),xaxis=:log,line=true,m=:circle,ms=2.5,c=3,label="σ² = $(Vars[21])")
    plot!(times[1:end],log10.(1 .+ nb_vortices_high_avg[30,1,1:end]),xaxis=:log,line=true,m=:circle,ms=2.5,c=3,label="σ² = $(Vars[21])")
    plot!(times[1:end],log10.(1 .+ nb_vortices_low_avg[9,21,1:end]),xaxis=:log,line=true,m=:circle,ms=2.5,c=5,label="σ² = $(Vars[21])")
    plot!(times[1:end],log10.(1 .+ nb_vortices_high_avg[9,21,1:end]),xaxis=:log,line=true,m=:circle,ms=2.5,c=5,label="σ² = $(Vars[21])")
    xticks!([0.1,1,10,100,1000],["10^{-1}","10^0","10^1","10^2","10^3"])

## Time to reach SS with Vars
@unpack L,init,BC,tmax,Ts,Vars,times,nb_vortices = load("data/phase_space_number_vortices_hightemp.jld")
    nb_vortices_low_avg = mean(nb_vortices,dims=4)[:,:,:,1]

times_half = NaN*zeros(length(Vars))
temp = 20
for j in eachindex(times_half)
    maxi = nb_vortices_low_avg[temp,j,end]
    times_half[j] = times[findfirst(x->x≥maxi/2,nb_vortices_low_avg[temp,j,:])]
end
p = plot()
    plot!(Vars[1:end],nb_vortices_low_avg[temp,:,end])
    plot!(Vars[1:end],30 ./ (Vars[1:end]))
    xlims!(0.01,0.3)


## Mean Lifetime Phase space
L = 200
    Ts = collect(0:0.025:0.6)
    Vars = collect(0.0:0.01:0.2)
    Ts = [0.]
    Vars = [0.1,0.11,0.12]
    init = "pair"
    BC = "periodic"
    transients = 200 ; tmax = 1000
    track_every = 2 ; times_tracking = Float64.(collect(0:track_every:tmax))
    q = +1 ; r0 = 20

    lifetime = Matrix{Union{Float64,Missing}}(missing,length(Ts),length(Vars))# dummy initialisation

z = @elapsed for i in eachindex(Ts) , j in eachindex(Vars)
    T = Ts[i] ; Var = Vars[j]
    println("i = $i/$(length(Ts)) et j = $j/$(length(Vars))")
        thetas,omegas,dt = initialize(L,T,Var,init,q,r0)

        t = 0.0 ;
        while t < transients t += dt ; update(thetas,omegas,L,T,dt) end

        t = 0.0 ; token_t = 1 ;
        fiche_vortices,fiche_antivortices,vortices_old,antivortices_old = initialize_fiches_defects(thetas,BC,t)
        while t < tmax
            t += dt
            update(thetas,omegas,L,T,dt)
            if round(t,digits=2) > times_tracking[token_t]
                vortices_new,antivortices_new = spot_defects_separated(thetas,BC)
                fiche_vortices,fiche_antivortices = defect_tracker(thetas,BC,fiche_vortices,fiche_antivortices,vortices_old,antivortices_old,times_tracking[token_t])

                if fiche_vortices[1][4] ≠ nothing
                    lifetime[i,j,r] = round(t,digits=2)
                    break
                end

                vortices_old,antivortices_old = vortices_new,antivortices_new
                token_t += 1
            end
        end
    if lifetime[i,j] === missing
        lifetime[i,j] = tmax
    end
end
JLD2.jldsave("data\\TV_lifetime.jld2";L,Ts,Vars,init,BC,tmax,transients,track_every,times_tracking,q,r0,R,lifetime)
println("Runtime : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")

lifetime_avg = mean(lifetime,dims=3)[:,:,1]
heatmap(Vars,Ts,lifetime_avg)
heatmap(lifetime_avg)


##
TT = 0.0
VVar = 0.1
thetas,omegas,dt = initialize(L,TT,VVar,init,q,r0)

t = 0.0 ;
while t < 100 t += dt ; update(thetas,omegas,L,TT,dt) end

t = 0.0 ; token_t = 1 ;
fiche_vortices,fiche_antivortices,vortices_old,antivortices_old = initialize_fiches_defects(thetas,BC,t)
while t < tmax
    t += dt
    update(thetas,omegas,L,TT,dt)
    if round(t,digits=2) > times_tracking[token_t]
        vortices_new,antivortices_new = spot_defects_separated(thetas,BC)
        fiche_vortices,fiche_antivortices = defect_tracker(thetas,BC,fiche_vortices,fiche_antivortices,vortices_old,antivortices_old,times_tracking[token_t])

        if fiche_vortices[1][4] ≠ nothing
            lifetime[i,j,r] = round(t,digits=2)
            break
        end

        vortices_old,antivortices_old = vortices_new,antivortices_new
        token_t += 1
    end
end
if lifetime[i,j] === missing
lifetime[i,j] = tmax
end

## Tracking et identification des vortex IN REAL TIME
    # Le but est de connaître le temps de vie des vortex pour des systèmes
    # avec plusieurs (mais pas trop) de vortex.
L = 200
    T = 0.7; Var = 0.
    init = "lowtemp"
    BC = "periodic"
    transients = 100 ; tmax = 500
    track_every = 1 ; times_tracking = Float64.(collect(0:track_every:tmax))
    q = +1 ; r0 = 30

thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
    # thetas = thetas .+ initialize(L,T,Var,init,q,5r0)[1]'

    t = 0.0 ; token_t = 1 ;
    ztr = @elapsed while t < transients t += dt ; update(thetas,omegas,L,T,dt) end
    println("Transients time : $(round(Int,ztr)) seconds = $(round(ztr/60,digits=2)) minutes.")
    # p = display_thetas(thetas)
    #     highlight_defects(p,thetas)

    fiche_vortices,fiche_antivortices,vortices_old,antivortices_old = initialize_fiches_defects(thetas,BC,t)

    z = @elapsed while t < transients + tmax
        t += dt
        update(thetas,omegas,L,T,dt)
        if round(t,digits=2) > times_tracking[token_t] + transients
            vortices_new,antivortices_new = spot_defects_separated(thetas,BC)
            fiche_vortices,fiche_antivortices = defect_tracker(thetas,BC,fiche_vortices,fiche_antivortices,vortices_old,antivortices_old,times_tracking[token_t])
            vortices_old,antivortices_old = vortices_new,antivortices_new
            token_t += 1
        end
    end
println("Runtime : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")
p = display_thetas(thetas)
    highlight_defects(p,thetas)

xx = 0:track_every:30
    # plot(xlabel="Lifetime",ylabel="")
    histogram!(lifetimes(fiche_vortices,tmax),normalize=true,yaxis=:log,bins=xx)
    plot!(xx,exp.(-0.3xx),c=:black)

## Benchmark Tracking
L = 200
    T = 0.1; Var = 0.
    init = "lowtemp"
    BC = "periodic"
    transients = 100 ; tmax = 500.
    track_every = 0.5 ; times_tracking = collect(transients:track_every:tmax)
    q = +1 ; r0 = 12


thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
# thetas = thetas .+ initialize(L,T,Var,init,q,5r0)[1]'
t = 0.0
fiche_vortices,fiche_antivortices,vortices_old,antivortices_old = initialize_fiches_defects(thetas,BC,t)
for i in 1:2000 update(thetas,omegas,L,T,dt) end
vortices_new,antivortices_new = spot_defects_separated(thetas,BC)
antivortices_new == antivortices_old
fiche_vortices,fiche_antivortices = defect_tracker(thetas,BC,fiche_vortices,fiche_antivortices,vortices_old,antivortices_old,tmax)
p = display_thetas(thetas)
    highlight_defects(p,thetas)

## Evolution temporelle L = 500,1000 pour comprendre cette g(r) trop piquée
L = 1000
    T = 0.7 ; Var = 0
    init = "lowtemp"
    BC = "periodic"
    tmax = 2000
    q = +1 ; r0 = Int(L/2)

thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
t = 0.0
    z = @elapsed while t < tmax
        t += dt
        update(thetas,omegas,L,T,dt)
    end
println("Runtime : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")

save("data/L1000_tmax$(tmax)_$init.jld","thetas",thetas,"omegas",omegas,"tmax",tmax,"T",T,"Var",Var,"init",init,"dt",dt)
thetas = load("data/L1000_tmax2000.jld","thetas") ; L = 1000 ;
vortices,antivortices = spot_defects_separated(thetas,"periodic")
p2 = plot()
    scatter!(vortices,label="+1")
    scatter!(antivortices,label="-1")

p1 = display_thetas(thetas)


m = length(vortices)
println() ; println("Nombre vortices = $m")
distances_opposite_charge = Vector{Float64}()
for i in 1:m , j in 1:m
    push!(distances_opposite_charge,dist(vortices[i],antivortices[j],L))
end

dr = 1
    r = 0:dr:L/2
    aire_couronne = π*(2 .* r[1:end-1] * dr .+ dr^2)

    normalisation = 1 ./ aire_couronne * L^2 / length(distances_opposite_charge) # so that g(r->\inf) = 1

    grt_opposite_charge = fit(Histogram,distances_opposite_charge,r).weights .* normalisation
    p3 = plot(r[2:end-1],grt_opposite_charge[2:end],axis=:log)

plot(p1,p2,p3,layout=(1,3),size=(1200,400),title="T=0.7, LowTemp")
savefig("figures\\vortices_T=0.7_LowTemp.png")

## Generation de données pour les g(r) (fonctionnel)
L = 200
    init = "hightemp"
    init == "isolated" ? BC = "free" : BC = "periodic"
    R = 20
    tmax = 200 ; t_relax = 200
    times = collect(10:10:tmax)
    q = +1 ; r0 = Int(L/2) # dummy

    # TV = [(0.7,0.0) , (1.0,0.0) , (0.2,0.05) , (0.3,0.1) , (0.6,0.15) , (1.0,0.1) , (0.0,0.3)]
    TV = [(0.2,0.0)]


file = jldopen("data/data_gr_0810.jld2","w")
    file["metadata"] = Dict("L"=>2L , "init"=>init , "t_relax"=>t_relax , "BC"=>BC , "R"=>R , "tmax"=>tmax , "TV"=>TV , "times"=>times)
    z = @elapsed for (T,Var) in TV
    println("Simulation for T = $T and σ² = $Var started.")
    group_name = "T$(T)_Var$Var"

    distances_same_charge     = Vector{Float64}()
    distances_opposite_charge = Vector{Float64}()

    Threads.@threads for real in 1:R
        thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
        t = 0.0
        while t < tmax - t_relax
            t += dt
            update(thetas,omegas,L,T,dt)
        end

        thetas_ext = zeros(2L,2L)
        thetas_ext[1:L,1:L] = thetas_ext[L+1:end,1:L] = thetas_ext[1:L,L+1:end] = thetas_ext[L+1:end,L+1:end] = thetas
        omegas_ext = zeros(2L,2L)
        omegas_ext[1:L,1:L] = omegas_ext[L+1:end,1:L] = omegas_ext[1:L,L+1:end] = omegas_ext[L+1:end,L+1:end] = omegas

        while t < tmax
            t += dt
            update(thetas_ext,omegas_ext,2L,T,dt)
        end
        if real == 1 file[group_name*"/thetas_r1"] = thetas_ext end

        vortices,antivortices = spot_defects_separated(thetas_ext,BC)
        m = length(vortices)
        println() ; println("Number of vortices = $m for realisation $real")

        # Same charge
        distances = Float64[]
        for i in 1:m , j in 1+i:m
            push!(distances,dist(vortices[i],vortices[j],2L))
            push!(distances,dist(antivortices[i],antivortices[j],2L))
        end
        push!(distances_same_charge,filter!(x->x≤2L/2,distances)...) # les "..." sont là pour ajouter les valeurs et non le vecteur

        # Opposite charge
        distances = Float64[]
        for i in 1:m , j in 1:m
            push!(distances,dist(vortices[i],antivortices[j],2L))
        end
        push!(distances_opposite_charge,filter!(x->x≤2L/2,distances)...)
    end # realisations

        file[group_name*"/same_charge"]     = distances_same_charge
        file[group_name*"/opposite_charge"] = distances_opposite_charge

        println("Simulation for T = $T and σ² = $Var finished.")
    end
    close(file)
println("Runtime : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")


## Progress Bar Atom
iter = ProgressBar(1:100)
       for i in iter
          sleep(0.1)
          loss = 1/i
          set_description(iter, string(println("Loss: %.2f", loss)))
       end

## Generation des données pour les films de la nucléation
L = 200
    # TV = [(0.1,0.0),(0.2,0.0),(0.3,0.0),(0.2,0.05),(0.2,0.1),(0.2,0.15)]
    TV = [(0.1,0.1)]
    init = "hightemp"
    BC = "periodic"
    tmax = 400
    q = +1 ; r0 = Int(L/2) # dummy

file = jldopen("data/films_thetas_$(init).jld2","w")
   file["metadata"] = Dict("L"=>L , "init"=>init , "BC"=>BC , "tmax"=>tmax , "TV"=>TV)
   zs = @elapsed for i in eachindex(TV)
       T,Var = TV[i]
       println("Simulation for T = $T , Var = $Var started.")
       thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
       thetas_hist = zeros(Float16,L,L,1+round(Int,tmax/dt))
       thetas_hist[:,:,1] = Float16.(thetas)

       # Time series
       for i in 1:round(Int,tmax/dt)
           update(thetas,omegas,L,T,dt)
           thetas_hist[:,:,i+1] = Float16.(thetas)
       end
       file["thetas_history_T$(T)_Var$Var"] = thetas_hist
   end
   close(file)

println("Runtime : $(round(Int,zs)) seconds = $(round(zs/60,digits=2)) minutes.")

## Analyse de la nucléation pour films mp4
filename = "data/films_T0.05_lowtemp.jld"
TV,tmax,init,BC,tmax = JLD2.load(filename,"metadata").values

i = 4 ; println("TV = $(TV[i])") ; T,Var = TV[i] ; dt = determine_dt(T,Var,L) ;
    thetas_history = Float64.(JLD2.load(filename,"thetas_history_T$(T)_Var$Var"))
    # Compute theta dot
    thetas_dot_history = NaN*zeros(size(thetas_history))
    thetas_history_ema = NaN*zeros(size(thetas_history))
    window_ema = 50 # arbitrary, robustness to be proven
    tau = 50 # dt in finite difference scheme for first order derivative
    z = @elapsed for j in 1:L , i in 1:L
        tmp = ema(thetas_history[i,j,:],window_ema,wilder=true)
        thetas_history_ema[i,j,1:length(tmp)] = tmp
        for k in 1:length(tmp)-tau
            thetas_dot_history[i,j,k] = arclength(thetas_history_ema[i,j,k],thetas_history_ema[i,j,k+tau])/tau/dt
        end
        tmp = ema(thetas_dot_history[i,j,:],window_ema,wilder=true)
        thetas_dot_history[i,j,1:length(tmp)] = tmp
        # println(thetas_dot_history[1,1,1])
    end
    println("Computation time : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")

# Visualisation temporal series
# i,j = rand(1:100),rand(1:100)
# tt = collect(1:length(thetas_history[i,j,1:end]))*dt
#     p1 = plot(xlabel="t",ylabel="θ(t)",legend=:best)
#     plot!(tt,thetas_history[i,j,1:end],rib=0,label="Original Signal")
#     plot!(tt,thetas_history_ema[i,j,1:end],rib=0,label="Smoothed Signal")
#     # ylims!(-0.5,1.5)
#     p2 = plot(xlabel="t",ylabel=L"\dot{\theta}(t)")
#     plot!(tt,thetas_dot_history[i,j,1:end],c=:black)
#     # ylims!(-0.008,0.004)
#     plot(p1,p2,size=(800,400))
# savefig("figures\\sample_trajectory_time_theta_thetadot_6.pdf")

# plot(thetas_dot_history[i,j,1:end])
#     plot!(ema(thetas_dot_history[i,j,1:end],50,wilder=true))

tt = 45 ;
pAtheta = heatmap(mod.(thetas_history[:,:,round(Int,tt/determine_dt(0,0))]',2pi),c=cols,clims=(0,2π),colorbar=false,ticks=false)
pAthetadot = heatmap(thetas_dot_history[:,:,round(Int,tt/determine_dt(0,0))]',colorbar=false,clims=(-0.02,0.02),ticks=false)
tt = 447 ;
pBtheta = heatmap(mod.(thetas_history[:,:,round(Int,tt/determine_dt(0,0))]',2pi),c=cols,clims=(0,2π),colorbar=false,ticks=false)
pBthetadot = heatmap(thetas_dot_history[:,:,round(Int,tt/determine_dt(0,0))]',colorbar=false,clims=(-0.02,0.02),ticks=false)
tt = 568 ;
pCtheta = heatmap(mod.(thetas_history[:,:,round(Int,tt/determine_dt(0,0))]',2pi),c=cols,clims=(0,2π),colorbar=true,ticks=false,colorbar_title="θ",size=(450,400))
highlight_defects(pCtheta,thetas_history[:,:,round(Int,tt/determine_dt(0,0))])
pCthetadot = heatmap(thetas_dot_history[:,:,round(Int,tt/determine_dt(0,0))]',colorbar=true,clims=(-0.02,0.02),ticks=false,colorbar_title=L"d\theta/dt",size=(480,400))

ptheta = plot(pAtheta,pBtheta,pCtheta,layout=(1,3),size=2.(1250,400))
savefig("figures\\figures_paper\\snapshots_theta.png")
savefig("figures\\figures_paper\\snapshots_theta.svg")
pthetadot = plot(pAthetadot,pBthetadot,pCthetadot,layout=(1,3),size=(1250,400))
savefig("figures\\figures_paper\\snapshots_thetadot.pdf")
savefig("figures\\figures_paper\\snapshots_thetadot.svg")


# Film
# fin = findfirst(isnan,thetas_dot_history[1,1,:])-1
#     if fin == nothing fin = length(thetas_dot_history) end
#     zm = @elapsed anim = @animate for i in 10:10:fin
#         p1 = heatmap(mod.(thetas_history[:,:,i]',2π),c=cgrad([:black,:blue,:green,:orange,:red,:black]),clims=(0,2π),size=(600,512))
#         # p2 = heatmap(thetas_dot_history[:,:,i]',size=(600,512))
#         p2 = heatmap(thetas_dot_history[:,:,i]',size=(600,512),clims=(-0.05,0.05))
#         p = plot(p1,p2,size=(1200,512))
#         title!("L = $L , σ² = $Var , T=$T , t=$(round(dt*i,digits=2))")
#     end
#     println("Time Movie : $(round(Int,zm)) seconds = $(round(zm/60,digits=2)) minutes.")
#     # mp4(anim,"figures/films/test_ema_$(init)_T$(T)_Var$(Var)_freescale.mp4",fps=10);
#     mp4(anim,"figures/films/exp_moving_avg_$(init)_T$(T)_Var$(Var)_fixedscale.mp4",fps=10);

## Benchmark : robustness of tau and window_ema
@unpack L,init,BC,tmax,TV = JLD.load("data/films_T0_lowtemp.jld","metadata")

i = 3 ; println("TV = $(TV[i])") ; T,Var = TV[i] ; dt = determine_dt(T,Var,L) ;
thetas_history = Float64.(JLD.load("data/films_T0_lowtemp.jld","thetas_history_T$(T)_Var$Var"))
thetas_dot_history = NaN*zeros(size(thetas_history))
thetas_history_ema = NaN*zeros(size(thetas_history))
window_ema = 50 # arbitrary, robustness to be proven
    tau = 60 # dt in finite difference scheme for first order derivative

    for j in 1:10 , i in 1:10
        tmp = ema(thetas_history[i,j,:],window_ema,wilder=true)
        thetas_history_ema[i,j,1:length(tmp)] = tmp
        for k in 1:length(tmp)-tau
            thetas_dot_history[i,j,k] = arclength(thetas_history_ema[i,j,k],thetas_history_ema[i,j,k+tau])/tau/dt
        end
        tmp = ema(thetas_dot_history[i,j,:],window_ema,wilder=true)
        thetas_dot_history[i,j,1:length(tmp)] = tmp
    end

    i = rand(1:10) ; j = rand(1:10)
        p1 = plot(thetas_history[i,j,1:end])
        plot!(thetas_history_ema[i,j,1:end])
        # ylims!(-0.5,1.5)
        p2 = plot(thetas_dot_history[i,j,1:end],c=:black)
        # ylims!(-0.008,0.004)
        plot(p1,p2,layout=(2,1),size=(400,800))

## Technique pour dériver des données bruitées
dt = 0.01
t = 0:dt:2pi ;
w = 2pi*fftfreq(length(t),1/dt)
# signal = cos.(t).*exp.(-t.^2) + 0.0*randn(length(t))
signal = cos.(t) + 0.1*randn(length(t))

plot(w,real((fft(signal))))
signal_filtered = fftshift(fft(signal))
signal_filtered[1:round(Int,length(signal_filtered)/4)] .= 0
signal_filtered[round(Int,length(signal_filtered)*3/4):end] .= 0
signal_filtered = ifftshift(signal_filtered)
plot!(w,real(lowpass(fft(signal),w)))
xlims!(-5,5)

function lowpass(fft_signal,w,cutoff_w=2π/3(diff(w)[1]*2pi/length(w)))
    shifted_signal = fftshift(fft_signal)
    shifted_w = collect(fftshift(w))
    ind_neg = findfirst(x->x>-cutoff_w,shifted_w)-1
    ind_pos = findfirst(x->x>cutoff_w,shifted_w)
    shifted_signal[1:ind_neg] .= 0
    shifted_signal[ind_pos:end] .= 0
    return ifftshift(shifted_signal)
end

plot(t,signal)
    plot!(t,real(ifft(lowpass(fft(signal),w))))

plot(t,signal)
    plot!(t,real(ifft(im * w .* fft(signal))))

## Structure factor S(q) = 1 + \rho FT[g(r)]
i = 6
    xx = grt_opposite_charge[i,1:100]*mean([length(distances_same_charge[i]) for i in eachindex(distances_same_charge)])/L^2
    Sq = 1 .+ real.(fft(xx))
# Sq = 1 .+ real.(fft(ones(100)))
    plot(Sq[2:end],title="T = $(TV[i][1]) , σ² = $(TV[i][2])")

## g(r,t=tmax), Test
L = 200
    T = 0.1
    Var = 0.1
    init = "hightemp"
    init == "isolated" ? BC = "free" : BC = "periodic"
    q = +1
    r0 = Int(L/2)
    R = 12


    dmax = floor(Int,L/2) # maximum distance
    dr = 0.25 ; r = 1:dr:dmax

grt_p = zeros(length(r)-1,R)

token = 0
Threads.@threads for real in 1:R
    thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
    t = 0.0 ; tmax = 300
        z = @elapsed while t < tmax
            t += dt
            update(thetas,omegas,L,T,dt)
        end
    #     display_thetas(thetas)

    vortices,antivortices = spot_defects_separated(thetas,BC)
    m = length(vortices)

    count_p = zeros(m,length(r)-1) # g(r,t) for vortices with same charge
    for i in 1:m
        distances = Float64[]
        for j in 1:m
            # push!(distances,dist(vortices[i],vortices[j],L))
            # push!(distances,dist(antivortices[i],antivortices[j],L))
            push!(distances,dist(vortices[i],antivortices[j],L))
        end
        h = fit(Histogram,distances,r)
        count_p[i,:] = h.weights
    end
    g_gas = m/L^2
    grt_p[:,real] = mean(count_p,dims=1)[1,:]/π ./ (2 .* r[1:end-1] * dr .+ dr^2) / g_gas

    token += 1
    println("r = $token/$R")
end
plot(r,vec(mean(grt_p,dims=2)),title="T=$T, σ²=$Var")
# plot(r,smooth(grt_p,over=1),title=L"g_{\!\!+}"*" for T=$T, σ²=$Var")

vec(mean(grt_p,dims=2))

## Basic Tests
L = 200
    T = 0.1
    Var = 0.1
    init = "hightemp"
    init == "isolated" ? BC = "free" : BC = "periodic"
    q = +1
    r0 = Int(L/2)

thetas,omegas,dt = initialize(L,T,Var,init,q,r0)
# display_thetas(thetas)
t = 0.0 ; tmax = 1500
z = @elapsed while t < tmax
    t += dt
    update(thetas,omegas,L,T,dt)
end
p = display_thetas(thetas)
highlight_defects(display_thetas(thetas),thetas)

spot_defects(thetas,BC)

#= Tested :
    * 4 different initialisations
    * display_thetas (with and without title)
    * highlight_defects (OK, pas de doublon)
    * update (0.01 second / update (L=200) )

=#
