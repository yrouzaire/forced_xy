cd("D:/Documents/Research/projects/forced_XY")
    using JLD2,StatsBase,BenchmarkTools,Statistics,Distributions,LinearAlgebra,Parameters,Random,SpecialFunctions
    include("methods.jl")
    using Plots,ColorSchemes,LaTeXStrings
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5)
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&

## Brute Force SAW
T = 2000
p = 1/3
M = Int(1E8)
SD = [Float64[] for m in 1:T]
z = @elapsed for m in 1:M
    walk = SAW(T,p,"BF")
    for t in 1:length(walk)
        push!(SD[t],walk[t]...)
    end
end
prinz(z)
# jldsave("data/Gxt_SAW_BruteForce.jld2";SD,T)

# MSD = [mean(abs2,SD[t]) for t in 1:T]
# plot!(1:T,remove_negative(MSD),axis=:log)
#     plot!(x->0.3x^1.5)


@unpack SD,T = load("data/Gxt_SAW_BruteForce.jld2")
pgxtSAW=plot()
    for i in 10:5:30
        data = abs.(SD[i])
        N = length(data)
        h = fit(Histogram,data,0:floor(Int,maximum(data)))
        r = Array(h.edges[1])
        plot!(r[2:end-1],(remove_negative(h.weights[2:end]/2N)),rib=0,yaxis=:log)
    end
    annotate!(1,3,text(L"G(x,t)",15,:left))
    annotate!(20.2,2E-5,text(L"x",15,:center))
    xlims!(0,22)
    ylims!(8E-6,10)
    yticks!([1E-5,1E-4,1E-3,1E-2,1E-1,1,10])
    pgxtSAW
# savefig("figures/Gxt_SAW.svg")


p2=plot(xticks=0:2:6,size=(200,200))
    expo = 3
    for i in 10:5:30
        data = abs.(SD[i])
        N = length(data)
        h = fit(Histogram,data,0:floor(Int,maximum(data)))
        # h = normalize(h,mode=:density) # if bins not equal
        r = Array(h.edges[1])
        # plot!(r[2:end-1].^2/(i*every)^expo,(i*every)^(expo/2)*smooth(remove_negative(h.weights[2:end]/2N),1),rib=0,yaxis=:log)
        plot!(r[2:end-1].^4/(i)^expo,(i)^(expo/4)*smooth(remove_negative(h.weights[2:end]/N),1),rib=0,yaxis=:log)
        #weights/TWO N because here only one sided so we have the double of the actual distrib
    end
    xx = 0:0.1:6
    plot!(xx,0.7exp.(-0.85xx),c=:black,line=:dash)
    # annotate!(8,0.005,text(L"\sim f(y) = e^{-0.8\,y}",10,:left))
    annotate!(0,1.9,text(L"t^{3/4}\,G(x,t)",10,:left))
    annotate!(4.3,0.0035,text(L"x^4/t^3",10,:left))
    # ylabel!(L"t^{3/4}\,G(x,t)")
    # xlabel!(L"x^4/t^3")
    ylims!(0.002,4)
    xlims!(-0.3,6.5)
    p2
savefig("figures/Gxt_SAW_inset.svg")

# p2=plot(xlabel=L"x^4/t^3",ylabel=L"t^{3/4}\,G",yticks=[1E-5,1E-4,1E-3,1E-2,1E-1,1])
#     expo = 3*3/4
#     for i in 15:3:25
#         data = abs.(SD[i])
#         N = length(data)
#         h = fit(Histogram,data,0:floor(Int,maximum(data)))
#         # h = normalize(h,mode=:density) # if bins not equal
#         r = Array(h.edges[1])
#         # plot!(r[2:end-1].^2/(i*every)^expo,(i*every)^(expo/2)*smooth(remove_negative(h.weights[2:end]/2N),1),rib=0,yaxis=:log)
#         # plot!(r[2:end-1].^3/(i)^expo,(i)^(expo/4)*smooth(remove_negative(h.weights[2:end]/N),1),rib=0,yaxis=:log)
#         plot!(r[2:end-1].^3.5/(i)^2.5,i^(5/7)*smooth(remove_negative(h.weights[2:end]/N),1),rib=0,yaxis=:log)
#         #weights/TWO N because here only one sided so we have the double of the actual distrib
#     end
#     xx = 0:0.1:15
#     # plot!(xx,2exp.(-0.8xx),c=:black,line=:dash)
#     # annotate!(8,0.005,text(L"\sim f(y) = e^{-0.8\,y}",10,:left))
#     p2

plot(p1,p2,size=(800,400))
savefig("figures/Gxt_SAW.pdf")

## Simple Sampling SAWs and testing MSD
T = 2000
M = Int(1E5)
every = 10
Gx = Matrix{Int}(undef,Int(T/every),M)*NaN
Gy = Matrix{Int}(undef,Int(T/every),M)*NaN
W = Matrix{Float64}(undef,Int(T/every),M)*NaN

z = @elapsed for m in 1:M
    walk,weight = SAW_weighted(T,100)
    # for t in 1:round(Int,length(walk)/every,RoundDown)
    #     push!(Gx[t],walk[every*t][1])
    #     push!(Gy[t],walk[every*t][2])
    #     push!(W[t],weight[every*t])
    # end
    Gx[:,m] = [el[1] for el in walk[every:every:T]]
    Gy[:,m] = [el[2] for el in walk[every:every:T]]
    W[:,m] = weight[every:every:T]
end
prinz(z)

average_distance = zeros(Int(T/every))
    @time for t in 1:Int(T/every)
        Sum = sum(W[t,:])
        average_distancex = sum((Gx[t,:]).^2 .*W[t,:])/ Sum
        average_distancey = sum((Gy[t,:]).^2 .*W[t,:])/ Sum
        average_distance[t] = (average_distancey + average_distancex)/2
    end
plot(every:every:T,average_distance,axis=:log,m=true)
plot!(x->0.4x^1.5)


MSD = [mean(abs2,G[i]) for i in each(G)]
plot(every:every:T,MSD,axis=:log,m=true)
    plot!(x->0.4x^1.5)
## Stats of SAW
T = 2000
M = Int(1E6)
every = 100
G = [Int[] for i in 1:Int(T/every)]
@unpack G = load("data/SAW_Gxt.jld2")
plot(every:every:T,remove_negative([length(G[x]) for x in eachindex(G)]),yaxis=:log,m=true)
[length(G[x]) for x in eachindex(G)]

p = plot()
z = @elapsed for n in 1:30
    println("n = $n")
        # walk_long = SAW(T)
        # for t in Int.(1:length(walk_long)/every)
        #     push!(G[t],walk_long[every*t]...)
        # end
        for m in 1:round(Int,M/n)
            walk = SAW(T,n)
            for t in 1:round(Int,length(walk)/every,RoundDown)
                push!(G[t],walk[every*t]...)
            end
        end
        display(plot!(every:every:T,remove_negative([length(G[x]) for x in eachindex(G)]),yaxis=:log))
    end
    prinz(z)
    jldsave("data/SAW_mixed_Gxt.jld2";T,M,every,G)

p1=plot(xlabel="x",ylabel="G(x,t)")
    for i in each(G)
        data = abs.(G[i])
        N = length(data)
        h = fit(Histogram,data,0:floor(Int,maximum(data)))
        r = Array(h.edges[1])
        plot!(r[2:end-1],smooth(remove_negative(h.weights[2:end]/2N),1),rib=0,yaxis=:log)
        #weights/TWO N because here only one sided so we have the double of the actual distrib
    end
    # varr = 200
    # plot!(0:50,1/sqrt(4*π*varr)*exp.(-Array(0:50).^2/4/varr);c=:black,line=:dot)
    p1


p2=plot(xlabel=L"x^3/t^2",ylabel=L"t^{2/3}\,G\,\Gamma(4/3)",yticks=[1E-5,1E-4,1E-3,1E-2,1E-1,1])
    expo = 3
    for i in each(G)
        data = abs.(G[i])
        N = length(data)
        h = fit(Histogram,data,0:floor(Int,maximum(data)))
        # h = normalize(h,mode=:density) # if bins not equal
        r = Array(h.edges[1])
        # plot!(r[2:end-1].^2/(i*every)^expo,(i*every)^(expo/2)*smooth(remove_negative(h.weights[2:end]/2N),1),rib=0,yaxis=:log)
        plot!(r[2:end-1].^4/(i*every)^expo,(i*every)^(expo/3)*gamma(4/3)*smooth(remove_negative(h.weights[2:end]/N),1),rib=0,yaxis=:log)
        #weights/TWO N because here only one sided so we have the double of the actual distrib
    end
    xx = 0:0.1:15
    plot!(xx,2exp.(-0.8xx),c=:black,line=:dash)
    annotate!(8,0.005,text(L"\sim f(y) = e^{-0.8\,y}",10,:left))
    p2

plot(p1,p2,size=(800,400))
savefig("figures\\SAW_Gxt.pdf")
# p2=plot(xlabel=L"x^2/t^{3/2}",ylabel=L"t^{3/4}\,G")
#     expo = 2
#     for i in 5:5:80
#         data = abs.(G[i])
#         N = length(data)
#         h = fit(Histogram,data,0:floor(Int,maximum(data)))
#         # h = normalize(h,mode=:density) # if bins not equal
#         r = Array(h.edges[1])
#         # plot!(r[2:end-1].^2/(i*every)^expo,(i*every)^(expo/2)*smooth(remove_negative(h.weights[2:end]/2N),1),rib=0,yaxis=:log)
#         plot!(r[2:end-1].^3/(i*every)^expo,(i*every)^(expo/3)*smooth(remove_negative(h.weights[2:end]/2N),1),rib=0,yaxis=:log)
#         #weights/TWO N because here only one sided so we have the double of the actual distrib
#     end
#     p2
