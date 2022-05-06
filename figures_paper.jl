cd("D:/Documents/Research/projects/forced_XY")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Plots,ColorSchemes,LaTeXStrings
    include("methods.jl")
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    colors5 = [:purple,ColorSchemes.tab10[1],ColorSchemes.tab10[3],:darkgoldenrod1,ColorSchemes.tab10[4],ColorSchemes.tab10[5]]
    pathfig = "C:/Users/Ylann Rouzaire/Dropbox/Ylann/PhD/article_long/figures/"
## Snapshots of the systems (FIGURE 1)
L = 200
    init = "hightemp" ; BC = "periodic"
    q = +1 ; r0 = Int(L/2) # dummy

#
# T = 0.2 ; Var = 0.
#     t = 0.0 ; thetasXY,omegasXY,dt = initialize(L,T,Var,init,q,r0)
# z = @elapsed while t < 800
#         t += dt ; update(thetasXY,omegasXY,L,T,dt)
#     end
# pXY = heatmap(mod.(thetasXY,2π)',c=cols,xlims=(0,L),ylims=(0,L),colorbar=false) ; highlight_defects(pXY,thetasXY)
#
# thetasXY_moved = copy(thetasXY)
#     for i in 1:L
#         thetasXY_moved[i,:] = circshift(thetasXY_moved[i,:],20)
#     end
#     heatmap(mod.(thetasXY,2π)',c=cols,xlims=(0,L),ylims=(0,L),colorbar=false)
# heatmap(mod.(thetasXY_moved,2π)',c=cols,xlims=(0,L),ylims=(0,L),colorbar=false)
# pXY = heatmap(mod.(thetasXY_moved,2π)',c=cols,xlims=(0,L),ylims=(0,L),colorbar=false) ; highlight_defects(pXY,thetasXY_moved)

#
# T = 0.2 ; Var = 0.1
#     t = 0.0 ; thetasNKM,omegasNKM,dt = initialize(L,T,Var,init,q,r0)
# z = @elapsed while length(spot_defects(thetasNKM,BC)) > 8 # on pourrait aussi faire while nb vortices \ne 4
#         t += 25dt ; for i in 1:25 update(thetasNKM,omegasNKM,L,T,dt) end
#         println("n = ",length(spot_defects(thetasNKM,BC)))
#     end
# pNKM = heatmap(mod.(thetasNKM,2π)',c=cols,xlims=(0,L),ylims=(0,L),colorbar=false) ; highlight_defects(pNKM,thetasNKM)
#
# thetasNKM_moved = copy(thetasNKM)
#     for i in 1:L
#         thetasNKM_moved[i,:] = circshift(thetasNKM_moved[i,:],-20)
#     end
# pNKM_moved = heatmap(mod.(thetasNKM_moved,2π)',c=cols,xlims=(0,L),ylims=(0,L),colorbar=false) ; highlight_defects(pNKM_moved,thetasNKM_moved)
# pNKM = pNKM_moved

#
# JLD2.jldsave("data\\figures_snapshots_systems.jld2";L,init,BC,T,VarNKM=0.1,pXY,pNKM)
# JLD2.jldsave("data\\snapshots_systems.jld2";L,init,BC,thetasNKM,thetasXY_moved,T,VarNKM=0.1,pXY,pNKM)
# @unpack pXY,pNKM = JLD2.load("data\\snapshots_systems.jld2")
# @unpack thetasNKM,L,BC = JLD2.load("data\\snapshots_systems.jld2")
# pNKM = heatmap(mod.(thetasNKM,2π)',c=cols,xlims=(0,L),ylims=(0,L),size=(430,400),colorbar_title="θ",clims=(0,2π)) ; highlight_defects(pNKM,thetasNKM)
#
# plot(pXY,pNKM,size=(830,400))
#     savefig("C:/Users/Ylann Rouzaire/Dropbox/Ylann/PhD/article_long/figures/snapshots.pdf")
    # savefig("figures\\snapshots.svg")
# heatmap(zeros(200,200),clims=(0,2pi),c=cols,colorbar_title="θ")
# savefig("figures\\colorbar.svg")

## FIGURE 2 C(r,t) et collapse
@unpack L,Ts,Vars,tmax,init,C,BC,times,R = load("data/Crt_paper_Var0.jld2")
C_avg_Var0 = mean(C,dims=4)[:,:,:,1]
@unpack L,Ts,Vars,tmax,init,C,BC,times,R = load("data/Crt_paper_Var0.1.jld2")
C_avg_Var01 = mean(C,dims=5)[:,:,:,:,1]
p2a = plot(xlabel="r")
    for tt in 5:4:length(times)
        plot!(1:Int(L/2),remove_negative(C_avg_Var0[1,:,tt]),axis=:log,alpha=0.96^tt,c=:purple,lw=1.7)
        plot!(1:Int(L/2),remove_negative(C_avg_Var01[1,1,:,tt]),axis=:log,alpha=0.97^tt,c=colors5[3],lw=1.7)
    end
    ylims!(0.03,1.8)
    plot!(1:Int(L/2),exp(-1)*ones(length(1:Int(L/2))),c=:grey,line=:dot,lw=1)
    annotate!(1.6,1.3,text(L"C(r,t)",12))
    annotate!(1.6,0.3,text("1/e",10,:gray30))


p2b = plot(xlabel=L"r\,\sqrt{\ln t\,/\,t}")
    for tt in 12:2:length(times)
        plot!(collect(1:Int(L/2))*sqrt(log(times[tt])/times[tt]),collect(1:Int(L/2)).^(T/2pi).*remove_negative(C_avg_Var0[1,:,tt]),axis=:log,alpha=0.96^tt,c=:purple,lw=1.7)
        plot!(collect(1:Int(L/2))*sqrt(log(times[tt])/times[tt]),collect(1:Int(L/2)).^(T/2pi).*remove_negative(C_avg_Var01[1,1,:,tt]),axis=:log,alpha=0.96^tt,c=colors5[3],lw=1.7)
    end
    ylims!(0.03,1.8)
    # plot!(1:Int(L/2),exp(-1)*ones(length(1:Int(L/2))),c=:grey,line=:dot,lw=1)
    annotate!(0.3,1.3,text(L"r^\eta\,C(r\,\sqrt{\ln t\,/\,t},t)",12))
    # annotate!(1.6,0.3,text("1/e",10,:gray30))

plot(p2a,p2b,size=(800,400),background_color=:transparent)
# savefig(pathfig*"Crt.svg")

## FIGURE 3 \xi et collapse
@unpack init,tsave,tmax,Ts,Vars = JLD2.load("data/scalingXI_Var_L500_HighTemp.jld")

# load("data/XY_saturationXI_L200_HighTemp.jld") # violet
# load("data/scalingXI_Var_L100_HighTemp.jld")
# load("data/scalingXI_Var_L200_HighTemp.jld")

seuil = exp(-1)
variances = [0.0,0.02,0.05,0.1,0.2] ; ind_variances = [1,3,6,11,21]
variances = [0.0,0.02,0.05,0.1,0.15,0.2] ; ind_variances = [1,3,6,11,16,21]
Ls = [200,500]
ξ = zeros(length(Ls),length(variances),39)
for i in eachindex(Ls)
    println(i)
    L = Ls[i] ; Lover2 = round(Int,L/2,RoundDown) ; r_vector = sort(vcat([n for n in 1:Lover2],[n*sqrt(2) for n in 1:Lover2]))
    C = load("data/scalingXI_Var_L$(L)_HighTemp.jld","C")
    C_avg = mean(C,dims=5)[1,ind_variances,:,:,1]
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

colors = [:purple,:darkorange,ColorSchemes.tab10[1],ColorSchemes.tab10[3],:darkgoldenrod1,ColorSchemes.tab10[4]]
    p3a = plot(xlabel="t",legend=(0.452,0.12))
    plot!([NaN],[NaN],c=:grey,m=:star,label="L = 200",line=false)
    plot!([NaN],[NaN],c=:grey,m=:utriangle,label="L = 500",line=false)
    for j in eachindex(variances)
        # plot!([NaN],[NaN],c=colors[j],rib=0,label="σ² = $(variances[j])")
        plot!(tsave[1:39],ξ[1,j,:],m=:star,axis=:log,c=colors[j];line=false)
        plot!(tsave[1:39],ξ[2,j,:],m=:utriangle,axis=:log,c=colors[j];line=false)
    end
    annotate!(110,17,text(L"\sim \sqrt{t\,/\,\ln t}",40.0,10))
    annotate!(1.7,21,text(L"ξ(t)",12))
    ylims!(0.8,1.2maximum(ξ))
    xlims!(0.8,1000)
    yticks!(vcat(1:10,20),string.([1,2,3,4,5,"","","","",10,20]))
    xticks!([1,10,100,1000],[L"10^0",L"10^1",L"10^2",L"10^3"])
    plot!(tsave[20:38],2.7 .* sqrt.(tsave[20:38] ./ log.(tsave[20:38])),c=:black)
    p3a
# savefig("figures/dynamic_coarsening_rest_legend.svg")

# INSET
ind_vars = [1,3,6,11,16,21]
p3b = plot(xlabel=L"σ²\,t")
    # plot!([NaN],[NaN],c=:grey,m=:star,label="L = 200",line=false)
    # plot!([NaN],[NaN],c=:grey,m=:utriangle,label="L = 500",line=false)
    for j in 2:length(variances)
        plot!(tsave[1:39]*variances[j],ξ[2,j,:]*variances[j]^0.5,axis=:log,m=:utriangle,c=colors[j];line=false)
        # plot!(tsave[1:39]exp(-8lambertw(-0.5variances[j]^0.5)),ξ[2,j,:]*variances[j]^0.5,axis=:log,m=:utriangle,c=colors[j];line=false)
        # plot!(tsave[1:39]*variances[j],ξ[2,j,:]*variances[j]^0.5,axis=:log,m=:utriangle,c=colors[j];line=false)
    end
    a = 2.3 ; b = 0.65  # on veut ab = 1.2 et a = 2.3
    plot!(0.01:0.01:210,x->a*(1-exp(-b*(x)^.5)),c=:black,lw=1.7)
    annotate!(0.02,2.3,text(L"σ\,ξ",12))
    # annotate!(70,.15,text(L"σ²t",12))
    ylims!(0.1,3)
    yticks!(vcat(.1:0.1:1,2,3),string.(Any[0.1,"","","","","","","","",1,"",""]))
    xticks!([.01,0.1,1,10,100],[L"10^{-2}","",L"10^0","",L"10^2"])
    p3b

plot(p3a,p3b,size=(800,400),background_color=:transparent)
savefig(pathfig*"xi.svg")
&

# Inset Time to reach SS versus sigma
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
        end
    end
    Varss = 0.01:0.01:0.25
    esp = zeros(length(Varss))
    for i in eachindex(Varss)
        σ = sqrt(Varss[i])
        g(x) = cos(x)*x*exp(-(sin(x)/σ)^2)*2/σ/sqrt(pi)
        esp[i] = quadgk(g,0,pi/2)[1]
    end

    # calcul du temps nécessaire pour que le cas XY atteigne la valeur \ell = pi/2/esp
    tau = Vector(undef,length(Vars))
    for i in 1:length(Vars)
        ind_tmp = findfirst(x->x≥π/2/esp[i],ξ[1,:])
        if ind_tmp ≠ nothing tau[i] = tsave[ind_tmp] end
    end

    p2c = plot(xlabel="σ²",size=(200,200))
    plot!(Vars[3:end],tau[3:end],axis=:log,m=:circle,ms=5,line=false,c=:grey)
    scatter!((Vars[3],tau[3]),axis=:log,m=:circle,ms=5,c=colors[2])
    scatter!((Vars[6],tau[6]),axis=:log,m=:circle,ms=5,c=colors[3])
    scatter!((Vars[11],tau[11]),axis=:log,m=:circle,ms=5,c=colors[4])
    scatter!((Vars[16],tau[16]),axis=:log,m=:circle,ms=5,c=colors[5])
    scatter!((Vars[21],tau[21]),axis=:log,m=:circle,ms=5,c=colors[6])
    plot!(Varss[2:end],6 ./Varss[2:end],c=:black)
    xticks!(vcat(0.01:0.01:0.1,0.2:0.1:0.4),string.([L"10^{-2}","","","","","","","","",L"10^{-1}","","",""]))
    yticks!(vcat(10:10:100,200:100:1000),string.([L"10^{1}","","","","","","","","",L"10^{2}","","","","","","","","",L"10^{3}"]))
    xlims!(0.01,0.3)
    ylims!(10,1500)
    annotate!(0.016,800,text(L"τ",12))
    annotate!(0.123,200,text(L"∼1/σ²",12))

# savefig(pathfig*"scaling_transients.svg")

## FIGURE 4 number of vortices defects
@unpack L,init,BC,tmax,Ts,Vars,times,nb_vortices = load("data/phase_space_number_vortices_hightemp.jld")
nb_vortices_avg = mean(nb_vortices,dims=4)[:,:,:,1]
colors = [:purple,ColorSchemes.tab10[1],ColorSchemes.tab10[3],:darkgoldenrod1,ColorSchemes.tab10[4],ColorSchemes.tab10[5]]
p4a=plot(xlabel="t",legend=false)
    token = 1
    plot!(times[2:end],nb_vortices_avg[6,1,2:end],axis=:log,linewidth=2,c=:purple,label="σ² = 0")
    plot!(times[2:end],NaN*nb_vortices_avg[6,1,2:end],axis=:log,line=false,ms=2.5,m=:circle,c=:darkorange,label="σ² = 0.02")
    for j in ([6,11,16,21,26])
        plot!(times[2:end],nb_vortices_avg[9,j,2:end],axis=:log,line=false,m=:circle,ms=2.5,c=colors[token+1],label="σ² = $(Vars[j])")
        token += 1
    end
    p2d
    ylims!(1,20000)
    xticks!([1,10,100,1000],[L"10^0",L"10^1",L"10^2",L"10^3"])
    yticks!([1,10,100,1000,10000],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"])
    plot!(times[17:end],1E3 * log.(times[17:end]) ./times[17:end],c=:black)
    annotate!(370,8.5,text(L"\sim \ln t\,/\,t",-40.0,10))
    annotate!(2,10000,text(L"n(t)",15))


@unpack L,init,BC,tmax,C,ξ,nb_defects,R,times,T,Vars = load("data\\xi_sqrtn_holds.jld2")

colors = [:purple,:darkorange,ColorSchemes.tab10[1],ColorSchemes.tab10[3],ColorSchemes.tab10[4],ColorSchemes.tab10[4]]
p4b=plot(legend=false,xlabel="t")
    for i in [1,2,3]
        plot!(times,smooth(mean(ξ,dims=3)[:,i,1].* sqrt.(mean(nb_defects,dims=3)[:,i,1])/L,over=1),
        xaxis=:log,label="σ² = $(Vars[i])",c=colors[i],rib=0)
    end
    plot!(times[1:end-5],smooth(mean(ξ,dims=3)[1:end-5,4,1].* sqrt.(mean(nb_defects,dims=3)[1:end-5,4,1])/L,over=1),
    xaxis=:log,label="σ² = $(Vars[4])",c=colors[4],rib=0)
    plot!(times,smooth(mean(ξ,dims=3)[:,6,1].* sqrt.(mean(nb_defects,dims=3)[:,6,1])/L,over=1),
    xaxis=:log,label="σ² = $(Vars[6])",c=colors[6],rib=0)
    ylims!(0.28,0.58)
    xlims!(0.9,300)
    annotate!(4,0.554,text(L"ξ(t)⋅\sqrt{n(t)}/L",12))

# plot(p4a,p4b,size=(800,400),background_color=:transparent)
# savefig(pathfig*"n.svg")


## Nombre total de vortex : phase space and over time (FIGURE 1b et FIGURE 6a)
@unpack L,init,BC,tmax,Ts,Vars,times,nb_vortices = load("data/phase_space_number_vortices_lowtemp.jld")
nb_vortices_avg = mean(nb_vortices,dims=4)[:,:,:,1]

# Final time
@unpack Ts,Vars,nb_vortices = load("data/phase_space_number_vortices_hightemp.jld")
nb_vortices_complement1 = load("data/phase_space_number_vortices_complement_hightemp.jld","nb_vortices")
nb_vortices_complement2 = load("data/phase_space_number_vortices_complement2_hightemp.jld","nb_vortices")
tmax = 15000
nb_final_time = nb_vortices[:,:,end,:]
nb_final_time[1:5,:,:] = nb_vortices_complement1[:,:,end,:]
nb_final_time[6:25,:,:] = nb_vortices_complement2[:,:,end,:]
nb_vortices_avg = mean(nb_final_time,dims=3)[:,:,1]
p1 = heatmap(Vars,Ts,log10.(1 .+ nb_vortices_avg),c=cgrad([:blue,:green,:orange,:red]),clims=(0,log10.(1 .+ maximum(nb_vortices_avg[:,:,end]))),
    xlabel="σ²",ylabel="T",colorbar_title=L"\log_{10} (1+n)",size=(500,400))
    # scatter!((0.03,0.1),m=:dtriangle,c=:white,ms=7)
    # scatter!((0.1,0.3),m=:circle,c=:white,ms=6)
    # scatter!((0.03,0.3),m=:star5,c=:white,ms=6)
    # scatter!((0.1,0.1),m=:utriangle,c=:white,ms=5)
    xlims!(0,0.25)
    ylims!(0,0.8)
# savefig("figures/figures_paper/n_phasespace.pdf")

# savefig(pathfig*"n_phasespace.pdf")

# n(t)
# @unpack L,init,BC,tmax,Ts,Vars,times,nb_vortices = load("data/phase_space_number_vortices_hightemp.jld")
# nb_vortices_avg = mean(nb_vortices,dims=4)[:,:,:,1]
# colors = [:purple,ColorSchemes.tab10[1],ColorSchemes.tab10[3],:darkgoldenrod1,ColorSchemes.tab10[4],ColorSchemes.tab10[5]]
# p2d=plot(xlabel="t",legend=false)
#     token = 1
#     plot!(times[2:end],nb_vortices_avg[6,1,2:end],axis=:log,linewidth=2,c=:purple,label="σ² = 0")
#     plot!(times[2:end],NaN*nb_vortices_avg[6,1,2:end],axis=:log,line=false,ms=2.5,m=:circle,c=:darkorange,label="σ² = 0.02")
#     for j in ([6,11,16,21,26])
#         plot!(times[2:end],nb_vortices_avg[9,j,2:end],axis=:log,line=false,m=:circle,ms=2.5,c=colors[token+1],label="σ² = $(Vars[j])")
#         token += 1
#     end
#     p2d
#     ylims!(1,20000)
#     xticks!([1,10,100,1000],[L"10^0",L"10^1",L"10^2",L"10^3"])
#     yticks!([1,10,100,1000,10000],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"])
#     plot!(times[17:end],1E3 * log.(times[17:end]) ./times[17:end],c=:black)
#     annotate!(370,8.5,text(L"\sim \ln t\,/\,t",-40.0,10))
#     annotate!(2,10000,text(L"n(t)",15))
#     # lens!([10,500],[1.5,10], inset = (1, bbox(0.55, 0.02, 0.45, 0.4)))

    # savefig("figures/decay_number_defects.pdf")

## Legend FIGURE 2
# p=plot(legend=true)
#     for i in 1:length(colors5)
#         plot!(rand(10),rib=0,c=colors5[i],label="zefr")
#     end
#     plot!(rand(10),rib=0,c=:darkorange,label="zefr")
#     p
# savefig("figures\\figures_paper\\colors_legend_fig2.svg")


## Impact Thermal Noise on the system FIGURE SM
# T = 0.2  ; Var = 0.2 ; init = "lowtemp"
# # At short times, it does not change much the domains
# tmax = 200
# @unpack plots = load("data/impact_noise_$(init)_fixed_omegas_fixed_T$(T)_Var$(Var)_tmax$(tmax).jld2")
# for plot in plots
# plot(p...,size=(800,800))
# # savefig("figures/figures_paper/impact_noise_$(init)_fixed_omegas_fixed_T$(T)_Var$(Var)_tmax$(tmax).png")
#
# # At long times, it does change a lot, meaning that the defects' motion is affected by thermal noise
# tmax = 1000
# @unpack p = load("data/impact_noise_$(init)_fixed_omegas_fixed_T$(T)_Var$(Var)_tmax$(tmax).jld2")
# plot(p...,size=(800,800))
# savefig("figures/figures_paper/impact_noise_$(init)_fixed_omegas_fixed_T$(T)_Var$(Var)_tmax$(tmax).png")


## Sample_trajectory_time_theta_thetadot FIGURE SM
@unpack L,thetas_history,thetas_dot_history,tt_interesting,TV,init,BC,tmax,pAtheta,pAthetadot,pBtheta,pBthetadot,pCtheta,pCthetadot,ptheta,pthetadot = load("data/snapshots_theta_thetadot.jld2")

function visualisation_theta_thetadot_ij(i,j,thetas_history,thetas_dot_history,dt,legendparam=false)
    tt = collect(1:length(thetas_history[i,j,1:end]))*dt
    p1 = plot(xlabel="t",ylabel="θ(t)",legend=legendparam)
    plot!(tt,thetas_history[i,j,1:end],rib=0,label="Original Signal")
    plot!(tt,thetas_history_ema[i,j,1:end],rib=0,label="Smoothed Signal")
    # ylims!(-0.5,1.5)
    p2 = plot(xlabel="t",ylabel=L"d\theta/dt")
    plot!(tt,thetas_dot_history[i,j,1:end],c=:black)
    plot!(tt,0*tt,c=:grey,lw=1)
    # ylims!(-0.008,0.004)
    return plot(p1,p2,layout=(2,1),size=(400,800),grid=true)
end

using MarketTechnicals
t = 0:0.01:10
x = sin.(t)+ 0.5randn(length(t))
y1 = ema(x,50,wilder=true)
y2 = ema_maison(x,50)
    plot(t,x)
    plot!(t[1:end-49],y1)
    plot!(t,1 .+circshift(y2,-50))

function ema_maison(x,window)
    y = zeros(length(x))
    for t in 1:length(y)
        y[t] = sum([x[t-i]*((window-1)/window)^(i) for i in 0:t-1])/window
    end
    return y
end

# Visualisation temporal series
# (61, 60) , (49, 78) , (8, 76) , (20, 97) , (28, 45)
pa = visualisation_theta_thetadot_ij(28,45,thetas_history,thetas_dot_history,determine_dt(0,0),:topleft)
pb = visualisation_theta_thetadot_ij(20,97,thetas_history,thetas_dot_history,determine_dt(0,0))
pc = visualisation_theta_thetadot_ij(8,76,thetas_history,thetas_dot_history,determine_dt(0,0))
plot(pa,pb,pc,layout=(1,3),size=(1200,800))
savefig("figures\\figures_paper\\sample_trajectory_time_theta_thetadot.pdf")

# plot(thetas_dot_history[i,j,1:end])
#     plot!(ema(thetas_dot_history[i,j,1:end],50,wilder=true))

# tt_interesting = [45,447,568]
# tt = 45 ;
# pAtheta = heatmap(mod.(thetas_history[:,:,round(Int,tt/determine_dt(0,0))]',2pi),c=cols,clims=(0,2π),colorbar=false,ticks=false)
# pAthetadot = heatmap(thetas_dot_history[:,:,round(Int,tt/determine_dt(0,0))]',colorbar=false,clims=(-0.02,0.02),ticks=false)
# tt = 447 ;
# pBtheta = heatmap(mod.(thetas_history[:,:,round(Int,tt/determine_dt(0,0))]',2pi),c=cols,clims=(0,2π),colorbar=false,ticks=false)
# pBthetadot = heatmap(thetas_dot_history[:,:,round(Int,tt/determine_dt(0,0))]',colorbar=false,clims=(-0.02,0.02),ticks=false)
# tt = 568 ;
# pCtheta = heatmap(mod.(thetas_history[:,:,round(Int,tt/determine_dt(0,0))]',2pi),c=cols,clims=(0,2π),colorbar=true,ticks=false,colorbar_title="θ",size=(450,400))
# highlight_defects(pCtheta,thetas_history[:,:,round(Int,tt/determine_dt(0,0))])
# pCthetadot = heatmap(thetas_dot_history[:,:,round(Int,tt/determine_dt(0,0))]',colorbar=true,clims=(-0.02,0.02),ticks=false,colorbar_title=L"d\theta/dt",size=(480,400))

# ptheta = plot(pAtheta,pBtheta,pCtheta,layout=(1,3),size=(1250,400))
# # savefig("figures\\figures_paper\\snapshots_theta.png")
# # savefig("figures\\figures_paper\\snapshots_theta.svg")
# pthetadot = plot(pAthetadot,pBthetadot,pCthetadot,layout=(1,3),size=(1250,400))
# # savefig("figures\\figures_paper\\snapshots_thetadot.png")
# # savefig("figures\\figures_paper\\snapshots_thetadot.svg")

## Films nucléation thetadot->theta FILM SM
filename = "data/films_T0.05_lowtemp.jld"
TV,tmax,init,BC,L = JLD2.load(filename,"metadata").values

i = 4 ; dt = determine_dt(T,Var,L) ; T,Var = TV[i]
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
    end
    println("Computation time : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes.")


fin = findfirst(isnan,thetas_dot_history[1,1,:])-1-tau
    if fin == nothing fin = length(thetas_dot_history) end
    zm = @elapsed anim = @animate for i in 10:10:fin
        p1 = heatmap(mod.(thetas_history[:,:,i+tau]',2π),c=cgrad([:black,:blue,:green,:orange,:red,:black]),clims=(0,2π),size=(600,512),colorbar_title="θ")
        highlight_defects(p1,(mod.(thetas_history[:,:,i],2π)))
        p2 = heatmap(thetas_dot_history[:,:,i]',size=(600,512),clims=(-0.05,0.05),colorbar_title="dθ/dt")
        p = plot(p1,p2,size=(1200,512))
        title!("t=$(round(dt*i,digits=2))")
    end
    println()
    println("Time Movie : $(round(Int,zm)) seconds = $(round(zm/60,digits=2)) minutes.")
    # mp4(anim,"figures/figures_paper/film_test.mp4",fps=10);
    mp4(anim,"figures/figures_paper/film_$(init)_T$(T)_Var$(Var)_fixedscale.mp4",fps=10);

## Lieux de création et annihilation comparaison XY et NKM : FIGURE 4
L = 1000
    T = 0.5
    Var = 0
    init = "lowtemp"
    BC = "periodic"
    q = +1 ; r0 = 50 # dummy
    nb_creation_wanted = 1000
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
while length(fiche_vortices) < nb_creation_wanted && t < tmax
    println("Already $(length(fiche_vortices)) vortex creations.")
    t += dt ; update(thetas,omegas,L,T,dt)
    vortices_new,antivortices_new = spot_defects_separated(thetas,BC)
    fiche_vortices,fiche_antivortices = defect_tracker(thetas,BC,fiche_vortices,fiche_antivortices,vortices_old,antivortices_old,t)
    vortices_old,antivortices_old = vortices_new,antivortices_new
end


# function skinnycross(w = 0.1)
#     x = [w,w,1]; rx = reverse(x)
#     y = [1,w,w]; ry = reverse(y)
#     Shape(vcat(x,rx,-x,-rx), vcat(y,-ry,-y,ry))
# end
# skinnyx(w=0.1) = rotate(skinnycross(w), 0.25π)

dm,dv = density_creation(fiche_vortices,L,20)
heatmap(dm)
histogram(dv)
plot(xlims=(0,L),ylims=(0,L))
    scatter!(location_creation(fiche_vortices),m=:circle,c=:limegreen,ms=6)
    scatter!(location_annihiliation(fiche_vortices),m=skinnyx(0.2),c=:red)
    # display_thetas(thetas)
    # lens!([145, 160], [180,188], inset = (1, bbox(0.065, 0.05, 0.2, 0.2)))
# savefig("figures\\figures_paper\\locations_creations_NKM_T$(T)_Var$(Var).svg")
# savefig("figures\\figures_paper\\locations_creations_XY_T$(T).svg")

## Snaps to show the progression of the defects on the boundaries FIGURE 5
@unpack L,T,Var,init,BC,transients,tmax,times,save_thetas = load("data\\defects_surf.jld2")

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
# savefig("figures\\figures_paper\\defects_surf.pdf")


## FIGURE 6a Vortices are Free
using LsqFit,SpecialFunctions,LambertW,Roots
L,R,r0s,TV,tmax,save_every1,locations_pair_vortices,distance_pair_vortices = JLD.load("archives\\data_end_pdm\\data_vortices/pair_defaults0.jld","L","R","r0s","TV","tmax","save_every","locations","distance")
L,R,r0s,TV,tmax2,save_every2,locations,distance = JLD.load("../vortices_study/data_vortices/pair_defaults_PBC_L200_T0.1_Var0.1.jld","L","R","r0s","TV","tmax","save_every","locations","distance")
dt = determine_dt(0,0)
    distance_pair_vortices_avg  = Array{Float64}(undef,(length(0:save_every1:tmax/dt),length(TV),length(r0s))) ; distance_pair_vortices_std  = Array{Float64}(undef,(length(0:save_every1:tmax/dt),length(TV),length(r0s)))
    for n in 1:size(distance_pair_vortices)[3]
        for i in 1:size(distance_pair_vortices)[2]
            for t in 1:size(distance_pair_vortices)[1]
                distance_pair_vortices_avg[t,i,n] = mean(skipmissing(distance_pair_vortices[t,i,n,:]))
                distance_pair_vortices_std[t,i,n] = std(skipmissing(distance_pair_vortices[t,i,n,:]))
            end
        end
    end
    distance_avg  = Array{Float64}(undef,(length(0:save_every2:tmax2/dt),length(TV),nR0)) ; distance_std  = Array{Float64}(undef,(length(0:save_every2:tmax2/dt),length(TV),nR0))
    for n in 1:size(distance)[3]
        for i in 1:size(distance)[2]
            for t in 1:size(distance)[1]
                distance_avg[t,i,n] = mean(skipmissing(distance[t,i,n,:]))
                distance_std[t,i,n] = std(skipmissing(distance[t,i,n,:]))
            end
        end
    end

tv = 3
    tt1 = Array(0:dt*save_every1:tmax)
    p0 = [.01]
    # lastindex = [52,204,300]
    # plotlastindex = [53,208,413]
    lastindex = [60,223,340]
    plotlastindex = [60,223,391]
    p6a  = plot(xlabel=L"t",legend=nothing,size=(400,400))
    for n in 1:length(r0s)
        r0 = r0s[n].-2
        plot!(tt1,distance_pair_vortices_avg[1:637,1,n],rib=0,c=n)
        @. model(x, p) = r0*exp.(- (erfinv.(p[1]*x/r0)).^2)
        ind = lastindex[n]
        fitt = curve_fit(model, tt1[1:ind], distance_pair_vortices_avg[1:ind,1,n], p0)
        ind = plotlastindex[n]
        println(coef(fitt)[1],"   ",confidence_interval(fitt, 0.05)[1])
        plot!(tt1[1:ind],model(tt1[1:ind],coef(fitt)),line=:dash,c=:black)
    end

# L,R,r0s,TV,tmax2,save_every2,locations,distance = JLD.load("../vortices_study/data_vortices/pair_defaults_PBC_L200_T0.1_Var0.1.jld","L","R","r0s","TV","tmax","save_every","locations","distance")
# nR0 = length(r0s)

    tt2 = Array(0:dt*2save_every2:tmax2)
    plot!(tt2,distance_avg[1:2:end,1,1],line=:dot,c=1)
    plot!(tt2,distance_avg[1:2:end,1,2],line=:dot,c=2)
    plot!(tt2,distance_avg[1:2:end,1,3],line=:dot,c=3)
    xlims!((-50,1200))
    ylims!((-1,38))
    annotate!(-30,35.5,text(L"⟨R(t)⟩",15,:left,:black))
    annotate!(-30,38*4/100,text("(a)",15,:left,:black))

savefig("figures\\figures_paper\\fig6a.svg")
## FIGURE 6b MSD
# XY Model T = 0:0.025:0.7 et Var = 0
L,R,save_every,tmax,dt,MSD_avg,MSD_std,Ts = JLD2.load("data\\MSD_single_default_L200_Var0.0_R480.jld","L","R","save_every","tmax","dt","MSD_avg","MSD_std","Ts")
    Vars = [0.0]
    p6binset =plot(xlabel=L"T\,t\,/\,\ln t",axis=:log,size=(200,200),legend=false)
    couleur = 1
    styles = [:solid,:dashdot,:dash]
    t = Array(0:save_every:tmax/dt)*dt
    tlogt = t ./ log.(t)
    for j in 1:length(Vars)
        for i in [9,13,17,21]
            plot!(tlogt[3:end]*Ts[i],remove_negative(MSD_avg[3:end,i,j],1/R),rib=0,c=couleur,label="T = $(Ts[i])")
            global couleur += 1
        end
    end
    # plot!([2,150],1.075*[2,150],line=:dash,c=:black)
    plot!([4,150],2.575*[4,150],line=:dash,c=:black)
    # plot!([5,1000],0.3*[5,1000],line=:dash,c=:black)
    annotate!(4,3,text(L"\sim t\,/\,\ln t",12,:black,:left))
    # annotate!(10,350,text(L"\sim t\,/\,\ln t",12,:black,:left))
    # annotate!(0.49,800,text(L"⟨(r(t)-r_0)^2⟩",15,:black,:left))
    annotate!(0.69,400,text(L"MSD",15,:black,:left))
    ylims!(0.25,1500)
    # xlims!(0.08,200)
    # annotate!(0.09,0.48,text("(c)",15,:left))
# savefig("figures\\figures_paper\\inset_MSD_XY.svg")


# Forced Model T = 0.05 et Var = 0
L,R,save_every,tmax,dt,MSD_avg,MSD_std,Ts,Vars = load("data\\MSD_single_default_L200_T0.05_R396.jld","L","R","save_every","tmax","dt","MSD_avg","MSD_std","Ts","Vars")
p6b = plot(xlabel=L"\sigma t",axis=:log,size=(400,400),legend=:bottomright)
    couleur = 1
    styles = [:solid,:dashdot,:dash]
    t = Array(0:save_every:tmax/dt)*dt
    tlogt = t ./ log.(t)
    for j in 2:2:length(Vars)
        for i in 1:length(Ts)
            plot!(t[2:500]*Vars[j]^0.5,remove_negative(MSD_avg[2:500,i,j],1/R),c=couleur,rib=0,label="σ² = $(Vars[j])")
            global couleur += 1
        end
    end
    plot!([10,200],0.35*[10,200].^1.5,line=:dash,c=:black)
    annotate!(30,500,text(L"\sim  (σt)^{3/2}",12,:black,:center))
    annotate!(0.1,1700,text(L"MSD",15,:black,:left))
    # annotate!(0.1,1700,text(L"⟨(r(t)-r_0)^2⟩",15,:black,:left))
    ylims!(0.3,3000)
    xlims!(0.08,500)
    xticks!([0.1,1,10,100],[L"10^{-1}",L"10^{0}",L"10^{1}",L"10^{2}"])
    yticks!([1,10,100,1000],[L"10^{0}",L"10^{1}",L"10^{2}",L"10^{3}"])

# savefig("figures\\figures_paper\\MSD_forced.svg")


## FIG UNDERDAMPED
# xi(t) Underdamped compared to overdamped
styles = [:solid,:dash,:dot,:dash,:solid]
@unpack L,TV,ms,init,BC,tmax,times,ξ,C,n = load("data/comparaison_xin_T.jld")
    R = size(C,5)
    C_avg = mean(C,dims=5)[:,:,:,:,1]
    ξ_avg = mean(ξ,dims=4)[:,:,:,1]
    n_avg = mean(n,dims=4)[:,:,:,1]


# ξ(t) XY
pxixy = plot(legend=false)
    for i in 5
        plot!(times,ξ_avg[i,1,:],c=1,axis=:log,label="m=0",rib=0)
        for j in 2:length(ms)
            plot!(times*sqrt(ms[j]),ξ_avg[i,j,:],c=j,axis=:log,label="m=$(ms[j])",rib=0)
        end
    end
    plot!(times[10:end-5],3.5sqrt.(times[10:end-5] ./ log.(10times[10:end-5])),line=:dashdot,label="t/log t",c=:black)
    annotate!(5000,1.2,text(L"t",15,:black,:left))
    annotate!(.15,30,text(L"ξ(t)",15,:black,:left))
    annotate!(.5,15,text(L"\sim\, \sqrt{t \,/ \,\ln t\ }",12,:black,:left))
    annotate!(5000,32,text("(a)",15,:black,:center))
    xlims!(.11,1.2E4)
    pxixy
# savefig("figures\\comparaison_under_over.pdf")

# ξ(t) Forced
pxi = plot(legend=false)
    for i in [1,4]
        plot!(times,ξ_avg[i,1,:],c=1,axis=:log,ls=styles[i])
        # plot!(times,NaN*times,c=:grey,axis=:log,ls=styles[i],label="T=$(TV[i][1]) , σ²=$(TV[i][2])")
        for j in 2:4#length(ms)
            plot!(times*sqrt(ms[j]),ξ_avg[i,j,:],c=j,axis=:log,ls=styles[i])
        end
        plot!(times[1:end-3]*sqrt(ms[5]),ξ_avg[i,5,1:end-3],c=5,axis=:log,ls=styles[i])
    end
    plot!(times[10:end-5],3.5sqrt.(times[10:end-5] ./ log.(10times[10:end-5])),line=:dashdot,c=:black)
    # scatter!((7800,17.8),m=:dtriangle,c=:black,ms=6)
    # scatter!((7800,14.7),m=:star5,c=:black,ms=7)
    # scatter!((7800,9.),m=:utriangle,c=:black,ms=6)
    # scatter!((7800,7.7),m=:circle,c=:black,ms=5)
    annotate!(5000,1.2,text(L"t",15,:black,:left))
    annotate!(.35,23,text(L"ξ(t)",15,:black,:left))
    annotate!(4500,24,text("(b)",15,:black,:center))
    # annotate!(.5,15,text(L"\sim\, \sqrt{t \,/ \,\log t\ }",12,:black,:left))
    pxi
    ylims!(1,30)
    xlims!(.25,1E4)


# n(t) XY
pnxy = plot(legend=false)
    for i in 5
        plot!(times,(1 .+ n_avg[i,1,:]),c=1,axis=:log,label="m=0",rib=0)
        for j in 2:length(ms)
            plot!(times*sqrt(ms[j]),(1 .+ n_avg[i,j,:]),c=j,axis=:log,label="m=$(ms[j])",rib=0)
        end
    end
    plot!(times[10:end],1.1E3*(log.(10times[10:end]) ./ times[10:end]),line=:dashdot,label="t/log t",c=:black)
    annotate!(2000,1.6,text(L"t",15,:black,:left))
    annotate!(.013,22000,text(L"n(t)",15,:black,:left))
    annotate!(2,15,text(L"\sim\, \ln t\,/\, t",12,:black,:left))
    xlims!(.01,1E4)
    xticks!([0.01,0.1,1,10,100,1000,1E4],[L"10^{-2}","",L"10^{0}","",L"10^{2}","",L"10^{4}"])
    ylims!(1,4E4)
    annotate!(3500,19000,text("(c)",15,:black,:center))
    pnxy
# savefig("figures\\comparaison_under_over.pdf")

# n(t) Forced
pn = plot(legend=:bottomleft)
    plot!([NaN,NaN],ls = :solid, c=:grey,label="T=0.1, σ²=0.03")
    plot!([NaN,NaN],ls = :dash, c=:grey,label="T=0.3, σ²=0.1")
    for i in [1,4]
        plot!(times,(1 .+ n_avg[i,1,:]),c=1,axis=:log,ls=styles[i])
        for j in 2:5 #length(ms)
            plot!(times*sqrt(ms[j]),(1 .+ n_avg[i,j,:]),c=j,axis=:log,ls=styles[i])
        end
    end
    plot!(times[10:end],1.1E3*(log.(10times[10:end]) ./ times[10:end]),line=:dashdot,c=:black)
    annotate!(2000,1.6,text(L"t",15,:black,:left))
    annotate!(.013,22000,text(L"n(t)",15,:black,:left))
    ylims!(1,4E4)
    xticks!([0.01,0.1,1,10,100,1000,1E4],[L"10^{-2}","",L"10^{0}","",L"10^{2}","",L"10^{4}"])
    annotate!(7000,19000,text("(d)",15,:black,:center))
    pn
    # ylims!(1,30)
    # xlims!(.25,1E4)

plot(pxixy,pxi,pnxy,pn,layout=(2,2),size=(800,800))
# savefig(pathfig*"xin_comparaison.png")


## MSD Underdamped
@unpack L,TV,ms,init,BC,tmax,times,history_locations = load("data/locations_MSD_different_m.jld")
transients = 300
R = size(history_locations,4)
MSD = Array{Union{Missing,Float64}, 4}(missing,length(TV),length(times),length(ms),R)
    for r in 1:R , j in eachindex(ms) , i in eachindex(times) , ii in eachindex(TV)
        MSD[ii,i,j,r] = dist(history_locations[ii,i,j,r],history_locations[ii,1,j,r],L).^2
    end
MSD_avg = Array{Union{Missing,Float64}, 3}(missing,length(TV),length(times),length(ms))
    for j in eachindex(ms) , i in eachindex(times) , ii in eachindex(TV)
        MSD_avg[ii,i,j] = mean(skipmissing(MSD[ii,i,j,:]))
    end

# XY
pXY = plot(legend=false)
    plot!((times[2:end].-transients),2.5MSD_avg[2,2:end,3];c=1,axis=:log,label="m=0",rib=0)
    for j in 2:length(ms)
        plot!((times[2:end].-transients)*sqrt(ms[j]),MSD_avg[2,2:end,j];c=j,axis=:log,label="m=$(ms[j])",rib=0)
    end
    plot!((times[10:end].-transients),0.4 * (times[10:end].-transients)./ log.(times[10:end].-transients);c=:black,line=:dash)
    annotate!(6500,0.15,text(L"t",15,:black,:left))
    annotate!(1.1,140,text(L"MSD",15,:black,:left))
    annotate!(20,50,text(L"\sim\, t \,/ \,\log t",12,:black,:left))
    pXY

# KM
pKM = plot()
    plot!((times[2:end].-transients),4MSD_avg[1,2:end,3];c=1,axis=:log)
    for j in 2:length(ms)
        plot!((times[2:end].-transients)*sqrt(ms[j]),MSD_avg[1,2:end,j];c=j,axis=:log)
    end
    plot!((times[10:end].-transients),0.04 * (times[10:end].-transients).^1.5;c=:black,line=:dash)
    annotate!(6500,0.3,text(L"t",15,:black,:left))
    annotate!(1,5000,text(L"MSD",15,:black,:left))
    annotate!(70,600,text(L"\sim\, t^{3/2}",12,:black,:left))
    pKM

plot(pXY,pKM,size=(800,400))
# savefig("figures\\figures_paper\\MSD_comparaison.svg")

# plot()
#     plot!(rand(10),label=L"0",rib=0)
#     plot!(rand(10),label=L"10^{-2}",rib=0)
#     plot!(rand(10),label=L"10^{-1}",rib=0)
#     plot!(rand(10),label=L"10^{0}",rib=0)
#     plot!(rand(10),label=L"10^{1}",rib=0)
# savefig("figures\\figures_paper\\legend_colors_masses2.svg")
