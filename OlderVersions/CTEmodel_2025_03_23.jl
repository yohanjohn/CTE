using CairoMakie, Random, LinearAlgebra

Ncortices = 5

Ctx = zeros(3,Ncortices)

infra = 1
supra = 3

# Presence of tau in first cortical area 
Ctx[infra,1] = 0.6
Ctx[supra,1] = 0.0

# Range of cortical types
Ctype = LinRange(0,1,Ncortices)

pset = []


# Connection strength function
function W(receiver::Float64, sender::Float64)
    M = zeros(3,3)

    #amp = exp(-(receiver - sender)^2)
    amp = 1.0
   # M[infra,supra] = 0.02*amp * 1.0/(1.0 + 2.5*max(receiver - sender,0))
   # M[supra,infra] = amp * 1.0/(1.0 + 2.5*max(sender - receiver,0))

    M[infra,supra] = 0.1 * amp * (1 - max(receiver - sender,0))
    M[supra,infra] = 0.1 * amp * (1 - max(sender - receiver,0))
    #M[infra,infra] = amp * (max(receiver - sender,0))
    #M[supra,supra] = amp * (max(sender - receiver,0))
    #M[infra,infra] = 0.1*amp  
    #M[supra,supra] = 0.1*amp 

    return M
end

ylab = ["infra", "gran","supra"]

fig = Figure(size = (1200, 1000), fonts = (; regular = "Arial"))


toshow = [1 5 10 30]
t = 0
hmaps = []
trends = []
densities = []
cmax = 0.0

for i = 1:toshow[end]
    global t, cmax, hmaps, trends, densities

    if i in toshow
        t = t + 1
        ax = Axis(fig[1, t], yticks = (1:3, ["L5/6", "L4", "L2/3"]))
        heatmap!(ax,  Ctx', colormap = :heat)
        #heatmap!(ax,  Ctx', colorrange=(0.0,100.0), colormap = :heat)
        cmax = max(maximum(Ctx),cmax)

        push!(hmaps, Ctx')

        if t==1
            ax2 = Axis(fig[2, t], ylabel = "SGI")
        else
            ax2 = Axis(fig[2, t])
        end

        trend = Ctx[3,:] ./ (Ctx[1,:] .+ Ctx[3,:] .+ 0.001)
        push!(trends,trend)

        lines!(ax2, LinRange(0.5,5.5,Ncortices),0.5 .* ones(Ncortices),linestyle = :dash, color = :red)
        scatterlines!(ax2,trend)
        ylims!(ax2,-0.1,1.1)
        xlims!(ax2,0.5,5.5)
        if t==1
            ax3 = Axis(fig[3, t], xlabel = "Cortical type", ylabel = "Density (a.u.)")
        else
            ax3 = Axis(fig[3, t], xlabel = "Cortical type")
        end
        push!(densities,Ctx[3,:] .+ Ctx[1,:])
        scatterlines!(ax3,Ctx[3,:] .+ Ctx[1,:])
        ylims!(ax3,-0.1,2.2)

        ax.title = "A-" * string(t)
        ax2.title = "B-" * string(t)
        ax3.title = "C-" * string(t)
    end

    # Spread of tau per time step
    for r = 1:Ncortices
        for s = 1:Ncortices
            if s!=r
                #Ctx[:,r] = Ctx[:,r] .+ 0.2.*(1.0 .- Ctx[:,r]) .* (W(Ctype[r],Ctype[s])*max.(Ctx[:,s] .- 0.2,0))
                #Ctx[:,r] = Ctx[:,r] .+ 0.005.*(100.0 .- Ctx[:,r]) .* (W(Ctype[r],Ctype[s]) * Ctx[:,s] )
                Ctx[:,r] = Ctx[:,r] .+ 0.5.* (W(Ctype[r],Ctype[s]) * Ctx[:,s] )
                #Ctx[:,r] = Ctx[:,r] .+ 0.015.*(100.0 .- Ctx[:,r]) .* W(Ctype[r],Ctype[s]) * max.(Ctx[:,s] .- 0.5,0)     

                #Ctx[:,r] = Ctx[:,r] .+ 0.10.*((1 .- Ctx[:,r]) .* W(Ctype[r],Ctype[s]) * Ctx[:,s] - 0.3*Ctx[:,r])
                #Ctx[:,r] = Ctx[:,r] .+ 0.2.* (W(Ctype[r],Ctype[s]) * Ctx[:,s] )
            end
        end
    end
 

end

i = 0
for t in toshow
    global i
    i = i + 1
    ax = Axis(fig[1, t], yticks = (1:3, ["L5/6", "L4", "L2/3"]))
    #heatmap!(ax,  Ctx', colormap = :heat)
    #heatmap!(ax,  Ctx', colorrange=(0.0,100.0), colormap = :heat)
  

    if t==1
        ax2 = Axis(fig[2, t], ylabel = "SGI")
    else
        ax2 = Axis(fig[2, t])
    end

    trend = Ctx[3,:] ./ (Ctx[1,:] .+ Ctx[3,:] .+ 0.001)
    lines!(ax2, LinRange(0.5,5.5,Ncortices),0.5 .* ones(Ncortices),linestyle = :dash, color = :red)
    scatterlines!(ax2,trend)
    ylims!(ax2,-0.1,1.1)
    xlims!(ax2,0.5,5.5)
    if t==1
        ax3 = Axis(fig[3, t], xlabel = "Cortical type", ylabel = "Density (a.u.)")
    else
        ax3 = Axis(fig[3, t], xlabel = "Cortical type")
    end
    scatterlines!(ax3,Ctx[3,:] .+ Ctx[1,:])
    ylims!(ax3,-0.1,2.2)

    ax.title = "A-" * string(t)
    ax2.title = "B-" * string(t)
    ax3.t
end
mx = maximum(Ctx)
cb = Colorbar(fig[1, 5]; colorrange=(0.0,mx), label = "Proportion of tau", colormap = :heat)
fig

#save("CTE_SimResult_Original2.png", fig; pt_per_unit=1)
#save("CTE_SimResult.tif", fig; pt_per_unit=1)