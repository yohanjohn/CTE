using CairoMakie, Random, LinearAlgebra

Ncortices = 5

Ctx = zeros(3,Ncortices)

# Presence of tau in first cortical area 
Ctx[1,1] = 1.0
Ctx[3,1] = 1.0

# Range of cortical types
Ctype = LinRange(0,1,Ncortices)

pset = []


infra = 1
supra = 3

# Connection strength function
function W(receiver::Float64, sender::Float64)
    M = zeros(3,3)

    amp = exp(-(receiver - sender)^2)

    M[infra,supra] = amp * (1 - max(receiver - sender,0))
    M[supra,infra] = amp * (1 - max(sender - receiver,0))

    return M
end

ylab = ["infra", "gran","supra"]

fig = Figure(size = (1200, 1000), fonts = (; regular = "Arial"))


toshow = [1 5 10 30]
t = 0

for i = 1:30
    global t

    if i in toshow
        t = t + 1
        ax = Axis(fig[1, t], yticks = (1:3, ["L5/6", "L4", "L2/3"]))

        heatmap!(ax,  Ctx', colormap = :heat)

        if t==1
            ax2 = Axis(fig[2, t], ylabel = "SGI")
        else
            ax2 = Axis(fig[2, t])
        end

        trend = Ctx[3,:] ./ (Ctx[1,:] .+ Ctx[3,:] .+ 0.001)
        lines!(ax2, LinRange(0.5,5.5,Ncortices),0.5 .* ones(Ncortices),linestyle = :dash, color = :red)
        scatterlines!(ax2,trend)
        ylims!(ax2,-0.1,0.9)
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
        ax3.title = "C-" * string(t)
    end

    # Spread of tau per time step
    for r = 1:Ncortices
        for s = 1:Ncortices
            if s!=r
                Ctx[:,r] = Ctx[:,r] .+ 0.1.*(1.0 .- Ctx[:,r]) .* (W(Ctype[r],Ctype[s])*Ctx[:,s])
            end
        end
    end
 

end

cb = Colorbar(fig[1, 5]; limits=(0,1), label = "Proportion of tau", colormap = :heat)
fig

save("CTE_SimResult.eps", fig; pt_per_unit=1)