using Plots
pyplot()
cmap =reverse(cgrad(:gist_earth)) 
markercolors_vesicle = [:seagreen :firebrick ]

function PDEplotsingle(c,v,x,(p,q),t; clim=(0,15), title = string("t=", string(round(t, digits=2))))
    (; domain, dx) = p


    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dx:domain[2,2]


    subp = heatmap(x_arr,y_arr, c, title = title, c=cmap, clim=clim, legend=false)

    # plot vesicle
    bound = v[1]
    scatter!(subp, [x[1]],[x[2]],markercolor = markercolors_vesicle[round(Int,bound+1)],markersize=10)  
    return subp
end

PDEgifsingle(sol, (p,q), args...; kwargs...) = PDEgifsingle([sol], [(p,q)], args...; kwargs...)
function PDEgifsingle(sols::Vector, Ps::Vector, dt=0.1; save=true, name = "")
    T = 0
    anim = Animation()
    for (sol, P) in zip(sols, Ps)
        for t in 0:dt:sol.t[end]
            c,v,x = sol2cvx(sol, t)
            plt = PDEplotsingle(c,v,x,P,t+T)
            frame(anim, plt)
        end
        T += sol.t[end]
    end
    if save==true
        Plots.gif(anim, string("src/img/pde",name,".gif"), fps = 10)
    end
end