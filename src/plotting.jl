using Plots
pyplot()
cmap =reverse(cgrad(:gist_earth)) 
markercolors_vesicle = :binary

function PDEplot(c,w,x,(p,q),t; clim=(0,3), title = string("t=", string(round(t, digits=2))))
    (; domain, dx) = p

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dx:domain[2,2]

    subp = heatmap(x_arr,y_arr, c, title = title, c=cmap, clim=clim)
    
    # plot vesicle
    scatter!(subp, x[:,1],x[:,2],markerstrokecolor=:white, markersize=10,c=:black,legend=false)  
    return subp
end

function PDEgif(sol, (p,q), dt=0.01; save=true, name = "")
    T = 0
    anim = Animation()
    for t in 0:dt:sol.t[end]
        c,w,x = sol2cwx(sol, t)
        plt = PDEplot(c,w,x,(p,q),t+T)
        frame(anim, plt)
    end

    if save==true
        Plots.gif(anim, string("src/img/pde",name,".gif"), fps = 10)
    end
end

function particleplot(y,s, x, w,(p,q),t; binnumber = 10, clim=(0,3),title = string("t=", string(round(t, digits=2))))
    (; domain) = p
    (;N ) = q

    hist,xrange,yrange = particlehistogram(y,s,domain,binnumber)
    dV= (xrange[2]-xrange[1])* (yrange[2]-yrange[1])
    hist=1/(N*dV)*hist
    subp = heatmap(xrange, yrange,hist', title = title, c=cmap, clim=clim)
    
    # plot vesicle
    scatter!(subp, x[:,1],x[:,2],markerstrokecolor=:white, markersize=10,c=:black,legend=false)  
    return subp
end

function particlegif(ys, xs, ss, ws, (p,q); dN=10, save=true, name = "")
    (; dt) = p
    anim = Animation()
    for n in 1:dN:size(ws,1)
            y = ys[n]
            x = xs[n]
            w = ws[n]
            s = ss[n]

            plt = particleplot(y, s, x, w,(p,q),n*dt)
            frame(anim, plt)
    end

    if save==true
        Plots.gif(anim, string("src/img/particle",name,".gif"), fps = 10)
    end
end

function particleoccupancy(ws,(p,q); save=true, name = "")
    (; dt) = p
    (; M) = q
    ws = reduce(vcat,transpose.(ws))
    subp = plot()
    for m in 1:M
        plot!(subp,range(0,size(ws,1)-1)*dt,ws[:,m],ylim=(0,1.1),label="relative occupancy",xlabel ="t",ylabel="w(t)")
    end
    
    if save==true
        savefig(string("src/img/particlemodel_occupancy_",name,".png"))
    end
end

function PDEoccupancy(sol,(p,q); dt=0.02, save=true, name = "")
    (;M) = q
    times = collect(0:dt:sol.t[end])
    ws=zeros(size(times,1),M)
    i=1
    for t in times
        c,w,x = sol2cwx(sol, t)
        ws[i,:] = w
 
        i+=1
    end

    subp = plot()
    for m in 1:M
        plot!(subp,times,ws[:,m],ylim=(0,1.1),label="relative occupancy",xlabel ="t",ylabel="w(t)")
    end

    if save==true
        savefig(string("src/img/PDE_occupancy_",name,".png"))
    end
    
end 