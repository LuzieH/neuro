using Plots
pyplot()
cmap =  :YlGnBu #:roma #  cmap =:PuBuGn   cmap = :GnBu
markercolors_vesicle = :binary
clim = (0,2)

function PDEplot(c,w,x,(p,q),t; clim=clim, title = string("t=", string(round(t, digits=2))),ylabel="",legend=false)
    (; domain, dx) = p

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dx:domain[2,2]

    subp = heatmap(x_arr,y_arr, c, title = title, c=cmap, clim=clim)
    
    # plot vesicle
    scatter!(subp, x[:,1],x[:,2],markerstrokecolor=:white, markersize=10,c=:black,legend=legend,ylabel=ylabel,label="Vesicles",xlim = (domain[1,1],domain[1,2]),ylim = (domain[2,1],domain[2,2]))  
    return subp
end

function PDEgif(sol, (p,q); dt=0.01, save=true, name = "")
    T = 0
    anim = Animation()
    for t in 0:dt:sol.t[end]
        c,w,x = sol2cwx(sol, t)
        plt = PDEplot(c,w,x,(p,q),t+T)
        frame(anim, plt)
    end

    if save==true
        Plots.gif(anim, string("img/pde",name,".gif"), fps = 10)
    end
end


function PDEoccupancy(sol,(p,q); dt=0.02, save=true, name = "")
    (;M) = q
    @assert M>=2
    labels = [L"W_1", L"W_2"]
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
        plot!(subp,times,ws[:,m],ylim=(0,1.1),label=labels[m],xlabel ="t")
    end

    if save==true
        savefig(string("img/PDEoccupancy",name,".png"))
    end
    
end 


function PDEsnapshots(sol, (p,q), ts; save = true, name="",clim = clim)
    nsnapshots = length(ts)
    plotarray = Any[]  
    legend = true
    for i in 1:nsnapshots
        c,w,x = sol2cwx(sol, ts[i])

        ylabel =string("t = ", string(round(ts[i], digits=2)))
        if i>1
            legend=false
        end
        subp=  PDEplot(c,w,x,(p,q),ts[i],legend=legend,ylabel=ylabel,title="",clim=clim)
        
        push!(plotarray, subp)
    end    
    gridp=plot(plotarray..., layout=(nsnapshots,1),size=(95*5,nsnapshots*50*5),link=:all)
 

    for k=1:nsnapshots-1
        plot!(gridp[k],xformatter=_->"")
    end

    if save==true
        savefig(string("img/PDEsnapshots",name,".png"))
        savefig(string("img/PDEsnapshots",name,".pdf"))
    end
    return gridp
end


###

function particleplot(y,s, x, w,(p,q),t; binnumber = 10, clim=clim,title = string("t=", string(round(t, digits=2))))
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


function particleplotscatter(y,s, x, w,(p,q),t; binnumber = 10, clim=clim,title = string("t=", string(round(t, digits=2))),legend=false,ylabel="")
    (; domain) = p
    (;N ) = q

    indices = findall(s.==0)
    subp = scatter(y[indices,1],y[indices,2],markerstrokecolor=:white, markersize=5,title=title,label="Calcium ions",xlim = (domain[1,1],domain[1,2]),ylim = (domain[2,1],domain[2,2]))
    
    # plot vesicle
    scatter!(subp, x[:,1],x[:,2],markerstrokecolor=:white, markersize=10,c=:black,label="Vesicles",legend=legend,ylabel=ylabel)  
    return subp
end

function particlegif(ys, xs, ss, ws, (p,q); dN=10, save=true, name = "",plotfunction = particleplotscatter)
    (; dt) = p
    anim = Animation()
    for n in 1:dN:size(ws,1)
            y = ys[n]
            x = xs[n]
            w = ws[n]
            s = ss[n]

            plt = plotfunction(y, s, x, w,(p,q),n*dt)
            frame(anim, plt)
    end

    if save==true
        Plots.gif(anim, string("img/particle",name,".gif"), fps = 10)
    end
end

function particleoccupancy(ws,(p,q); save=true, name = "")
    (; dt) = p
    (; M) = q
    @assert M>=2
    labels = [L"W_1", L"W_2"]
    ws = reduce(vcat,transpose.(ws))
    subp = plot()
    for m in 1:M
        plot!(subp,range(0,size(ws,1)-1)*dt,ws[:,m],ylim=(0,1.1),label=labels[m],xlabel ="t")
    end
    
    if save==true
        savefig(string("img/particleoccupancy",name,".png"))
    end
end


function particlesnapshots(ys, xs, ss, ws, (p,q), ts; save = true, name="",plotfunction = particleplotscatter)
    (;dt) = p
    nsnapshots = length(ts)
    plotarray = Any[]  
    legend = true
    for i in 1:nsnapshots
        t=ts[i]
        y=ys[t]
        x=xs[t]
        s=ss[t]
        w=ws[t]
        ylabel =string("t = ", string(round((t-1)*dt, digits=2)))
        if i>1
            legend=false
        end
        subp= plotfunction(y, s, x, w,(p,q),t; binnumber = 20, clim=clim,title = "",legend=legend,ylabel=ylabel)
        push!(plotarray, subp)
    end    
    gridp=plot(plotarray..., layout=(nsnapshots,1),size=(95*5,nsnapshots*50*5),link=:all)
 

    for k=1:nsnapshots-1
        plot!(gridp[k],xformatter=_->"")
    end

    if save==true
        savefig(string("img/particlesnapshots",name,".png"))
        savefig(string("img/particlesnapshots",name,".pdf"))
    end
    return gridp
end
