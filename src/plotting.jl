using Plots
pyplot()
grid=false
cmap =  :YlGnBu #:roma #  cmap =:PuBuGn   cmap = :GnBu
markercolors_vesicle = :binary
clim = (0,2)
bs =0.05
radius = 0.1
arcshape(θ1, θ2) = Shape(vcat(Plots.partialcircle(θ1, θ2, 50, radius),
                      reverse(Plots.partialcircle(θ1, θ2, 50, 0.01*radius))))
fullcircle = Shape(Plots.partialcircle(0, 2*π, 50, radius))
plotxsize = 390
plotxsizeparticles = 330
plotysize = 200
plotxsizeocc = 330
plotysizeocc = 250
plotxsizeoccrates = plotxsizeocc*1.1
plotysizeoccrates = plotysizeocc*1.1
dpi=300

function PDEplot(c,w,x,(p,q),t; clim=clim, title = string("t=", string(round(t, digits=2))),ylabel="",legend=false)
    (; domain, dx) = p
    (;M ) = q

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dx:domain[2,2]

    subp = heatmap(x_arr,y_arr, c, title = title, c=cmap, clim=clim,grid=grid,dpi=dpi)
    
    # plot vesicle
    #scatter!(subp, x[:,1],x[:,2],markerstrokecolor=:white, markersize=10,c=:black,legend=legend,ylabel=ylabel,label="Vesicles",xlim = (domain[1,1]-bs,domain[1,2]+bs),ylim = (domain[2,1]-bs,domain[2,2]+bs))  
    for m in 1:M
        if m==1
            plot!(subp,[x[m,1]],[x[m,2]],marker = (8, fullcircle),linewidth = 1, c=:white,ylabel=ylabel,label="Vesicles",legend=legend)
        else
            plot!(subp,[x[m,1]],[x[m,2]],marker = (8, fullcircle),linewidth = 1, c=:white,ylabel=ylabel,label=false,legend=legend)
        end
        if w[m]>0
            plot!(subp,[x[m,1]],[x[m,2]],marker = (7, arcshape(0,w[m]*2*π)),ylabel=ylabel,c=:black,fillcolor = :black)
        end
    end

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
    labels = [L"w_1", L"w_2"]
    times = collect(0:dt:sol.t[end])
    ws=zeros(size(times,1),M)
    i=1
    for t in times
        c,w,x = sol2cwx(sol, t)
        ws[i,:] = w
        i+=1
    end

    subp = plot(grid=grid,size=(plotxsizeocc,plotysizeocc),dpi=dpi)
    for m in 1:M
        plot!(subp,times,ws[:,m],ylim=(0,1.1),label=labels[m],xlabel ="t")
    end

    if save==true
        savefig(string("img/PDEoccupancy",name,".png"))
        savefig(string("img/PDEoccupancy",name,".pdf"))
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
    gridp=plot(plotarray..., layout=(nsnapshots,1),size=(plotxsize,nsnapshots*plotysize),link=:all,dpi=dpi)
 

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
    subp = heatmap(xrange, yrange,hist', title = title, c=cmap, clim=clim,grid=grid,dpi=dpi)
    
    # plot vesicle
    scatter!(subp, x[:,1],x[:,2],markerstrokecolor=:white, markersize=10,c=:black,legend=false)  
    return subp
end


function particleplotscatter(y,s, x, w,(p,q),t; binnumber = 10, clim=clim,title = string("t=", string(round(t, digits=2))),legend=false,ylabel="")
    (; domain) = p
    (;M ) = q

    indices = findall(s.==0)
    subp = scatter(y[indices,1],y[indices,2],markerstrokewidth=0, markersize=5,title=title,label="Calcium ions",xlim = (domain[1,1]-bs,domain[1,2]+bs),ylim = (domain[2,1]-bs,domain[2,2]+bs),grid=grid)
    
    # plot vesicle
    for m in 1:M
        if m==1
            plot!(subp,[x[m,1]],[x[m,2]],marker = (8, fullcircle),linewidth = 1, c=:white,ylabel=ylabel,label="Vesicles",legend=legend)
        else
            plot!(subp,[x[m,1]],[x[m,2]],marker = (8, fullcircle),linewidth = 1, c=:white,ylabel=ylabel,label=false,legend=legend)
        end
        if w[m]>0
            plot!(subp,[x[m,1]],[x[m,2]],marker = (7, arcshape(0,w[m]*2*π)),ylabel=ylabel,c=:black,fillcolor = :black)
        end
    end

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
    labels = [L"w_1", L"w_2"]
    ws = reduce(vcat,transpose.(ws))
    subp = plot(grid=grid,size=(plotxsizeocc,plotysizeocc),dpi=dpi)
    times = range(0,size(ws,1)-1)*dt
    for m in 1:M
        plot!(subp,times,ws[:,m],ylim=(0,1.1),label=labels[m],xlabel ="t")
    end
    
    if save==true
        savefig(string("img/particleoccupancy",name,".png"))
        savefig(string("img/particleoccupancy",name,".pdf"))
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
    gridp=plot(plotarray..., layout=(nsnapshots,1),size=(plotxsizeparticles,nsnapshots*plotysize),link=:all,dpi=dpi)
 

    for k=1:nsnapshots-1
        plot!(gridp[k],xformatter=_->"")
    end

    if save==true
        savefig(string("img/particlesnapshots",name,".png"))
        savefig(string("img/particlesnapshots",name,".pdf"))
    end
    return gridp
end



"""plot ensemble of particle-dynamics"""
function ensembleplot(meanhist, wsaverage, xsaverage, xrange, yrange, ts, (p,q);clim=(minimum(meanhist)-0.1,maximum(meanhist)+0.1), name="", save=true)
    (;dt, domain) = p 
    (;M) = q
    nsnapshots = size(ts,1)
    particleoccupancy(wsaverage,(p,q),name="ensemble")
    

    # plot snapshots
    plotarray = Any[]  
    
    for i in eachindex(ts)
        t=ts[i]
        x=xsaverage[t]
        w=wsaverage[t]
        ylabel =string("t = ", string(round((t-1)*dt, digits=2)))
        subp = heatmap(xrange, yrange,meanhist[i,:,:]', c=cmap, clim=clim,xlim = (domain[1,1]-bs,domain[1,2]+bs),ylim = (domain[2,1]-bs,domain[2,2]+bs),grid=grid,dpi=dpi)

        if i>1
            labelscatter = ""
            legendscatter = false
        else
            legendscatter = true
            labelscatter = "Vesicles"
        end
        # plot vesicle
        for m in 1:M
            if m==1
                plot!(subp,[x[m,1]],[x[m,2]],marker = (8, fullcircle),linewidth = 1, c=:white,ylabel=ylabel,label=labelscatter,legend=legendscatter)
            else
                plot!(subp,[x[m,1]],[x[m,2]],marker = (8, fullcircle),linewidth = 1, c=:white,ylabel=ylabel,label=false,legend=legendscatter)
            end
            if w[m]>0
                plot!(subp,[x[m,1]],[x[m,2]],marker = (7, arcshape(0,w[m]*2*π)),ylabel=ylabel,c=:black,fillcolor = :black)
            end
        end
        push!(plotarray, subp)
    end    
    gridp=plot(plotarray..., layout=(nsnapshots,1),size=(plotxsize,nsnapshots*plotysize),link=:all,dpi=dpi)
 

    for k in 1:nsnapshots-1
        plot!(gridp[k],xformatter=_->"")
    end

    if save==true
        savefig(string("img/ensemblesnapshots",name,".png"))
        savefig(string("img/ensemblesnapshots",name,".pdf"))
    end
end