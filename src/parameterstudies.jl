using LinearAlgebra

  
function comparerates(;T=1, Ns=[100, 1000], gpluss = [0.5, 1,2,4], gminus = 2, alphas = [0.1,1], bs = [0.25,0.75], Nsims= [5_000, 500])
    n_snapshots = size(Ns,1)

    p1 = PDEconstruct()
    p2 = particleconstruct()
    
    function compute(q,gplus,gminus,fplus,fminus,T,p1,p2,dt_PDE,NT,Nsim)
        q2  = merge(q, (;gplus=gplus,gminus=gminus,fplus = fplus,fminus=fminus))
        sol, _ = PDEsolve(T; p = p1, q= q2)
        times = collect(0:dt_PDE:sol.t[end])
        ws_PDE=zeros(size(times,1))
        j=1
        for t in times
            c,w,x = sol2cwx(sol, t)
            ws_PDE[j] = w[1]
            j+=1
        end
        wsaverage_PB = zeros(NT+1)
        for s in 1:Nsim
            ys, ss, xs, ws  = particlesolve(NT, p=p2, q=q2, chosenseed=s)
            wsaverage_PB .+=  1/Nsim*reduce(vcat,transpose.(ws))
        end
        return ws_PDE,wsaverage_PB
    end

    (;dt) = p2
    dt_PB=dt
    (;dt) = p1
    dt_PDE = dt
    NT = Int(round((T/dt_PB)))
    
    plotarray = Any[]
    for j in eachindex(Ns)
        N = Ns[j]
        Nsim = Nsims[j]
        q = parameters(N=N,M=1,initial="init1")
        if j == n_snapshots
            legendspec = :outerright
            ylabel = ""
        elseif j == 1
            legendspec = false
            ylabel="w(t)"
        else
            legendspec = false
            ylabel = ""
        end
        subp=plot(ylim=(0,1.05),xlim=(0,T),xlabel ="t",ylabel=ylabel,legend=legendspec, title = string("n = ",string(N)))
        i=1
        for gplus in gpluss

            function fplus(x) # determines binding rate depending on vesicle occupancy
                if x<1
                    return (1-x)
                else
                    return 0
                end
            end
            fminus(x)=1
            global ws_PDE,wsaverage_PB = compute(q,gplus,gminus,fplus,fminus,T,p1,p2,dt_PDE,NT,Nsim)
            plot!(subp,0:dt_PB:T,wsaverage_PB,label=string(L"\gamma^+ = ",string(gplus)),color=i)
            plot!(subp,0:dt_PDE:T,ws_PDE,label="",color=i,linestyle=:dash)
            i+=1

        end
        push!(plotarray, subp)
    end
    gridp=plot(plotarray..., layout=(1,n_snapshots),size=(n_snapshots*70*5,50*5),link=:all)
 
    for k=2:n_snapshots
        plot!(gridp[k],yformatter=_->"")
    end
    savefig(string("src/img/uncoop.png"))
    savefig(string("src/img/uncoop.pdf"))


    plotarray = Any[]
    for j in eachindex(Ns)
        N = Ns[j]
        Nsim = Nsims[j]
        q = parameters(N=N,M=1,initial="init1")
        if j == n_snapshots
            legendspec = :outerright
            ylabel = ""
        elseif j == 1
            legendspec = false
            ylabel="w(t)"
        else
            legendspec = false
            ylabel = ""
        end
        subp=plot(ylim=(0,1.05),xlim=(0,T),xlabel ="t",ylabel=ylabel,legend=legendspec, title = string("n = ",string(N)))
        i=1
        for gplus in gpluss

            for alpha in alphas
                function fplus(x) # determines binding rate depending on vesicle occupancy
                    if x<1
                        return (1-x)
                    else
                        return 0
                    end
                end
                fminus(x)=1-x+alpha
                ws_PDE,wsaverage_PB = compute(q,gplus,gminus,fplus,fminus,T,p1,p2,dt_PDE,NT,Nsim)
                plot!(subp,0:dt_PB:T,wsaverage_PB,label=string(L"\gamma^+ = ",string(gplus),L", \alpha^-=",string(alpha)),color=i)
                plot!(subp,0:dt_PDE:T,ws_PDE,label="",color=i,linestyle=:dash)
                i+=1
            end

        end
        push!(plotarray, subp)
    end
    gridp=plot(plotarray..., layout=(1,n_snapshots),size=(n_snapshots*70*5,50*5),link=:all)
 
    for k=2:n_snapshots
        plot!(gridp[k],yformatter=_->"")
    end
    savefig(string("src/img/cooplinearunbinding.png"))
    savefig(string("src/img/cooplinearunbinding.pdf"))



    plotarray = Any[]
    for j in eachindex(Ns)
        N = Ns[j]
        Nsim = Nsims[j]
        q = parameters(N=N,M=1,initial="init1")
        if j == n_snapshots
            legendspec = :outerright
            ylabel = ""
        elseif j == 1
            legendspec = false
            ylabel="w(t)"
        else
            legendspec = false
            ylabel = ""
        end
        subp=plot(ylim=(0,1.05),xlim=(0,T),xlabel ="t",ylabel=ylabel,legend=legendspec, title = string("n = ",string(N)))
        i=1
        for gplus in gpluss

            for alpha in alphas
                function fplus(x) # determines binding rate depending on vesicle occupancy
                    if x<1
                        return (1-x)*(x + alpha)
                    else
                        return 0
                    end
                end
                fminus(x)=1
                ws_PDE,wsaverage_PB = compute(q,gplus,gminus,fplus,fminus,T,p1,p2,dt_PDE,NT,Nsim)
                plot!(subp,0:dt_PB:T,wsaverage_PB,label=string(L"\gamma^+ = ",string(gplus),L", \alpha^+=",string(alpha)),color=i)
                plot!(subp,0:dt_PDE:T,ws_PDE,label="",color=i,linestyle=:dash)
                i+=1
            end

        end
        push!(plotarray, subp)
    end
    gridp=plot(plotarray..., layout=(1,n_snapshots),size=(n_snapshots*70*5,50*5),link=:all)
 
    for k=2:n_snapshots
        plot!(gridp[k],yformatter=_->"")
    end
    savefig(string("src/img/cooplinearbinding.png"))
    savefig(string("src/img/cooplinearbinding.pdf"))
    

    plotarray = Any[]
    for j in eachindex(Ns)
        N = Ns[j]
        Nsim = Nsims[j]
        q = parameters(N=N,M=1,initial="init1")
        (;a) = q
        if j == n_snapshots
            legendspec = :outerright
            ylabel = ""
        elseif j == 1
            legendspec = false
            ylabel="w(t)"
        else
            legendspec = false
            ylabel = ""
        end
        subp=plot(ylim=(0,1.05),xlim=(0,T),xlabel ="t",ylabel=ylabel,legend=legendspec, title = string("n = ",string(N)))
        
        i=1
        for gplus in gpluss

            for b in bs
                function fplus(x) # determines binding rate depending on vesicle occupancy
                    if x<1
                        return (1-x) 
                    else
                        return 0
                    end
                end
                fminus(x)=b^(a*N*x-1)
                ws_PDE,wsaverage_PB = compute(q,gplus,gminus,fplus,fminus,T,p1,p2,dt_PDE,NT,Nsim)
                plot!(subp,0:dt_PB:T,wsaverage_PB,label=string(L"\gamma^+ = ",string(gplus),L", \beta=",string(b)),color=i)
                plot!(subp,0:dt_PDE:T,ws_PDE,label="",color=i,linestyle=:dash)
                i+=1
            end

        end
        push!(plotarray, subp)
    end
    gridp=plot(plotarray..., layout=(1,n_snapshots),size=(n_snapshots*70*5,50*5),link=:all)
 
    for k=2:n_snapshots
        plot!(gridp[k],yformatter=_->"")
    end
    savefig(string("src/img/coopexpunbinding.png"))
    savefig(string("src/img/coopexpunbinding.pdf")) 
end

