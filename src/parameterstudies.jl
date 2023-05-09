using LinearAlgebra

  
function comparerates(T=1, N=100; gpluss = [0.5, 1,2,4], gminus = 2, alphas = [0.1,1], bs = [0.25,0.75], Nsim= 1000)
    q = parameters(N=N,M=1,initial="init1")
    p1 = PDEconstruct()
    p2 = particleconstruct()
    legendloc = :outerright
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
    
    pl=plot(ylim=(0,1.05),xlim=(0,T),xlabel ="t",ylabel="w(t)",legend=legendloc)
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
        plot!(pl,0:dt_PB:T,wsaverage_PB,label=string(L"\gamma^+ = ",string(gplus)),color=i)
        plot!(pl,0:dt_PDE:T,ws_PDE,label="",color=i,linestyle=:dash)
        i+=1

    end
    savefig(string("src/img/uncoop",string(N),".png"))
    savefig(string("src/img/uncoop",string(N),".pdf"))

    pl=plot(ylim=(0,1.05),xlim=(0,T),xlabel ="t",ylabel="w(t)",legend=legendloc)
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
            plot!(pl,0:dt_PB:T,wsaverage_PB,label=string(L"\gamma^+ = ",string(gplus),L", \alpha^-=",string(alpha)),color=i)
            plot!(pl,0:dt_PDE:T,ws_PDE,label="",color=i,linestyle=:dash)
            i+=1
        end

    end
    savefig(string("src/img/cooplinearunbinding",string(N),".png"))
    savefig(string("src/img/cooplinearunbinding",string(N),".pdf"))

    pl=plot(ylim=(0,1.05),xlim=(0,T),xlabel ="t",ylabel="w(t)",legend=legendloc)
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
            plot!(pl,0:dt_PB:T,wsaverage_PB,label=string(L"\gamma^+ = ",string(gplus),L", \alpha^+=",string(alpha)),color=i)
            plot!(pl,0:dt_PDE:T,ws_PDE,label="",color=i,linestyle=:dash)
            i+=1
        end

    end
    savefig(string("src/img/cooplinearbinding",string(N),".png"))
    savefig(string("src/img/cooplinearbinding",string(N),".pdf"))
    
    (;a) = q
    pl=plot(ylim=(0,1.05),xlim=(0,T),xlabel ="t",ylabel="w(t)",legend=legendloc)
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
            plot!(pl,0:dt_PB:T,wsaverage_PB,label=string(L"\gamma^+ = ",string(gplus),L", \beta=",string(b)),color=i)
            plot!(pl,0:dt_PDE:T,ws_PDE,label="",color=i,linestyle=:dash)
            i+=1
        end

    end
    savefig(string("src/img/coopexpunbinding",string(N),".png"))
    savefig(string("src/img/coopexpunbinding",string(N),".pdf"))
end

