using LinearAlgebra



"""compare PDE simulations vs. average (Nsim Monte-Carlo simulations) particle-based simulations with different numbers of ions (given by Ns)"""
function comparison(NT=1000,Nsim=100,Ns=[100,1000,10000];  q = parameters(), p1 = PDEconstruct(), p2 = particleconstruct(),binnumber = 20,histogram=false)
    (;domain) = p1
    (;M ) = q
    ws_averages=Any[]
    hist_averages=Any[]

    for k in eachindex(Ns)
        N=Ns[k]
        q = merge(q, (;N=N))
        ws_average = zeros(NT+1,M)
        if histogram==true
            hist_average = Any[]
            hist_average = zeros(binnumber,binnumber)
        end
        for n in 1:Nsim
            ys, ss, xs, ws = particlesolve(NT, p=p2, q=q,chosenseed=n)
            ws = reduce(vcat,transpose.(ws))
            ws_average .+=1/Nsim* ws
            if histogram==true
                hist,xrange,yrange = particlehistogram(ys[end],ss[end],domain, binnumber) 
                hist_average += 1/Nsim*hist
            end
        end
         
        ws_averages=push!(ws_averages,ws_average)
        if histogram==true
            hist_averages=push!(hist_averages,hist_average)
            # plot average final concentration of particle-dynamics
            heatmap(hist_average',title=string("N = ",string(round(N))))
            savefig(string("src/img/densabm_",string(round(N)),".png"))
        end
        

    end
 
    (; dt) = p2
    sol, (p1,q) = PDEsolve(NT*dt; p = p1, q= q)
    times = collect(0:dt:sol.t[end])
    ws_pde=zeros(size(times,1),M)
    i=1
    for t in times
        c,w,x = sol2cwx(sol, t)
        ws_pde[i,:] = w
        i+=1
    end
    if histogram==true
        # plot final concentration of ions in PDE
        c,w,x = sol2cwx(sol, NT*dt)
        heatmap(c',title="PDE")
        savefig(string("src/img/denspde.png"))
    end

    # plot relative occupancy in time
    subp = plot()
    for m in 1:M
        plot!(subp,times,ws_pde[:,m],ylim=(0,1.1),label=string("PDE, m =",string(m)),linestyle=:dash,legend=:bottomright,xlabel ="t",ylabel="average w(t)")
        for k in eachindex(Ns)
            N=Ns[k]
            plot!(subp,range(0,size(ws_averages[k][:,m],1)-1)*dt,ws_averages[k][:,m],ylim=(0,1.1),label=string("n = ",string(round(N)),", m = ",string(m))) #range(0,size(ws_averages[k][:,m],1)-1)*dt,
        end
    end
    savefig(string("src/img/comparison.png"))  
    savefig(string("src/img/comparison.pdf")) 
end


function studydiscretization(T=0.5,Nsim=100,N=10000,dts = [0.001, 0.0005, 0.00025], dxs = [0.01, 0.025]; q = parameters(N=N), p1 = PDEconstruct(), p2 = particleconstruct())
    (; dt) = p1
    (;M ) = q
    dt_PDE = dt
    NT_PDE = Int(T/dt_PDE)
    lw = 0.5
    subp = plot() 
    # simulate particle based for different dt
    ws_averages=Any[]
    for dt in dts
        p2 = merge(p2, (;dt=dt))
        NT = Int(round(T/dt))
        ws_average = zeros(NT_PDE+1,M)  
        for n in 1:Nsim        
            ys, ss, xs, ws = particlesolve(NT, p=p2, q=q, chosenseed=n)
            ws = reduce(vcat,transpose.(ws))
            ws_average .+=1/Nsim* ws[1:Int(dt_PDE/dt):end,:]
        end
        ws_averages=push!(ws_averages,ws_average)
        for m in 1:M
            plot!(0:dt_PDE:T,ws_average[:,m],ylim=(0,1),label=string("PB, dt=",string(dt), ", m = ",string(m)),legend=:bottomright,xlabel ="t",ylabel="w(t)",linewidth=lw)
        end
    end

    PDE_wss =Any[]
    for dx in dxs
        p1 = merge(p1,(;dx=dx))
        sol, (p1,q) = PDEsolve(T;  p = p1, q= q)
        times = collect(0:dt_PDE:sol.t[end])
        PDE_ws=zeros(size(times,1),M)
        i=1
        for t in times
            c,w,x = sol2cwx(sol, t)
            PDE_ws[i,:] = w
            i+=1
        end
        PDE_wss=push!(PDE_wss,PDE_ws)
        for m in 1:M
            plot!(times,PDE_ws[:,m],ylim=(0,1),label=string("PDE, dx=",string(dx), ", m = ", string(m)),linewidth=lw)
        end
    end
    savefig(string("src/img/studydiscretization.png"))  
    savefig(string("src/img/studydiscretization.pdf"))  
    return ws_averages, PDE_wss
end

#studyparameters(2,200,100,[0.05, 0.075, 0.1, 0.1785, 0.25], [1, 2.5, 5, 10, 20], [10, 20, 30, 40,50], [1/50, 1/20, 1/10, 1/5]; q=parameters(N=100))
#studyparameters(1,200,100,[0.05,  0.1, 0.1785, 0.25], [ 2.5, 5, 10, 20], [10, 25, 35, 50], [1/50, 1/20, 1/10, 1/5]; q=parameters(N=100))

function studyparameters(T=1,Nsim=1000,N=100,eps=[0.05, 0.1, 0.1785, 0.25],gplus=[2.5, 5, 10, 20],gminus=[5, 10, 15, 20], as = [1/50, 1/20, 1/10]; q = parameters(N=N), p1 = PDEconstruct(), p2 = particleconstruct())
    (; dt) = p1
    dt_PDE = dt
    (; dt) = p2
    dt_particle = dt
    NT_PDE = Int(T/dt_PDE)
    NT = Int(round(T/dt_particle))
    (;M) = q
    
    function ws_compare(q,p1,p2,NT,NT_PDE,Nsim,T,dt_PDE,dt_particle,labelparticle=string("particlebased, eps=",string(e)),labelpde=string("PDE, eps=",string(e)),colorindex=1)
        particles_ws = zeros(NT_PDE+1,M)  
        #particle based
        for n in 1:Nsim        
            ys, ss, xs, ws = particlesolve(NT, p=p2, q=q, chosenseed=n)
            ws = reduce(vcat,transpose.(ws))
            particles_ws .+=1/Nsim*ws[1:Int(dt_PDE/dt_particle):end,:]
        end
        #pde based
        sol, (p1,q) = PDEsolve(T;  p = p1, q= q)
        times = collect(0:dt_PDE:sol.t[end])
        PDE_ws=zeros(size(times,1),M)
        j=1
        for t in times
            c,w,x = sol2cwx(sol, t)
            PDE_ws[j,:] = w
            j+=1
        end
        for m in 1:M
            plot!(0:dt_PDE:T,particles_ws[:,m],ylim=(0,1.05),label=labelparticle,legend=:bottomright,xlabel ="t",ylabel="w(t)",color=colorindex)
            plot!(0:dt_PDE:T,PDE_ws[:,m],ylim=(0,1.05),label=labelpde,color=colorindex,linestyle=:dash)
        end
    end


    subp = plot() 
    i=1
    for e in eps
        q1 = merge(q, (;eps=e))
        ws_compare(q1,p1,p2,NT,NT_PDE,Nsim,T,dt_PDE,dt_particle,string("PB, eps=",string(e)),string("PDE, eps=",string(e)),i)
        i+=1
    end
    savefig(string("src/img/studyeps.png")) 
    savefig(string("src/img/studyeps.pdf")) 

    subp = plot() 
    i=1
    for gp in gplus
        q2 = merge(q, (;gplus=gp))
        ws_compare(q2,p1,p2,NT,NT_PDE,Nsim,T,dt_PDE,dt_particle,string("PB, g+=",string(gp)),string("PDE, g+=",string(gp)),i)
        i+=1
    end
    savefig(string("src/img/studygammaplus.png")) 
    savefig(string("src/img/studygammaplus.pdf")) 

    subp = plot() 
    i=1
    for gm in gminus
        q3 = merge(q, (;gminus=gm))
        ws_compare(q3,p1,p2,NT,NT_PDE,Nsim,T,dt_PDE,dt_particle,string("PB, g-=",string(gm)),string("PDE, g-=",string(gm)),i)
        i+=1
    end
    savefig(string("src/img/studygammaminus.png")) 
    savefig(string("src/img/studygammaminus.pdf")) 

    subp = plot() 
    i=1
    for a in as
        q4 = merge(q, (;a=a))
        ws_compare(q4,p1,p2,NT,NT_PDE,Nsim,T,dt_PDE,dt_particle,string("PB, a=",string(a)),string("PDE, a=",string(a)),i)
        i+=1
    end
    savefig(string("src/img/studya.png")) 
    savefig(string("src/img/studya.pdf")) 
end


function compare1box(T,Nions,Nvesicle,Nsim=100;q=parameters(),dt_save=0.01,dt_solve=0.001,save=false)
    (;gplus,gminus,fplus,fminus)=q
    NT_save = Int(T/dt_save)
    NT_solve = Int(T/dt_solve)
    
    #particlebased
    Nbound_average =zeros(NT_save+1)
    for sim in 1:Nsim
        s=zeros(Nions)
        c=2
        for k in 2:NT_solve+1
            for i in shuffle(1:Nions)
                r=rand()
                if s[i]==0 && r<1-exp(-gplus*dt_solve*fplus(sum(s)/Nvesicle)) && sum(s)<Nvesicle
                    s[i]=1
                elseif s[i]==1 && r<1-exp(-gminus*fminus(sum(s)/Nvesicle)*dt_solve)
                    s[i]=0
                end
            end
            if k%Int(dt_save/dt_solve)==1
                Nbound_average[c]+=1/Nsim*sum(s)
                c+=1
            end
        end
    end

    #ode (reaction rate equation)
    Nbound_PDE = zeros(NT_save+1)
    c=2
    unbound=Nions
    for k in 2:NT_solve+1
        bound = Nions-unbound
        unbound = unbound + dt_solve*(gminus*bound*fminus(bound/Nvesicle) - gplus*unbound*fplus(bound/Nvesicle))
        if k%Int(dt_save/dt_solve)==1
            Nbound_PDE[c]+=Nions-unbound
            c+=1
        end
    end

    if save==true
        plot(0:dt_save:T,Nbound_average/Nvesicle,label="particlebased")
        plot!(0:dt_save:T,Nbound_PDE/Nvesicle,label="ODE")
        savefig(string("src/img/oneboxapprox.png"))
    end 
    return Nbound_average,Nbound_PDE
end

function oneboxparam(T=1, Nions = collect(2:1:100), Nvesicles = collect(2:1:50), Nsim = 100)
    final_w_particle = zeros(size(Nions,1),size(Nvesicles,1))
    final_w_pde = zeros(size(Nions,1),size(Nvesicles,1))
    for i in eachindex(Nions)
        for j in eachindex(Nvesicles)
            Nbound_average,Nbound_PDE = compare1box(T,Nions[i],Nvesicles[j],Nsim)
            final_w_particle[i,j] = mean(Nbound_average[end-4:end])/Nvesicles[j]
            final_w_pde[i,j] = mean(Nbound_PDE[end-4:end])/Nvesicles[j] 
        end
    end

    plot()
    heatmap(Nions,Nvesicles,final_w_particle',title="Final occupancy w of particle-based dynamics", c=cmap, clim=(0,1),xlabel="nbr ions",ylabel="nbr vesicle binding sites")
    savefig(string("src/img/oneboxparticles_newfplus.png"))

    plot()
    heatmap(Nions,Nvesicles,final_w_pde',title="Final occupancy w of ODE dynamics", c=cmap, clim=(0,1),xlabel="nbr ions",ylabel="nbr vesicle binding sites")
    savefig(string("src/img/oneboxaode_newfplus.png"))

    plot()
    heatmap(Nions,Nvesicles,(abs.(final_w_particle-final_w_pde)./abs.(final_w_particle))',title="Relative absolute error between final occupancies", c=cmap, clim=(0,1),xlabel="nbr ions",ylabel="nbr vesicle binding sites")
    savefig(string("src/img/oneboxaerror_newfplus.png"))

    return final_w_particle,final_w_pde 
end

function testrates(T=1, N=20, a=0.25, gpluss = [0.1,1], gminuss = [0.1,1], alphas = [0.1, 1], bs = [0.25, 0.75], Nsim= 1000)
    dt_save=0.01
    Nvesicle = a*N
    # uncoop
    p = plot()
    i = 1
    for gplus in gpluss
        for gminus in gminuss
            function fplus(x) # determines binding rate depending on vesicle occupancy
                if x<1
                    return (1-x)
                else
                    return 0
                end
            end
            fminus(x)=1
            q = (; N, gplus, gminus, a, fplus, fminus)

            Nbound_average,Nbound_PDE = compare1box(T,N,Nvesicle,Nsim;q=q,dt_save=dt_save,dt_solve=0.001,save=false)

            plot!(p,0:dt_save:T,Nbound_average/Nvesicle,ylim=(0,1.05),xlim=(0,T+1),label=string("PB, g+ = ",string(gplus), ", g-=",string(gminus)),legend=:bottomright,xlabel ="t",ylabel="w(t)",color=i)
            plot!(p,0:dt_save:T,Nbound_PDE/Nvesicle,ylim=(0,1.05),xlim=(0,T+1),label=string("ODE, g+ = ",string(gplus), ", g-=",string(gminus)),color=i,linestyle=:dash)
            i+=1
        end
    end
    savefig(string("src/img/uncoop.png"))

    # coop unbinding according to "Cooperative stochastic binding and unbinding explain synaptic size dynamics and statistics"
    p = plot()
    i = 1
    for gplus in gpluss
        for gminus in gminuss
            for alpha in alphas
                function fplus(x) # determines binding rate depending on vesicle occupancy
                    if x<1
                        return (1-x)
                    else
                        return 0
                    end
                end
                fminus(x)=1-x+alpha
                q = (; N, gplus, gminus, a, fplus, fminus)

                Nbound_average,Nbound_PDE = compare1box(T,N,Nvesicle,Nsim;q=q,dt_save=dt_save,dt_solve=0.001,save=false)

                plot!(p,0:dt_save:T,Nbound_average/Nvesicle,xlim=(0,T+1),ylim=(0,1.05),label=string("PB, g+ = ",string(gplus), ", g-=",string(gminus), ", a=",string(alpha)),legend=:bottomright,xlabel ="t",ylabel="w(t)",color=i)
                plot!(p,0:dt_save:T,Nbound_PDE/Nvesicle,xlim=(0,T+1),ylim=(0,1.05),label=string("ODE, g+ = ",string(gplus), ", g-=",string(gminus), ", a=",string(alpha)),color=i,linestyle=:dash)
                i+=1
            end
        end
    end
    savefig(string("src/img/coopunbindinglinear.png"))

    # coop binding  according to "Cooperative stochastic binding and unbinding explain synaptic size dynamics and statistics"
    p = plot()
    i = 1
    for gplus in gpluss
        for gminus in gminuss
            for alpha in alphas
                function fplus(x) # determines binding rate depending on vesicle occupancy
                    if x<1
                        return (1-x)*(x + alpha)
                    else
                        return 0
                    end
                end
                fminus(x)=1
                q = (; N, gplus, gminus, a, fplus, fminus)

                Nbound_average,Nbound_PDE = compare1box(T,N,Nvesicle,Nsim;q=q,dt_save=dt_save,dt_solve=0.001,save=false)

                plot!(p,0:dt_save:T,Nbound_average/Nvesicle,xlim=(0,T+1),ylim=(0,1.05),label=string("PB, g+ = ",string(gplus), ", g-=",string(gminus), ", a=",string(alpha)),legend=:bottomright,xlabel ="t",ylabel="w(t)",color=i)
                plot!(p,0:dt_save:T,Nbound_PDE/Nvesicle,xlim=(0,T+1),ylim=(0,1.05),label=string("ODE, g+ = ",string(gplus), ", g-=",string(gminus), ", a=",string(alpha)),color=i,linestyle=:dash)
                i+=1
            end
        end
    end
    savefig(string("src/img/coopbindinglinear.png"))

    # coop unbinding exponenetial
    p = plot()
    i = 1
    for gplus in gpluss
        for gminus in gminuss
            for b in bs
                function fplus(x) # determines binding rate depending on vesicle occupancy
                    if x<1
                        return (1-x) 
                    else
                        return 0
                    end
                end
                fminus(x)=b^(a*N*x-1)
                q = (; N, gplus, gminus, a, fplus, fminus)

                Nbound_average,Nbound_PDE = compare1box(T,N,Nvesicle,Nsim;q=q,dt_save=dt_save,dt_solve=0.001,save=false)

                plot!(p,0:dt_save:T,Nbound_average/Nvesicle,xlim=(0,T+0.5),ylim=(0,1.05),label=string("PB, g+ = ",string(gplus), ", g-=",string(gminus), ", b=",string(b)),legend=:bottomright,xlabel ="t",ylabel="w(t)",color=i)
                plot!(p,0:dt_save:T,Nbound_PDE/Nvesicle,xlim=(0,T+0.5),ylim=(0,1.05),label=string("ODE, g+ = ",string(gplus), ", g-=",string(gminus), ", b=",string(b)),color=i,linestyle=:dash)
                i+=1
            end
        end
    end
    savefig(string("src/img/coopunbindingexp.png"))
end

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



function plotratefunctions(; gpluss = [0.5, 1,2,4], gminus = 2, alphas = [0.1,1], bs = [0.25,0.75])
    legendloc = :outerright

    pl=plot(xlim=(0,1.05),xlabel ="w",ylabel="r+(w)",legend=legendloc)
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
        global w = collect(0:0.01:1)
        global r = [gplus*fplus(j) for j in w]
        plot!(pl,w,r,label=string(L"\gamma^+ = ",string(gplus)),color=i)
        i+=1

    end
    savefig(string("src/img/r_uncoop.png"))
    savefig(string("src/img/r_uncoop.pdf"))

    pl=plot(xlim=(0,1.05),xlabel ="w",ylabel="r-(w)",legend=legendloc)
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
            global w = collect(0:0.01:1)
            global r = [gminus*fminus(j) for j in w]
            plot!(pl,w,r,label=string(L"\gamma^+ = ",string(gplus),L", \alpha^-=",string(alpha)),color=i)


            i+=1
        end

    end
    savefig(string("src/img/r_cooplinearunbinding.png"))
    savefig(string("src/img/r_cooplinearunbinding.pdf"))

    pl=plot(xlim=(0,1.05),xlabel ="w",ylabel="r+(w)",legend=legendloc)
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
            global w = collect(0:0.01:1)
            global r = [gplus*fplus(j) for j in w]
            plot!(pl,w,r,label=string(L"\gamma^+ = ",string(gplus),L", \alpha^-=",string(alpha)),color=i)
            i+=1
        end

    end
    savefig(string("src/img/r_cooplinearbinding.png"))
    savefig(string("src/img/r_cooplinearbinding.pdf"))
    
    (;a) = q
    pl=plot(xlim=(0,1.05),xlabel ="w",ylabel="r-(w)",legend=legendloc)
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
            fminus(x)=b^(5*x-1)
            global w = collect(0:0.01:1)
            global r = [gminus*fminus(j) for j in w]
            plot!(pl,w,r,label=string(L"\gamma^+ = ",string(gplus),L", \beta=",string(b)),color=i)
            i+=1
        end

    end
    savefig(string("src/img/r_coopexpunbinding.png"))
    savefig(string("src/img/r_coopexpunbinding.pdf"))
end
