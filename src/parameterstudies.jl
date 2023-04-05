using LinearAlgebra



"""compare PDE simulations vs. average (Nsim Monte-Carlo simulations) particle-based simulations with different numbers of ions (given by Ns)"""
function comparison(NT=1000,Nsim=100,Ns=[100,1000,10000]; alg=Tsit5(), q = parameters(), p1 = PDEconstruct(), p2 = particleconstruct(),binnumber = 20,histogram=false)
    (;domain) = p1
    ws_averages=Any[]
    hist_averages=Any[]

    for k in eachindex(Ns)
        N=Ns[k]
        q = merge(q, (;N=N))
        ws_average = zeros(NT+1)
        if histogram==true
            hist_average = Any[]
            hist_average = zeros(binnumber,binnumber)
        end
        for n in 1:Nsim
            ys, ss, xs, ws = particlesolve(NT, p=p2, q=q,chosenseed=n)
            ws_average+=1/Nsim* ws
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
    sol, (p1,q) = PDEsolve(NT*dt; alg=alg, p = p1, q= q)
    times = []
    ws_pde=Any[]
    for t in sol.t
        c,w,x = sol2cwx(sol, t)
        times=push!(times,t)
        ws_pde=push!(ws_pde,w[1])
    end
    if histogram==true
        # plot final concentration of ions in PDE
        c,w,x = sol2cwx(sol, NT*dt)
        heatmap(c',title="PDE")
        savefig(string("src/img/denspde.png"))
    end

    # plot relative occupancy in time
    subp = plot(times,ws_pde,ylim=(0,1.1),label="PDE",legend=:bottomright,xlabel ="t",ylabel="average w(t)")
    for k in eachindex(Ns)
        N=Ns[k]
        plot!(range(0,size(ws_averages[k],1)-1)*dt,ws_averages[k],ylim=(0,1.1),label=string("n = ",string(round(N))))
    end
    savefig(string("src/img/comparison.png"))  
    savefig(string("src/img/comparison.pdf")) 
end


function studydiscretization(T=0.5,Nsim=100,N=10000,dts = [0.001, 0.0005, 0.00025], dxs = [0.01, 0.025, 0.05, 0.1]; q = parameters(N=N), p1 = PDEconstruct(), p2 = particleconstruct())
    (; dt) = p1
    dt_PDE = dt
    NT_PDE = Int(T/dt_PDE)
    lw = 0.5
    subp = plot() 
    # simulate particle based for different dt
    ws_averages=Any[]
    for dt in dts
        p2 = merge(p2, (;dt=dt))
        NT = Int(round(T/dt))
        ws_average = zeros(NT_PDE+1)  
        for n in 1:Nsim        
            ys, ss, xs, ws = particlesolve(NT, p=p2, q=q, chosenseed=n)
            ws_average+=1/Nsim*ws[1:Int(dt_PDE/dt):end]
        end
        ws_averages=push!(ws_averages,ws_average)
 
        plot!(0:dt_PDE:T,ws_average,ylim=(0,1),label=string("particlebased, dt=",string(dt)),legend=:bottomright,xlabel ="t",ylabel="w(t)",linewidth=lw)
    end

    PDE_wss =Any[]
    for dx in dxs
        p1 = merge(p1,(;dx=dx))
        sol, (p1,q) = PDEsolve(T;  p = p1, q= q)
        PDE_ws=Any[]
        for t in 0:dt_PDE:T
            c,w,x = sol2cwx(sol, t)
            PDE_ws=push!(PDE_ws,w[1])
        end
        PDE_wss=push!(PDE_wss,PDE_ws)
        plot!(0:dt_PDE:T,vec(PDE_ws),ylim=(0,1),label=string("PDE, dx=",string(dx)),linewidth=lw)
    end
    savefig(string("src/img/studydiscretization.png"))  
    savefig(string("src/img/studydiscretization.pdf"))  
    return ws_averages, PDE_wss
end

#studyparameters(2,200,100,[0.05, 0.075, 0.1, 0.1785, 0.25], [1, 2.5, 5, 10, 20], [10, 20, 30, 40,50], [1/50, 1/20, 1/10, 1/5]; q=parameters(N=100))

function studyparameters(T=1,Nsim=1000,N=100,eps=[0.05, 0.1, 0.1785, 0.25],gplus=[2.5, 5, 10, 20],gminus=[5, 10, 15, 20], as = [1/50, 1/20, 1/10]; q = parameters(N=N), p1 = PDEconstruct(), p2 = particleconstruct())
    (; dt) = p1
    dt_PDE = dt
    (; dt) = p2
    dt_particle = dt
    NT_PDE = Int(T/dt_PDE)
    NT = Int(round(T/dt_particle))
    
    function ws_compare(q,p1,p2,NT,NT_PDE,Nsim,T,dt_PDE,dt_particle,labelparticle=string("particlebased, eps=",string(e)),labelpde=string("PDE, eps=",string(e)),colorindex=1)
        particles_ws = zeros(NT_PDE+1)  
        #particle based
        for n in 1:Nsim        
            ys, ss, xs, ws = particlesolve(NT, p=p2, q=q, chosenseed=n)
            particles_ws+=1/Nsim*ws[1:Int(dt_PDE/dt_particle):end]
        end
        #pde based
        sol, (p1,q) = PDEsolve(T;  p = p1, q= q)
        PDE_ws=Any[]
        for t in 0:dt_PDE:T
            c,w,x = sol2cwx(sol, t)
            PDE_ws=push!(PDE_ws,w[1])
        end
        plot!(0:dt_PDE:T,particles_ws,ylim=(0,1.05),label=labelparticle,legend=:bottomright,xlabel ="t",ylabel="w(t)",color=colorindex)
        plot!(0:dt_PDE:T,PDE_ws,ylim=(0,1.05),label=labelpde,color=colorindex,linestyle=:dash)
    end


    subp = plot() 
    i=1
    for e in eps
        q1 = merge(q, (;eps=e))
        ws_compare(q1,p1,p2,NT,NT_PDE,Nsim,T,dt_PDE,dt_particle,string("particles, eps=",string(e)),string("PDE, eps=",string(e)),i)
        i+=1
    end
    savefig(string("src/img/studyeps.png")) 
    savefig(string("src/img/studyeps.pdf")) 

    subp = plot() 
    i=1
    for gp in gplus
        q2 = merge(q, (;gplus=gp))
        ws_compare(q2,p1,p2,NT,NT_PDE,Nsim,T,dt_PDE,dt_particle,string("particles, gamma+ =",string(gp)),string("PDE, gamma+ =",string(gp)),i)
        i+=1
    end
    savefig(string("src/img/studygammaplus.png")) 
    savefig(string("src/img/studygammaplus.pdf")) 

    subp = plot() 
    i=1
    for gm in gminus
        q3 = merge(q, (;gminus=gm))
        ws_compare(q3,p1,p2,NT,NT_PDE,Nsim,T,dt_PDE,dt_particle,string("particles, gamma- =",string(gm)),string("PDE, gamma- =",string(gm)),i)
        i+=1
    end
    savefig(string("src/img/studygammaminus.png")) 
    savefig(string("src/img/studygammaminus.pdf")) 

    subp = plot() 
    i=1
    for a in as
        q4 = merge(q, (;a=a))
        ws_compare(q4,p1,p2,NT,NT_PDE,Nsim,T,dt_PDE,dt_particle,string("particles, a =",string(a)),string("PDE, a =",string(a)),i)
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

    #pde
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