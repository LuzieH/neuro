using DifferentialEquations, LinearAlgebra
using Random

"""produces initial conditino for particle dynamics with ion positions uniformly distributed
and all ions unbound"""
function particleuniforminit((p,q))
    (; domain) = p
    (; N ) = q
    # distribution of calcium ions
    y0 = rand(N,2) .*(domain[1,2]-domain[1,1]) .+domain[1,1]
    # binding states of calcium ions
    s0 = zeros(N)
    # vesicle position
    x0 = [0.5 0.5]
    return y0,s0,x0
end

"""solve the particle-dynamics"""
function particlesolve(NT=100;  p = particleconstruct(), q= parameters(),chosenseed=1)
    Random.seed!(chosenseed)
    (; dt, domain) = p
    (; N, sigma, eps, a, fplus, fminus, gplus, gminus) = q
    y,s,x = particleuniforminit((p,q))

    # time series
    ys=[copy(y)]
    ss=[copy(s)]
    xs=[copy(x)]
    ws=[copy(sum(s)/(a*N))]

    for k in 2:NT+1
        v = sum(s)
        # binding and unbinding reactions
        for i in shuffle(1:N)
            r=rand()
            # if binding happens
            if s[i]==0 && norm(y[i,:]-x')<=eps && r<1-exp(-gplus*fplus(v/(a*N))*dt)
                s[i]=1
                v+=1 
            # unbinding happens
            elseif s[i]==1 && r<1-exp(-gminus*fminus(v/(a*N))*dt)
                s[i]=0
                v-=1 
                # place ion uniformly in ball of radius eps around x
                # https://mathworld.wolfram.com/DiskPointPicking.html
                radius = eps*sqrt(rand())
                theta = rand()*2*pi
                y[i,:] = x' + radius* [cos(theta) sin(theta)]'
            else
                # update positions
                y[i,:] = y[i,:] + (1 -s[i])* randn(2)*sqrt(dt)*sigma
            end
        end

        # reflective boundary conditions
        indx1 = findall(x->x>domain[1,2],y[:,1])
        indy1 = findall(x->x>domain[2,2],y[:,2])
        indx2 = findall(x->x<domain[1,1],y[:,1])
        indy2 = findall(x->x<domain[2,1],y[:,2])
        y[indx1,1] = - y[indx1,1] .+ 2* domain[1,2] 
        y[indy1,2] = - y[indy1,2] .+ 2* domain[2,2] 
        y[indx2,1] = - y[indx2,1] .+ 2* domain[1,1] 
        y[indy2,2] = - y[indy2,2] .+ 2* domain[2,1] 

        ys = push!(ys,copy(y))
        ss = push!(ss,copy(s))
        xs = push!(xs,copy(x))
        ws = push!(ws,deepcopy(sum(s)/(a*N)))
    end
    return ys, ss, xs, ws
end 

"""solve and plot particle-dynamics"""
function particlesolveplot(NT=100; chosenseed=1, p = particleconstruct(), q= parameters())
    ys, ss, xs, ws  = particlesolve(NT, p=p, q=q,chosenseed=chosenseed)
    particlegifsingle(ys, xs, ss, ws, (p,q),dN=5)
    particleoccupancy(ws,(p,q))
    return ys, ss, xs, ws, (p,q)
end

"""produce histogram from unbound ion positions"""
function particlehistogram(y,s,domain, binnumber=20)
    N=size(y,1)
    xrange = range(domain[1,1], domain[1,2],binnumber+1)
    yrange = range(domain[2,1], domain[2,2],binnumber+1)
    hist = zeros(binnumber,binnumber)
    for n in 1:N
        if s[n]==0 # if unbound
            for i in 1:binnumber 
                for j in 1:binnumber
                    if y[n,1]>=xrange[i] && y[n,1]<xrange[i+1] && y[n,2]>=yrange[j] && y[n,2]<yrange[j+1]
                        hist[i,j]+=1
                    end
                end
            end
        end
    end
    return hist,xrange,yrange
end

"""compare PDE simulations vs. average (Nsim Monte-Carlo simulations) particle-based simulations with different numbers of ions (given by Ns)"""
function comparison(NT=100,Nsim=100,Ns=[100,1000,10000,100000]; alg=Tsit5(), p1 = PDEconstruct(), p2 = particleconstruct(),binnumber = 20)
    (;domain) = p1
    ws_averages=Any[]
    hist_averages=Any[]

    for k in eachindex(Ns)
        N=Ns[k]
        q= parameters(N=N)
        ws_average=Any[]
        hist_average = Any[]
        ws_average = zeros(NT+1)
        hist_average = zeros(binnumber,binnumber)
        for n in 1:Nsim
            ys, ss, xs, ws = particlesolve(NT, p=p2, q=q,chosenseed=n)
            ws_average+=1/Nsim* ws
            hist,xrange,yrange = particlehistogram(ys[end],ss[end],domain, binnumber) 
            hist_average += 1/Nsim*hist
        end
        ws_averages=push!(ws_averages,ws_average)
        hist_averages=push!(hist_averages,hist_average)

        # plot average final concentration of particle-dynamics
        heatmap(hist_average',title=string("N = ",string(round(N))))
        savefig(string("src/img/densabm_",string(round(N)),".png"))
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
    # plot final concentration of ions in PDE
    c,w,x = sol2cwx(sol, NT*dt)
    heatmap(c',title="PDE")
    savefig(string("src/img/denspde.png"))

    # plot relative occupancy in time
    subp = plot(times,ws_pde,ylim=(0,1.1),label="PDE",legend=:bottomright)
    plt.xlabel("t")
    plt.ylabel("average w(t)")
    for k in eachindex(Ns)
        N=Ns[k]
        plot!(range(0,size(ws_averages[k],1)-1)*dt,ws_averages[k],ylim=(0,1.1),label=string("n = ",string(round(N))))
    end
    savefig(string("src/img/comparison.png"))  
end