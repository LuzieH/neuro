using DifferentialEquations, LinearAlgebra
using Random

function particleuniforminit((p,q))
    (; domain) = p
    (;N ) = q
    # distribution of calcium ions
    y0 = rand(N,2) .*(domain[1,2]-domain[1,1]) .+domain[1,1]
    # binding state of calcium ions
    s0 = zeros(N)
    # vesicle position
    x0 = [0.5 0.5]
    return y0,s0,x0
end


function particlesolve(NT=100;  p = particleconstruct(), q= parameters(),chosenseed=1)
    Random.seed!(chosenseed)
    (; dt, domain) = p
    (; N, sigma, eps, a, fplus, fminus, gplus, gminus) = q
    y,s,x = particleuniforminit((p,q))

    ys=[copy(y)]
    ss=[copy(s)]
    xs=[copy(x)]
    ws=[deepcopy(sum(s)/(a*N))]

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

#=             if norm(y[i,:]-x')<=eps &&  s[i]==0
                r = rand()
                if r<1-exp(-gplus*fplus(v/(a*N))*dt)
                    s[i]=1
                    v+=1
                else
                    y[i,:] = y[i,:] + randn(2)*sqrt(dt)*sigma
                end
            elseif s[i]==1
                r = rand()
                if r<1-exp(-gminus*fminus(v/(a*N))*dt)
                    s[i]=0
                    v-=1
                    radius = eps*sqrt(rand())
                    theta = rand()*2*pi
                    y[i,:] = x' + radius* [cos(theta) sin(theta)]'
                end
            else
                y[i,:] = y[i,:] + randn(2)*sqrt(dt)*sigma
            end
 =#
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

function particlesolveplot(NT=100; chosenseed=1, p = particleconstruct(), q= parameters())
    ys, ss, xs, ws  = particlesolve(NT, p=p, q=q,chosenseed=chosenseed)
    particlegifsingle(ys, xs, ss, ws, (p,q),dN=5)
    particleoccupation(ws,(p,q))
    return ys, ss, xs, ws, (p,q)
end

function particlehistogram(y,s,domain, binnumber=20)
    N=size(y,1)
    xrange = range(domain[1,1], domain[1,2],binnumber+1)
    yrange = range(domain[2,1], domain[2,2],binnumber+1)
    hist = zeros(binnumber,binnumber)
    for n in 1:N
        if s[n]==0
            for i in 1:binnumber 
                for j in 1:binnumber
                    if y[n,1]>=xrange[i] && y[n,1]<xrange[i+1] && y[n,2]>=yrange[j] && y[n,2]<yrange[j+1]
                        hist[i,j]+=1
                    end
                end
            end
        end
    end
    #hist=1/N*hist
    return hist,xrange,yrange
end

function comparison(NT=100,Nsim=100,Ns=[1000,10000]; alg=Tsit5(), p1 = PDEconstruct(), p2 = particleconstruct())
    (;domain) = p1
    ws_averages=Any[]
    hist_averages=Any[]

    for k in 1:size(Ns,1)
        N=Ns[k]
        q= parameters(N=N)
        ws_average=Any[]
        hist_average = Any[]
        for n in 1:Nsim
            ys, ss, xs, ws  =particlesolve(NT, p=p2, q=q,chosenseed=n)
            if n==1
                ws_average = 1/Nsim* ws
                hist,xrange,yrange = particlehistogram(ys[end],ss[end],domain, 20) 
                hist_average = 1/Nsim*hist
            elseif n>1
                ws_average+=1/Nsim* ws
                hist,xrange,yrange = particlehistogram(ys[end],ss[end],domain, 20) 
                hist_average += 1/Nsim*hist
            end
        end
        ws_averages=push!(ws_averages,ws_average)
        hist_averages=push!(hist_averages,hist_average)
        (; dt) = p2
    end

    sol, (p1,q) = PDEsolve(NT*dt; alg=alg, p = p1, q= q)
    times = []
    ws_pde=Any[]
    for t in sol.t
        c,w,x = sol2cwx(sol, t)
        times=push!(times,t)
        ws_pde=push!(ws_pde,w[1])
    end

    # plot occupation in time
    subp = plot(times,ws_pde,ylim=(0,1.1),label="occupation PDE",legend=:bottomright)
    for k in 1:size(Ns,1)
        N=Ns[k]
        plot!(range(0,size(ws_averages[k],1)-1)*dt,ws_averages[k],ylim=(0,1.1),label=string("occupation, N = ",string(round(N))))
    end
    savefig(string("src/img/comparison.png"))  

    # plot densities
    for k in 1:size(Ns,1)
        N=Ns[k]
        heatmap(hist_averages[k]',title=string("N = ",string(round(N))))
        savefig(string("src/img/densabm_",string(round(N)),".png"))
    end

    heatmap(c',title="PDE")
    savefig(string("src/img/denspde.png"))
    
end