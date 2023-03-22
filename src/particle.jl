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
    x0 = [0.8 0.8]
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
        yold = y
        y = yold .+ (1 .-s)  .* randn(N,2)*sqrt(dt)*sigma

        # reflective boundary conditions
        indx1 = findall(x->x>domain[1,2],y[:,1])
        indy1 = findall(x->x>domain[2,2],y[:,2])
        indx2 = findall(x->x<domain[1,1],y[:,1])
        indy2 = findall(x->x<domain[2,1],y[:,2])
        y[indx1,1] = - y[indx1,1] .+ 2* domain[1,2] 
        y[indy1,2] = - y[indy1,2] .+ 2* domain[2,2] 
        y[indx2,1] = - y[indx2,1] .+ 2* domain[1,1] 
        y[indy2,2] = - y[indy2,2] .+ 2* domain[2,1] 


        for i in shuffle(1:N)
            r=rand()
            # check binding
            if s[i]==0 && norm(yold[i,:]-x')<=eps && r<1-exp(-gplus*fplus(v/(a*N))*dt)
                s[i]=1
                v+=1 
            # check unbinding
            elseif s[i]==1 && r<1-exp(-gminus*fminus(v/(a*N))*dt)
                s[i]=0
                v-=1 
                # place ion uniformly in ball of radius eps around x
                # https://mathworld.wolfram.com/DiskPointPicking.html
                radius = eps*sqrt(rand())
                theta = rand()*2*pi
                y[i,:] = x' + radius* [cos(theta) sin(theta)]'
            end
        end
        yold = y
        ys = push!(ys,copy(y))
        ss = push!(ss,copy(s))
        xs = push!(xs,copy(x))
        ws = push!(ws,deepcopy(sum(s)/(a*N)))
    end

    return ys, ss, xs, ws
end 

function particlesolveplot(NT=100; alg=Tsit5(), p = particleconstruct(), q= parameters())
    ys, ss, xs, ws  =particlesolve(NT, p=p, q=q)
    particlegifsingle(ys, xs, ws, (p,q),dN=1)
    particleoccupation(ws,(p,q))
    return ys, ss, xs, ws, (p,q)
end

function particlehistogram(y,domain, binnumber=20)
    N=size(y,1)
    xrange = range(domain[1,1], domain[1,2],binnumber+1)
    yrange = range(domain[2,1], domain[2,2],binnumber+1)
    hist = zeros(binnumber,binnumber)
    for n in 1:N
        for i in 1:binnumber 
            for j in 1:binnumber
                if y[n,1]>=xrange[i] && y[n,1]<xrange[i+1] && y[n,2]>=yrange[j] && y[n,2]<yrange[j+1]
                    hist[i,j]+=1
                end
            end
        end
    end
    #hist=1/N*hist
    return hist,xrange,yrange
end

function comparison(NT=100,Nsim=100,Ns=[100,1000,10000,100000]; alg=Tsit5(), p1 = PDEconstruct(), p2 = particleconstruct())
    for k in 1:size(Ns,1)
        N=Ns[k]
        q= parameters(N=N)
        ws_average=Any[]
        for n in 1:Nsim
            ys, ss, xs, ws  =particlesolve(NT, p=p2, q=q,chosenseed=n)
            if n==1
                ws_average = 1/Nsim* ws
            elseif n>1
                ws_average+=1/Nsim* ws
            end
        end
        (; dt) = p2
        if k==1
            subp = plot(range(1,size(ws_average,1))*dt,ws_average,ylim=(0,1.1),label=string("occupation, N = ",string(round(N))))
        elseif k>1
            plot!(range(1,size(ws_average,1))*dt,ws_average,ylim=(0,1.1),label=string("occupation, N = ",string(round(N))))
        end
    end

    sol, (p,q) = PDEsolve(NT*dt; alg=alg, p = p1, q= q)
    times = []
    ws_pde=Any[]
    for t in sol.t
        c,w,x = sol2cwx(sol, t)
        times=push!(times,t)
        ws_pde=push!(ws_pde,w[1])
    end
    plot!(times,ws_pde,ylim=(0,1.1),label="occupation PDE")

    savefig(string("src/img/comparison.png"))  
    
end