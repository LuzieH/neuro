using DifferentialEquations, LinearAlgebra
using Random
using StatsBase

"""produces initial conditino for particle dynamics with ion positions uniformly distributed
and all ions unbound"""
function particleuniforminit((p,q))
    (; domain) = p
    (; N, M ) = q
    # distribution of calcium ions
    y0 = rand(N,2) .*(domain[1,2]-domain[1,1]) .+domain[1,1]
    # binding states of calcium ions
    s0 = zeros(N)
    # vesicle position
    x0 = rand(M,2)*(domain[1,2]-domain[1,1]) .+domain[1,1]
    return y0,s0,x0
end

function particleinit2((p,q))
    (; domain) = p
    (; N, M ) = q
    if M!=2
        print("M needs to be 2")
    end
    # distribution of calcium ions
    y0 = rand(N,2) .*(domain[1,2]-domain[1,1]) .+domain[1,1]
    # binding states of calcium ions
    s0 = zeros(N)
    # vesicle position
    x0 = [0.1 0.1; 0.5 0.5]*(domain[1,2]-domain[1,1]) .+domain[1,1]
    return y0,s0,x0
end

function particleinit1((p,q))
    (; domain) = p
    (; N, M ) = q
    if M!=1
        print("M needs to be 1")
    end
    # distribution of calcium ions
    y0 = rand(N,2) .*(domain[1,2]-domain[1,1]) .+domain[1,1]
    # binding states of calcium ions
    s0 = zeros(N)
    # vesicle position
    x0 = [0.5 0.5]*(domain[1,2]-domain[1,1]) .+domain[1,1]
    return y0,s0,x0
end

"""solve the particle-dynamics"""
function particlesolve(NT=100;  p = particleconstruct(), q= parameters(),chosenseed=1, initial="init1")
    Random.seed!(chosenseed)
    (; dt, domain) = p
    (; N, M, sigma, eps, a, fplus, fminus, gplus, gminus) = q
    if initial=="init2"
        y,s,x = particleinit2((p,q))
    elseif initial=="random"
        y,s,x = particleuniforminit((p,q))
    else
        y,s,x = particleinit1((p,q))
    end

    # time series
    ys = [copy(y)]
    ss = [copy(s)]
    xs = [copy(x)]
    v = [size(findall(x->x==i, s),1) for i in 1:M]
    ws = [copy(v/(a*N))]

    for k in 2:NT+1
        # v = [size(findall(x->x==i, s),1) for i in 1:M]
        # binding and unbinding reactions
        for i in shuffle(1:N)
            for m in shuffle(1:M)
                r=rand()
                # if binding happens
                if s[i]==0 && norm(y[i,:]-x[m,:])<=eps && r<1-exp(-gplus*fplus(v[m]/(a*N))*dt)
                    s[i]=m
                    v[m]+=1 
                # unbinding happens
                elseif s[i]==m && r<1-exp(-gminus*fminus(v[m]/(a*N))*dt)
                    s[i]=0
                    v[m]-=1 
                    # place ion uniformly in ball of radius eps around x
                    # https://mathworld.wolfram.com/DiskPointPicking.html
                    radius = eps*sqrt(rand())
                    theta = rand()*2*pi
                    y[i,:] = x[m,:] + radius* [cos(theta) sin(theta)]'
                else
                    # update positions
                    y[i,:] = y[i,:] + (s[i]>0)* randn(2)*sqrt(dt)*sigma
                end
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
        ws = push!(ws,copy(v/(a*N)))
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

    xrange = range(domain[1,1], domain[1,2],binnumber+1)
    yrange = range(domain[2,1], domain[2,2],binnumber+1)
    
    # only if s==0
    indices = findall(s.==0)
    y = y[indices,:]

    h=fit(Histogram,(y[:,1],y[:,2]),(xrange,yrange))

    return h.weights,xrange,yrange
end
