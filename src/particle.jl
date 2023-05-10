using DifferentialEquations, LinearAlgebra
using Random
using StatsBase
using JLD2
using Plots

"""Produces initial condition for particle dynamics with ion positions uniformly distributed
and all ions unbound, vesicle starting positions are randomly drawn"""
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

"""Produces initial condition for particle dynamics with ion positions uniformly distributed
and all ions unbound, positions of 2 vesicles are specified"""
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
    x0 = [0.1 0.1; 0.4 0.4]*(domain[1,2]-domain[1,1]) .+domain[1,1]
    return y0,s0,x0
end

"""Produces initial condition for particle dynamics with ion positions uniformly distributed
and all ions unbound, position of 1 vesicle is specified"""
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
function particlesolve(NT=100;  p = particleconstruct(), q= parameters(),chosenseed=1)
    Random.seed!(chosenseed)
    (; dt, domain) = p
    (; N, M, sigma,sigmav, eps, a, fplus, fminus, gplus, gminus,initial,force,intforce) = q
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
                    # place ion uniformly in ball of radius eps around x that overlaps with domain
                    # https://mathworld.wolfram.com/DiskPointPicking.html
                    validsample=false
                    while validsample==false 
                        radius = eps*sqrt(rand())
                        theta = rand()*2*pi
                        global sampledpos = x[m,:] + radius* [cos(theta) sin(theta)]' 
                        # check whether sampledpos lies in domain
                        if sampledpos[1]<domain[1,2] && sampledpos[1]>domain[1,1] && sampledpos[2]<domain[2,2] && sampledpos[2]>domain[2,1]
                            validsample=true
                        end
                    end
                    y[i,:] = sampledpos 
                else
                    # update positions
                    y[i,:] = y[i,:] + (s[i]<1)* randn(2)*sqrt(dt)*sigma
                end
            end
        end

        # reflective boundary conditions for Calcium-ions
        indx1 = findall(x->x>domain[1,2],y[:,1])
        indy1 = findall(x->x>domain[2,2],y[:,2])
        indx2 = findall(x->x<domain[1,1],y[:,1])
        indy2 = findall(x->x<domain[2,1],y[:,2])
        y[indx1,1] = - y[indx1,1] .+ 2* domain[1,2] 
        y[indy1,2] = - y[indy1,2] .+ 2* domain[2,2] 
        y[indx2,1] = - y[indx2,1] .+ 2* domain[1,1] 
        y[indy2,2] = - y[indy2,2] .+ 2* domain[2,1] 

        # move vesicles
        for m in 1:M
            x[m,:] = x[m,:] + dt*force(x[m,:]) +randn(2)*sqrt(dt)*sigmav
            if M>1 # interaction force between pairs of vesicles
                x[m,:] += dt* sum([intforce(x[m,:]-x[m2,:]) for m2 in vcat(1:m-1, m+1:M)])
            end
        end

        # reflective boundary conditions for vesicles
        indx1 = findall(x->x>domain[1,2],x[:,1])
        indy1 = findall(x->x>domain[2,2],x[:,2])
        indx2 = findall(x->x<domain[1,1],x[:,1])
        indy2 = findall(x->x<domain[2,1],x[:,2])
        x[indx1,1] = - x[indx1,1] .+ 2* domain[1,2] 
        x[indy1,2] = - x[indy1,2] .+ 2* domain[2,2] 
        x[indx2,1] = - x[indx2,1] .+ 2* domain[1,1] 
        x[indy2,2] = - x[indy2,2] .+ 2* domain[2,1] 

        ys = push!(ys,copy(y))
        ss = push!(ss,copy(s))
        xs = push!(xs,copy(x))
        ws = push!(ws,copy(v/(a*N)))
    end
    return ys, ss, xs, ws
end 

"""solve and plot particle-dynamics"""
function particlesolveplot(T=1; chosenseed=1, p = particleconstruct(), q= parameters())
    (;dt) = p
    NT = Int(round((T/dt)))
    ys, ss, xs, ws  = particlesolve(NT, p=p, q=q,chosenseed=chosenseed)
    particlegif(ys, xs, ss, ws, (p,q),dN=Int(round(NT/100)))
    particleoccupancy(ws,(p,q))
    particlesnapshots(ys, xs, ss, ws, (p,q), collect(1:Int(round(NT/4)):NT+1))
    return ys, ss, xs, ws, (p,q)
end

"""solve ensemble of particle-dynamics"""
function ensemblesolve(T=1, Nsim = 100_000;  p = particleconstruct(), q= parameters(), binnumber = 40, save=true, name="")
    (;dt,domain) = p
    (;N)=q
    NT = Int(round((T/dt)))
    ts = collect(1:Int(round(NT/4)):NT+1)
    nsnapshots = size(ts,1)
    meanhist = zeros(nsnapshots,binnumber,binnumber)
    for s in 1:Nsim
        ys, ss, xs, ws  = particlesolve(NT, p=p, q=q,chosenseed=s)
        for i in eachindex(ts)
            t=ts[i]
            global hist,xrange,yrange = particlehistogram(ys[t],ss[t], domain, binnumber)
            dV= (xrange[2]-xrange[1])* (yrange[2]-yrange[1])
            meanhist[i,:,:]+=1/(Nsim*N*dV)*hist
        end
        if s==1 
            global xsaverage = 1/Nsim*xs
            global wsaverage = 1/Nsim*ws
        else
            xsaverage += 1/Nsim*xs
            wsaverage += 1/Nsim*ws
        end
    end
 
    if save==true
        @save string("data/particleensemble.jld2") meanhist xrange yrange wsaverage xsaverage p q  ts
    end

    return meanhist, wsaverage, xsaverage, xrange, yrange, (p,q), ts 
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


