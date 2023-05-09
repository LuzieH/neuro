using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using JLD2
using Distances
using LaTeXStrings

 
"""computes second derivative matrix for finite-difference discretization"""
function secondderivative(Nx, dx)
    Mmatrix = Tridiagonal(ones(Nx-1), fill(-2., Nx), ones(Nx-1))
    # produces Neumann bcs
    Mmatrix[1,1] = -1.
    Mmatrix[end,end] = -1.
    Mmatrix .= (1/dx^2)*Mmatrix
    return Mmatrix
end


 
"""computes overlapp of all grid cells with ball around vesicle with Monte Carlo. This function is a bottleneck for solving 
the PDE and becomes very expensive when Nsamples is large"""
function ballslow(radius, center, dx, X, delta= 0.000001,Nsamples=10)
    Y = X'
    Nx = size(X,1)
    # give proportion of overlapp of ball and box belonging to certain index
    weights = zeros(Nx,Nx)
    for i in 1:Nx
        for j in 1:Nx
            #count proportion of uniform random samples lying in ball and particular box
            samples = rand(Nsamples,2)*dx .+ [X[i,j]-0.5*dx Y[i,j]-0.5*dx]
            weights[i,j] = 1/Nsamples * size(findall([norm(samples[i,:] - center) for i in 1:Nsamples].<=radius+delta),1)
        end
    end
    return weights
end

"""computes overlapp of all grid cells with ball around vesicle via a regular grid."""
function ball(radius, center, dx, X, delta= 0.000001,Nsamplesperdim=3)
    Y = X'
    Nx = size(X,1)
    Nsamples = Nsamplesperdim^2
    # give proportion of overlapp of ball and box belonging to certain index
    weights = zeros(Nx,Nx)
    stepsize= 1/Nsamplesperdim
    for i in 1:Nx
        for j in 1:Nx
            #count proportion of uniform random samples lying in ball and particular box
            grid = [x for x in 0.5*stepsize:stepsize:1-0.5*stepsize, y in 0.5*stepsize:stepsize:1-0.5*stepsize]
            unitsamples = [vec(grid) vec(grid')]
            samples = unitsamples*dx .+ [X[i,j]-0.5*dx Y[i,j]-0.5*dx]
            weights[i,j] = 1/Nsamples * size(findall([norm(samples[i,:] - center) for i in 1:Nsamples].<=radius+delta),1)
        end
    end
    return weights
end

"""produces initial condition for the PDE with uniform distribution of ions,
all ions still unbound, and random initial positions of vesicles"""
function uniforminit((p,q))
    (; Nx, dV, domain) = p
    (; M ) = q
    # distribution of calcium ions
    c0 = ones(Nx, Nx)
    c0 = c0/(sum(c0)*dV)

    # vesicle occupation proportion, gives proportion of attached ions
    w0 = zeros(M)
    # vesicle position
    x0 = rand(M,2)*(domain[1,2]-domain[1,1]) .+domain[1,1]
    return ArrayPartition(c0, w0, x0)
end

"""produces initial condition for the PDE with uniform distribution of ions,
all ions still unbound, and initial condition of 2 vesicles specified """
function init2((p,q))
    (; Nx, dV, domain) = p
    (; M ) = q
    if M!=2
        print("M needs to be 2")
    end
    # distribution of calcium ions
    c0 = ones(Nx, Nx)
    c0 = c0/(sum(c0)*dV)

    # vesicle occupation proportion, gives proportion of attached ions
    w0 = zeros(2)
    # vesicle position
    x0 = [0.1 0.1; 0.4 0.4]*(domain[1,2]-domain[1,1]) .+domain[1,1]
    return ArrayPartition(c0, w0, x0)
end

"""produces initial condition for the PDE with uniform distribution of ions,
all ions still unbound, and initial condition of 1 vesicle specified """
function init1((p,q))
    (; Nx, dV, domain) = p
    (; M ) = q
    if M!=1
        print("M needs to be 1")
    end
    # distribution of calcium ions
    c0 = ones(Nx, Nx)
    c0 = c0/(sum(c0)*dV)

    # vesicle occupation proportion, gives proportion of attached ions
    w0 = zeros(1)
    # vesicle position
    x0 = [0.5 0.5]*(domain[1,2]-domain[1,1]) .+domain[1,1]
    return ArrayPartition(c0, w0, x0)
end

function constructinitial((p,q); initial = "init1")
    if initial=="init2"
        u0 = init2((p,q))
    elseif initial == "random"
        u0 = uniforminit((p,q))
    elseif initial =="init1"
        u0 = init1((p,q))
    end    
    return u0
end

"""deterministic RHS of coupled PDE-dynamics"""
function f(du, u,(p,q),t)
    yield()
    (; sigma,  M, a, fplus, fminus,gplus,gminus,eps,force,intforce) = q
    (; dx, Mmatrix, Nx, dV,X) = p
    dx_grid = dx
    D = sigma^2 * 0.5
    c, w, x=  u.x
    dc, dw, dx = du.x
    
    # apply diffusion to density c
    dif = zeros(Nx, Nx)
    a_AB_BAT!(dif, D, Mmatrix, c) # same as dif = D*(M*c + c*M')
    
    # binding and unbinding reactions
    reacall = zeros(Nx, Nx)
    
    # loop over all vesicles
    for m in 1:M
        reacm = zeros(Nx, Nx)
        wghts = ball(eps, x[m,:], dx_grid, X)  
        alphafactor = 1/(dV*sum(wghts)) # normalization, approx given by 1/area(ball)
        reacm -= wghts .*c *fplus(w[m])*gplus
        reacm += wghts*fminus(w[m])*gminus*a*w[m]*alphafactor
        reacall += reacm
        # changes of relative vesicle occupancy
        dw[m] = -sum(reacm)*(dV/a) 
    end

    # total changes on concentration
    dc .= dif .+ reacall

    # vesicle movement
    for m in 1:M
        dx[m,:] = force(x[m])
        if M>1
            dx[m,:]+=sum([intforce(x[m,:]-x[m2,:]) for m2 in vcat(1:m-1, m+1:M)])
        end
    end
end

"""stochastic RHS of coupled PDE-dynamics"""
function g(du, u,(p,q),t)
    yield()
    (; sigmav) = q
    dc, dw, dx = du.x
    
    # no noise on concentration and vesicle occupancy
    dw .= 0
    dc .= 0

    # vesicle movement noise
    dx .= sigmav 
    
end


"""inplace Y = a*(A*B + B*A')"""
function a_AB_BAT!(Y, a, A, B)
    mul!(Y, A, B)
    mul!(Y, B, A', a, a)
end

"""Ensures reflecting boundary conditions of vesicle"""
function affect!(integrator)
    p,q = integrator.p
    (;domain) = p
    (;M)=q
    for m in 1:M
        if integrator.u.x[3][m,1] < domain[1,1] 
            integrator.u.x[3][m,1] = 2*domain[1,1] - integrator.u.x[3][m,1]
        elseif integrator.u.x[3][m,1] > domain[1,2]
            integrator.u.x[3][m,1] = 2*domain[1,2] - integrator.u.x[3][m,1]
        elseif integrator.u.x[3][m,2] < domain[2,1]
            integrator.u.x[3][m,2] = 2*domain[2,1] - integrator.u.x[3][m,2]
        elseif integrator.u.x[3][m,2] > domain[2,2]
            integrator.u.x[3][m,2] = 2*domain[2,2] - integrator.u.x[3][m,2]
        end
    end
end

"""Solves the PDE"""  
function PDEsolve(tmax=1.; p = PDEconstruct(), q= parameters())
    (;initial,sigmav,M) = q
    u0= constructinitial((p,q);initial=initial)
    (; dt,domain) = p
    if sigmav>0
        prob = SDEProblem(f,g, u0,(0.0,tmax),(p,q))
    else
        prob = ODEProblem(f,u0,(0.0,tmax),(p,q))
    end
    # boundary conditions for vesicle movement using callback function
    condition(u,t,integrator) = true
    cb = DiscreteCallback(condition,affect!;save_positions=(false,true))
    #condition(u,t,integrator) = prod([(u.x[3][m,1]-domain[1,1])*(u.x[3][m,1]-domain[1,2])*(u.x[3][m,2]-domain[2,1])*(u.x[3][m,2]-domain[2,2]) for m in 1:M]) #product of distances to boundary
    #cb = ContinuousCallback(condition,affect!;save_positions=(false,true))
    # Solve the PDE coupled to vesicle-SDEs
    sol = DifferentialEquations.solve(prob, callback = cb, save_start=true,saveat = 0:dt:tmax,dtmax = 0.05) #alg = alg
    return sol, (p,q)
end

"""Convert PDE solution to density, relative occupancy and vesicle position 
for a given time point"""
function sol2cwx(sol, t)
    c = sol(t).x[1]
    w = sol(t).x[2]
    x = sol(t).x[3]
    return c,w,x
end

"""Solve and plot PDE solution"""
function PDEsolveplot(tmax=1.;  p = PDEconstruct(), q= parameters(),clim=(0,2))
    sol, (p,q) = PDEsolve(tmax;  p = p, q= q)
    PDEgif(sol,(p,q),dt = tmax/50)
    PDEoccupancy(sol,(p,q))
    PDEsnapshots(sol, (p,q), collect(0:(tmax/4):tmax); save = true, name="", clim=clim)
    return sol, (p,q)
end
