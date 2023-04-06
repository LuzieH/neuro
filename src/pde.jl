using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using JLD2
using Distances
using LaTeXStrings

# TODOs:
# include several vesicles
# include vesicle dynamics

"""computes second derivative matrix for finite-difference discretization"""
function secondderivative(Nx, dx)
    Mmatrix = Tridiagonal(ones(Nx-1), fill(-2., Nx), ones(Nx-1))
    # produces Neumann bcs
    Mmatrix[1,1] = -1.
    Mmatrix[end,end] = -1.
    Mmatrix .= (1/dx^2)*Mmatrix
    return Mmatrix
end

"""computes overlapp of all grid cells with ball around vesicle"""
function ball(radius, center, dx, X, delta= 0.000001,Nsamples=100)
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

"""produces initial condition for the PDE with uniform distribution of ions,
and all ions still unbound"""
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
    x0 = [0.1 0.1; 0.5 0.5]*(domain[1,2]-domain[1,1]) .+domain[1,1]
    return ArrayPartition(c0, w0, x0)
end

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
    w0 = zeros(M)
    # vesicle position
    x0 = [0.5 0.5]*(domain[1,2]-domain[1,1]) .+domain[1,1]
    return ArrayPartition(c0, w0, x0)
end

function constructinitial((p,q); initial = "init1")
    if initial=="init2"
        u0 = init2((p,q))
    elseif initial == "random"
        u0 = uniforminit((p,q))
    else
        u0 = init1((p,q))
    end    
    
    return u0
end

"""RHS of coupled PDE-dynamics"""
function f(du, u,(p,q),t)
    yield()
    (; sigma, M, a, fplus, fminus,gplus,gminus,eps) = q
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
        weights = ball(eps, x[m,:], dx_grid, X)  
        alphafactor = 1/(dV*sum(weights)) # normalization, approx given by 1/area(ball)
        reacm -= weights .*c *fplus(w[m])*gplus
        reacm += weights*fminus(w[m])*gminus*a*w[m]*alphafactor
        reacall += reacm
        # changes of relative vesicle occupancy
        dw[m] = -sum(reacm)*(dV/a) 
    end

    # total changes on concentration
    dc .= dif .+ reacall

    # vesicle movement
    dx .= 0 #sigmav*randn(2)*sqrt(dt)
end


"""inplace Y = a*(A*B + B*A')"""
function a_AB_BAT!(Y, a, A, B)
    mul!(Y, A, B)
    mul!(Y, B, A', a, a)
end
 
"""solve the PDE"""  
function PDEsolve(tmax=1.; alg=Tsit5(), p = PDEconstruct(), q= parameters(), initial = "init1")

    u0= constructinitial((p,q);initial=initial)
    (; dt) = p
    
    # Solve the PDE
    prob = ODEProblem(f,u0,(0.0,tmax),(p,q))
    sol = DifferentialEquations.solve(prob, alg,save_start=true,saveat = 0:dt:tmax)
    return sol, (p,q)
end

"""convert PDE solution to density, relative occupancy and vesicle position 
for a given time point"""
function sol2cwx(sol, t)
    c = sol(t).x[1]
    w = sol(t).x[2]
    x = sol(t).x[3]
    return c,w,x
end

"""solve and plot PDE solution"""
function PDEsolveplot(tmax=1.; alg=Tsit5(), p = PDEconstruct(), q= parameters(),initial="init1")
    sol, (p,q) = PDEsolve(tmax; alg=alg, p = p, q= q, initial=initial)
    PDEgifsingle(sol,(p,q))
    PDEoccupancy(sol,(p,q))
    return sol, (p,q)
end
