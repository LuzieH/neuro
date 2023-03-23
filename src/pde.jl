using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using JLD2
using Distances
using LaTeXStrings

## todo set seed 

function second_derivative(Nx, dx)
    M = Tridiagonal(ones(Nx-1), fill(-2., Nx), ones(Nx-1))
    # produces Neumann bcs
    M[1,1] = -1.
    M[end,end] = -1.
    M .= (1/dx^2)*M

    return M
end

function ball(radius, center, dx, X, delta= 0.000001,Nsamples=10000)
    Y = X'
    Nx = size(X,1)
    # check indices of boxes where one corner lies in ball
    indices1 = findall([norm([X[i,j]-0.5*dx Y[i,j]-0.5*dx] - center) for i in 1:Nx, j in 1:Nx].<=radius+delta)
    indices2 = findall([norm([X[i,j]-0.5*dx Y[i,j]+0.5*dx] - center) for i in 1:Nx, j in 1:Nx].<=radius+delta)
    indices3 = findall([norm([X[i,j]+0.5*dx Y[i,j]-0.5*dx] - center) for i in 1:Nx, j in 1:Nx].<=radius+delta)
    indices4 = findall([norm([X[i,j]+0.5*dx Y[i,j]+0.5*dx] - center) for i in 1:Nx, j in 1:Nx].<=radius+delta)
    # merge indices to get all indices of boxes where at least one corner lies in ball
    allindices = zeros(Nx,Nx)
    allindices[indices1] .=1
    allindices[indices2] .=1
    allindices[indices3] .=1
    allindices[indices4] .=1
    # give proportion of overlapp of ball and box belonging to index
    weights = zeros(Nx,Nx)
    for i in 1:Nx
        for j in 1:Nx
            #count proportion of uniform random samples lying in ball and particular box
            samples = rand(Nsamples,2)*dx .+ [X[i,j]-0.5*dx Y[i,j]-0.5*dx]
            weights[i,j] = 1/Nsamples * size(findall([norm(samples[i,:] - center') for i in 1:Nsamples].<=radius+delta),1)
        end
    end
    return allindices, weights
end


function uniforminit((p,q))
    (; Nx , dV, gridpoints) = p
 
    # distribution of calcium ions
    c0 = ones(Nx, Nx)
    c0 = c0/(sum(c0)*dV)
    
    # vesicle occupation proportion, gives proportion of attached ions
    w0 = [0.]
    # vesicle position
    x0 = [0.5 0.5]
    return ArrayPartition(c0, w0, x0)
end


function constructinitial((p,q))
    u0 = uniforminit((p,q))
    return u0
end


function f(du, u,(p,q,weights),t)
    yield()
    (; sigma, a, fplus, fminus,gplus,gminus) = q
    (; dx, Nx, M, dV) = p
    D = sigma^2 * 0.5
    c, w, x=  u.x
    dc, dw, dx = du.x
    
    # apply diffusion to density c
    dif = zeros(Nx, Nx)
    a_AB_BAT!(dif, D, M, c) # same as dif = D*(M*c + c*M')
    

    # binding and unbinding reactions
    reac = zeros(Nx, Nx)
    #indices = ball(eps, x, X)
    alphafactor = 1/(dV*sum(weights)) # normalization, approx given by 1/area(ball)
    reac = -weights .*c *fplus(w...)*gplus
    reac += weights*fminus(w...)*gminus*a*w[1]*alphafactor 

    dc .= dif .+ reac
    dx .= 0 #vesicle doesnt change location
    dw .= -sum(reac)*(dV/a) 
end


# inplace Y = a*(A*B + B*A')
function a_AB_BAT!(Y, a, A, B)
    mul!(Y, A, B)
    mul!(Y, B, A', a, a)
end
 
  
function PDEsolve(tmax=1.; alg=Tsit5(), p = PDEconstruct(), q= parameters())

    u0= constructinitial((p,q))
    c,w,x = u0.x
    (; eps) = q
    (; X, dx,dt) = p
    allindices, weights = ball(eps, x, dx, X) #as long as vesicle doesn't move
    # Solve the ODE
    prob = ODEProblem(f,u0,(0.0,tmax),(p,q,weights))
    @time sol = DifferentialEquations.solve(prob, alg,save_start=true,saveat = 0:dt:tmax)
    return sol, (p,q)
end

function sol2cwx(sol, t)
    c = sol(t).x[1]
    w = sol(t).x[2]
    x = sol(t).x[3]
    return c,w,x
end

function PDEsolveplot(tmax=1.; alg=Tsit5(), p = PDEconstruct(), q= parameters())
    sol, (p,q) = PDEsolve(tmax; alg=alg, p = p, q= q)
    PDEgifsingle(sol,(p,q))
    PDEoccupation(sol)
    return sol, (p,q)
end
