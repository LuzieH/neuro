using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using JLD2
using Distances
using LaTeXStrings

## todo set seed 
## todo multiply binding rate in box by amoung of area of box that overlaps with interaction ball (e.g. by monte carlo, comparing number of corners that lie in circle, or numerically)

function second_derivative(Nx, dx)
    M = Tridiagonal(ones(Nx-1), fill(-2., Nx), ones(Nx-1))
    # produces Neumann bcs
    M[1,1] = -1.
    M[end,end] = -1.
    M .= (1/dx^2)*M

    return M
end

function ball(radius, center, X, delta= 0.000001)
    Y = X'
    Nx = size(X,1)
    indices = findall([norm([X[i,j] Y[i,j]] - center) for i in 1:Nx, j in 1:Nx].<=radius+delta)
    weights # give proportion of overlapp of ball and box belonging to index
    return indices
end

# 2D gaussian that integrates to 1 and centered at center
gaussian(x, center; sigma=0.3) = 1/(2*pi*sigma^2)* exp(-1/(2*sigma^2)*norm(x-center)^2)

function uniforminit((p,q))
    (; Nx , dV, gridpoints) = p
 
    # distribution of calcium ions
    c0 = ones(Nx, Nx)
    center = [-0.5 -0.5]
    #c0 += reshape([gaussian(gridpoints[j,:], center') for j in 1:Nx^2], Nx, Nx)
    c0 = c0/(sum(c0)*dV)
    
    # vesicle occupation proportion, gives proportion of attached ions
    w0 = [0.]
    # vesicle position
    x0 = [0.8 0.8]
    return ArrayPartition(c0, w0, x0)
end


function constructinitial((p,q))
    u0 = uniforminit((p,q))
    return u0
end


function f(du, u,(p,q),t)
    yield()
    (; sigma, eps, a, fplus, fminus,gplus,gminus) = q
    (; dx, Nx, M, X,dV) = p
    D = sigma^2 * 0.5
    c, w, x=  u.x
    dc, dw, dx = du.x
    
    # apply diffusion to density c
    dif = zeros(Nx, Nx)
    a_AB_BAT!(dif, D, M, c)

    # binding and unbinding reactions
    reac = zeros(Nx, Nx)
    indices = ball(eps, x, X)
    alphafactor = 1/(dV*size(indices,1)) # uniform redistribution in ball of radius eps around x
    reac[indices] = -fplus(w...)*gplus*c[indices] 
    reac[indices] .+= fminus(w...)*gminus*a*w*alphafactor 

    dc .= dif .+ reac
    dx .= 0 #vesicle doesnt change location
    dw .= -sum(reac)*dV*(1/a)  

end


# inplace Y = a*(A*B + B+A')
function a_AB_BAT!(Y, a, A, B)
    mul!(Y, A, B)
    mul!(Y, B, A', a, a)
end
 
  
function PDEsolve(tmax=0.1; alg=Tsit5(), p = PDEconstruct(), q= parameters())

    u0= constructinitial((p,q))
    # Solve the ODE
    prob = ODEProblem(f,u0,(0.0,tmax),(p,q))
    @time sol = DifferentialEquations.solve(prob, alg,save_start=true)
    return sol, (p,q)
end

function sol2cwx(sol, t)
    c = sol(t).x[1]
    w = sol(t).x[2]
    x = sol(t).x[3]
    return c,w,x
end

function PDEsolveplot(tmax=0.1; alg=Tsit5(), p = PDEconstruct(), q= parameters())
    sol, (p,q) = PDEsolve(tmax; alg=alg, p = p, q= q)
    PDEgifsingle(sol,(p,q))
    PDEoccupation(sol)
    return sol, (p,q)
end
