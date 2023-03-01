using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using JLD2
using Distances
using LaTeXStrings


### TODO use component arrays otherwise get error when reaction happens
## todo set seed 

function second_derivative(Nx, dx)
    M = Tridiagonal(ones(Nx-1), fill(-2., Nx), ones(Nx-1))
    # boundary conditions change in time to ensure flux balance, they are added when solving pde
    M[1,1] = -1.
    M[end,end] = -1.
    M .= (1/dx^2)*M

    return M
end

function centered_difference(Nx, dx)
    # centered first difference for force, doesnt work for different x, y grids
    C = 1/(2*dx)*Tridiagonal(-ones(Nx-1), zeros(Nx), ones(Nx-1))
    # at the boundary a one-sided difference scheme is used
    C[1,1:2] = 1/(dx)* [-1,1]
    C[end,end-1:end] = 1/(dx)* [-1,1]

    return C
end

function ball(radius, center, X)
    Y = X'
    Nx = size(X,1)
    indices = findall([norm([X[i,j] Y[i,j]] - center) for i in 1:Nx, j in 1:Nx].<=radius)
    return indices
end

#2D gaussian that integrates to 1 and centered at center
gaussian(x, center; sigma=0.3) = 1/(2*pi*sigma^2)* exp(-1/(2*sigma^2)*norm(x-center)^2)

function gaussianinit((p,q))
    (; Nx , dV, gridpoints) = p
    (; xv) = q

    # distribution of calcium ions
    c0 = zeros(Nx, Nx)
    center = [-0.5 -0.5]
    c0 += reshape([gaussian(gridpoints[j,:], center') for j in 1:Nx^2], Nx, Nx)
    c0 = c0/(sum(c0)*dV)
    
    # discrete vesicle dynamics, 0: no ion attached, 1: ion attached
    v0 = [0.]
    # vesicle position
    x0 = xv
    return ArrayPartition(c0, v0, x0)
end


function constructinitial((p,q))
    u0 = gaussianinit((p,q))
    return u0
end


function f(du, u,(p,q),t)
    yield()
    (; sigma) = q
    (; dx, Nx,  M) = p
    D = sigma^2 * 0.5
    c, v, x=  u.x
    dc, dv, dx = du.x

    dif = zeros(Nx, Nx)
    a_AB_BAT!(dif, D, M, c)

    dc .= dif 
    dx .= 0

end


# inplace Y = a*(A*B + B+A')
function a_AB_BAT!(Y, a, A, B)
    mul!(Y, A, B)
    mul!(Y, B, A', a, a)
end

function ratebind(u,(p,q),t) 
    c, v, x= u.x
    (; gplus, r, xv) = q
    (; dV, X) = p
    indices = ball(r, xv, X)
    return gplus*(1-v[1])*sum(c[indices])*dV
end



function affectbind!(integrator)
    p,q = integrator.p
    (; r, xv) = q
    (; dV, X) = p
    indices = ball(r, xv, X)
    U= integrator.u
    c, v, x=  U.x
    c[indices] .-=1/(n*dV*size(indices,1))
    v .= 1.
    nothing
end

function rateunbind(u,(p,q),t)
    c, v, x= u.x
    (; gminus) = q

    return gminus*v[1]
end

function affectunbind!(integrator)
    p,q = integrator.p
    (;  r, xv) = q
    (; dV, X) = p
    indices = ball(r, xv, X)
    U = integrator.u
    c, v, x=  U.x
    c[indices] .+=1/(n*dV*size(indices,1))
    v .= 0.
    nothing
end

function PDEsolve(tmax=0.1; alg=Tsit5(), p = PDEconstruct(), q= parameters())

    u0= constructinitial((p,q))
    # Solve the ODE
    prob = ODEProblem(f,u0,(0.0,tmax),(p,q))
    jump = VariableRateJump(ratebind, affectbind!)
    jump2 = VariableRateJump(rateunbind, affectunbind!)
    jump_prob = JumpProblem(prob, Direct(), jump, jump2)
    @time sol = DifferentialEquations.solve(jump_prob, alg,save_start=true)
    return sol, (p,q)
end

function sol2cvx(sol, t)
    solsize = size(sol(t),1)
    c = reshape(sol(t)[1:solsize-5], isqrt(solsize-5), :)
    v = sol(t)[solsize-4]
    x = sol(t)[solsize-3:solsize-2]
    return c,v,x
end

function PDEsolveplot(tmax=0.1; alg=Tsit5(), p = PDEconstruct(), q= parameters())
    sol, (p,q) = PDEsolve(tmax; alg=alg, p = p, q= q)
    PDEgifsingle(sol,(p,q))
    return sol, (p,q)
end
