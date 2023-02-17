using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using JLD2
using Distances
using LaTeXStrings


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

#2D gaussian that integrates to 1 and centered at center
gaussian(x, center; sigma=0.1) = 1/(2*pi*sigma^2)* exp(-1/(2*sigma^2)*norm(x-center)^2)

function gaussianinit((p,q))
    (; Nx , dV, gridpoints) = p
    (; xv) = q

    # distribution of calcium ions
    c0 = zeros(Nx, Ny)
    center = [-0.5 -0.5]
    c0 += reshape([gaussian(gridpoints[j,:], center) for j in 1:Nx^2], Nx, Nx)
    c0 = c0/(sum(c0)*dV)
    
    # discrete vesicle dynamics, 0: no ion attached, 1: ion attached
    v0 = 0
    # vesicle position
    x0 = xv
    return ArrayPartition(c0, v0, x0)
end


function constructinitial((p,q))
    X0 = gaussianinit((p,q))
    return X0
end


function f(dX, X,(p,q),t)
    yield()
    (;  N, gplus,gminus, r, xv, sigma) = q
    (; gridpoints, dx, dV, X, Nx, domain, C,  M) = p
    
    D = sigma^2 * 0.5
    c, v, x= X.x
    dc, dv, dx = dX.x

    dif  = zeros(Nx, Ny)
    a_AB_BAT!(dif, D, M, c)

    #balance fluxes at boundary (part of boundary conditions)
#=     dif[1,:]+= -D/dx * (cforce[1,:,1])
    dif[end,:]+= D/dx * (cforce[end,:,1])
    dif[:,1]+= -D/dy * (cforce[:,1,2])
    dif[:,end]+= D/dy * (cforce[:,end,2]) =#

    dc .=  dif 
    dx .= 0

end


# inplace Y = a*(A*B + B+A')
function a_AB_BAT!(Y, a, A, B)
    mul!(Y, A, B)
    mul!(Y, B, A', a, a)
end

function ratebind(X,(p,q),t) 
    c, v, x= X.x
    (;  N, gplus,gminus, r, xv, sigma) = q
    (; gridpoints, dx, dV, X, Nx, domain, C,  M) = p

    return gplus*(1-v)*c ##correctly adapt, replace c by integral of c over ball
end

function rateunbind(X,(p,q),t)
    c, v, x= X.x
    (; gminus) = q

    return gminus*v 
end


function affectbind!(integrator)
    integrator.u[1] = integrator.u[1] / 2
    nothing
end

function affectunbind!(integrator)
    integrator.u[1] = integrator.u[1] / 2
    nothing
end

function PDEsolve(tmax=0.1; alg=nothing, init="4inf", p = PDEconstruct(), q= parameters())

    X0= constructinitial((p,q))
    
    # Solve the ODE
    prob = ODEProblem(f,X0,(0.0,tmax),(p,q))
    jump = VariableRateJump(ratebind, affectbind!)
    jump2 = VariableRateJump(rateunbind, affectunbind!)
    jump_prob = JumpProblem(prob, Direct(), jump, jump2)
    @time sol = DifferentialEquations.solve(jump_prob, alg,  save_start=true)

    return sol, (p,q)
end

function sol2uyz(sol, t)
    u = sol(t).x[1]
    z = sol(t).x[2]
    y = sol(t).x[3]
    return u,z,y
end
