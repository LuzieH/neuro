function parameters(;
    N = 100, #nbr of calcium ions
    gplus = 10,
    gminus = 2, 
    r = 0.1, #radius of interaction ball around vesicle
    xv = [0 0.5], #vesicle position
    sigma = 1
    )

    q = (; N, gplus,gminus, r, xv, sigma)
    return q
end


function PDEconstruct(;
    # Define the constants for the PDE discretization
    dx = 0.05, # dx = dy
    domain = [-1 1; -1 1] #only allow square domains
    )
    Nx = Int(round((domain[1,2]-domain[1,1])/dx+1)) #Nx = Ny, nbr of grid cells per dimension
    dV = dx^2 # grid volume
    X = [x for x in domain[1,1]:dx:domain[1,2], y in domain[2,1]:dx:domain[2,2]]
    gridpoints = [vec(X) vec(X)]
    M = second_derivative(Nx, dx)
    C = centered_difference(Nx, dx)

    p = (; gridpoints, dx, dV, X, Nx, domain, C,  M)

    return p
end