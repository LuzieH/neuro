function parameters(;
    N = 100, # nbr of calcium ions
    M =  2, # number of vesicles
    gplus = 4, # rate of ion binding to vesicle
    gminus = 2, # rate of ion unbinding from vesicle
    eps = 0.2, # radius of interaction ball around vesicle, corresponds to approximately 10% of domain area
    sigma = 0.25, # noise strength of particles
    sigmav = 0., # noise strength of vesicle
    a = 1/20, # defines vesicle capacity = a*N, needs to be st a*N is an integer
    initial = "init2", # defines the initial conditions
    aforce = 0.25,
    aint = 0.05,
    bint = 5
    )

    function fplus(x) # determines binding rate depending on vesicle occupancy
        if x<1
            return (1-x)
        else
            return 0
        end
    end
    
    fminus(x)=1 # determines unbinding rate depending on vesicle occupancy

    # force for vesicle movement
    force(x) = aforce*[0, -1] # downward force

    # interatomic force between vesicles depending on (x_i-x_j) vector
    intforce(dist) = aint*bint*exp(-bint*norm(dist))*dist/norm(dist) # shortrange repulsion

    q = (; N, M, gplus, gminus, eps, sigma, sigmav, a, fplus, fminus,initial, force,intforce)
    return q
end

function PDEconstruct(;
    # Define the constants for the PDE discretization
    dx = 0.025, # dx = 0.05 is much faster! and also ok
    domain = [0 1; 0 1] #only allow square domains
    )
    Nx = Int(round((domain[1,2]-domain[1,1])/dx)) #Nx = Ny, nbr of grid cells per dimension
    dV = dx^2 # grid volume
    X = [x for x in domain[1,1]+0.5*dx:dx:domain[1,2]-0.5*dx, y in domain[2,1]+0.5*dx:dx:domain[2,2]-0.5*dx]
    # to each box x belongs the area [x-0.5*dx, x+0.5*dx] x [y-0.5*dx, y+0.5*dx]
    Y = X'
    gridpoints = [vec(X) vec(X')]
    Mmatrix = secondderivative(Nx, dx)
    dt = 0.01 # resolution for saving the numerical solution
    p = (; gridpoints, dx, dV, X,Y, Nx, domain,  Mmatrix,dt)
    return p
end

function particleconstruct(;
    # Define the constants for the particle-dynamics discretization
    dt = 0.0002, #0.0002 works perfect for convergence , dt = 0.0005 also works well
    domain = [0 1; 0 1] #only allow square domains
    )
    p = (; domain, dt)
    return p
end