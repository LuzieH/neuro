function parameters(;
    N = 100, # nbr of calcium ions
    gplus = 5, # rate of ion binding to vesicle
    gminus = 10, # rate of ion unbinding from vesicle
    eps = 0.1785, # radius of interaction ball around vesicle, corresponds to approximately 10% of domain area
    sigma = 1, # noise strength of particles
    sigmav = 0, # noise strength of vesicle
    a = 1/20 # defines vesicle capacity = a*N, needs to be st a*N is an integer
    )

    function fplus(x)
        if x<1
            return 1
        else
            return 0
        end
    end
    
    fminus(x)=1

    q = (; N, gplus,gminus, eps, sigma,sigmav,a,fplus,fminus)
    return q
end


function parametersnewf(;
    N = 100, # nbr of calcium ions
    gplus = 5, # rate of ion binding to vesicle
    gminus = 10, # rate of ion unbinding from vesicle
    eps = 0.1785, # radius of interaction ball around vesicle, corresponds to approximately 10% of domain area
    sigma = 1, # noise strength of particles
    sigmav = 0, # noise strength of vesicle
    a = 1/20, # defines vesicle capacity = a*N, needs to be st a*N is an integer
    b=0.5 #unbinding parameter
    )

    function fplus(x)
        if x<1
            return (1-x)
        else
            return 0
        end
    end

    fminus(x)=b^(a*N*x-1)

    q = (; N, gplus,gminus, eps, sigma,sigmav,a,fplus,fminus)
    return q
end



function PDEconstruct(;
    # Define the constants for the PDE discretization
    dx = 0.025, # dx = dy
    domain = [0 1; 0 1] #only allow square domains
    )
    Nx = Int(round((domain[1,2]-domain[1,1])/dx)) #Nx = Ny, nbr of grid cells per dimension
    dV = dx^2 # grid volume
    X = [x for x in domain[1,1]+0.5*dx:dx:domain[1,2]-0.5*dx, y in domain[2,1]+0.5*dx:dx:domain[2,2]-0.5*dx]
    # to each box x belongs the area [x-0.5*dx, x+0.5*dx] x [y-0.5*dx, y+0.5*dx]
    Y = X'
    gridpoints = [vec(X) vec(X')]
    M = secondderivative(Nx, dx)
    dt = 0.01 # resolution for saving the numerical solution
    p = (; gridpoints, dx, dV, X,Y, Nx, domain,  M,dt)

    return p
end

function particleconstruct(;
    # Define the constants for the particle-dynamics discretization
    dt = 0.0005, #0.0002 works perfect for convergence , dt = 0.0005 also works well
    domain = [0 1; 0 1] #only allow square domains
    )

    p = (; domain, dt)

    return p
end