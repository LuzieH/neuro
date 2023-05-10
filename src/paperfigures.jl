# functions to make figures for paper

function papersimulations(;T=1, Nsim = 100_000, p1 = particleconstruct(), p2 = PDEconstruct(), q = parameters())
    particlesolveplot(T; chosenseed=1, p = p1, q = q);
    @load string("data/particleensemble.jld2") meanhist xrange yrange wsaverage xsaverage p q  ts
    ensembleplot(meanhist, wsaverage, xsaverage, xrange, yrange, ts, (p,q); clim=(minimum(meanhist)-0.1,maximum(meanhist)+0.1));
    PDEsolveplot(T;  p = p2, q = q,clim=(minimum(meanhist)-0.1,maximum(meanhist)+0.1));
end

function paperratecomparison(;T=1, Ns = [100,1_000], gpluss = [0.5, 1,2,4], gminus = 2, alphas = [0.1,1], bs = [0.25,0.75], Nsims = [5_000,500], load=true, save=false)
    comparerates(T=T, Ns = Ns, gpluss =gpluss, gminus = gminus, alphas = alphas, bs = bs, Nsims= Nsims, load=load, save=save);
end
