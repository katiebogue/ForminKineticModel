anonrosen = @(z) SODFHOD(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5)));
anonrosen = @(z) SOS_err(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5)));
opts = optimoptions(@fmincon,"Display","iter",'FunctionTolerance',10^(-10),'OptimalityTolerance',10^(-8),Algorithm="interior-point");
opts = optimoptions(@fmincon,"Display","iter",'FunctionTolerance',10^(-10),'OptimalityTolerance',10^(-8));
opts = optimoptions(@fminunc,"Display","iter",'FunctionTolerance',10^(-10),'OptimalityTolerance',10^(-8));

rng default % For reproducibility
%problem = createOptimProblem("fmincon",...
problem = createOptimProblem("fminunc",...
    x0=rand(1,5),...
    objective=anonrosen,...
    lb=[0;0;0;0;0],...
    options=opts);

gs = GlobalSearch;
gs.MaxTime=120;
gs.NumStageOnePoints=100000;
gs.NumTrialPoints=200000;
[z2,fval2] = run(gs,problem)

ms=MultiStart;
[z2,fval2,exitflag,outpt,solutions] = run(ms,problem,1000)