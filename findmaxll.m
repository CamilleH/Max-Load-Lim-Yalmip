function results = findmaxll(mpc,dir)
% FINDMAXLL finds the maximum loadability limit point in a particular
% direction by using the YALMIP implementation that is formulated with the
% voltages expressed in rectangular coordinates.
%    DIR is the direction of load increase
%% Open the system file and creates the necessary matrices
define_constants;
n = size(mpc.bus,1);
idxVarPQ = find(dir);

%% Scale to per-unit
mpc0 = mpc;
tanphi0 = mpc.bus(idxVarPQ,QD)./mpc.bus(idxVarPQ,PD);
mpc.bus(:,[PD,QD]) = mpc.bus(:,[PD,QD])/mpc.baseMVA;
mpc.gen(:,[PG,QG]) = mpc.gen(:,[PG,QG])/mpc.baseMVA;
mpc.gen(:,[QMAX,QMIN]) = mpc.gen(:,[QMAX,QMIN])/mpc.baseMVA;
if isfield(mpc,'wind') && ~isempty(mpc.wind)
    mpc.wind(:,[PG,QG]) = mpc.wind(:,[PG,QG])/mpc.baseMVA;
    mpc.wind(:,4) = mpc.wind(:,4)/mpc.baseMVA; % wind farm capacity
end

%% Optimization problem
% Parameters
dirP = dir;

% optimization problem
x = sdpvar(2*n,1);
Ploads = sdpvar(n,1);
assign(x,ones(2*n,1));
assign(Ploads,mpc.bus(:,PD));
lambda = sdpvar(1,1);
Objective = -lambda;
Constraints = fmlp_buildConstraints(mpc,x,Ploads,[],idxVarPQ,[],[]);
Constraints = [Constraints,...
    Ploads == mpc.bus(:,PD)+lambda*dirP];
options = sdpsettings('verbose',0,'showprogress',1,'fmincon.Algorithm',...
    'interior-point','usex0',1,'savesolveroutput',1);
sol = optimize(Constraints,Objective,options);

% Extracting the result, scaling to MW
Ploads_val = value(Ploads);
x_val = value(x);
V_res = x_val(1:n)+1i*x_val(n+1:2*n);
results = mpc0; 
results.bus(idxVarPQ,PD) = Ploads_val(idxVarPQ)*mpc.baseMVA;
results.bus(idxVarPQ,QD) = tanphi0.*Ploads_val(idxVarPQ)*mpc.baseMVA;
results.bus(:,VM) = abs(V_res);
results.bus(:,VA) = angle(V_res)*180/pi;
results.stab_marg = value(lambda);
results.sol = sol;
% Extracting the multipliers
[~, pv, pq] = bustypes(results.bus, results.gen);
npv = length(pv);
npq = length(pq);
npvpq = npv+npq;
ngen = size(results.gen,1);
genpv = ismember(results.gen(:,GEN_BUS),pv);
ngenpv = sum(genpv);
results.nln.mu.u.Pmis = zeros(2*n,1);
results.nln.mu.u.Qmis = zeros(2*n,1);
results.var.mu.u.Qg = zeros(ngen,1);
results.var.mu.u.Vm = zeros(n,1);
results.nln.mu.u.Pmis([pv;pq]) = sol.solveroutput.lambda.eqnonlin(1:npvpq);
results.nln.mu.u.Qmis(pq) = sol.solveroutput.lambda.eqnonlin(npvpq+(1:npq));
results.var.mu.u.Qg(genpv) = sol.solveroutput.lambda.ineqnonlin(1:ngenpv);
results.var.mu.u.Qg(results.var.mu.u.Qg<1e-5) = 0; % numerical approximations can lead to close to zero Lag. multipliers
% Note: below, we take the Lagrange multiplier lambda_2 of the limit V_g^2 and
% transform it to the Lagrange multiplier lambda_1 of the limit V_g. The Lag. multi.
% are such that lambda_1 = - df^* / dV_g and lambda_2 = -df/d(V_g)^2. From
% these we get that lambda_1 = 2*V_g*lambda_2
results.var.mu.u.Vm(pv) = sol.solveroutput.lambda.ineqnonlin(ngenpv+(1:ngenpv))*...
    2.*results.gen(genpv,VG);
% Determining the generators that are at their limits
[gen_a,gen_b] = determineGenSetsAB(mpc);
idx_bus_sll = gen_a & gen_b;
if sum(idx_bus_sll) > 0
    results.bif.short_name = 'LIB';
    results.bif.full_name = 'limit-induced bifurcation';
    results.bif.gen_sll = find(idx_bus_sll);
else
    results.bif.short_name = 'SNB';
    results.bif.full_name = 'saddle-node bifurcation';
end
results.bif.gen_a = gen_a;
results.bif.gen_b = gen_b;
end