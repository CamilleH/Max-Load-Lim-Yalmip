function [sol,Ploadsopt,xopt] = findPtStabLim(systemName,caseName,dirP,Pgen0,Pwind0,Pload0)
% Open the system file and creates the necessary matrices
define_constants;
mpc = openCase(systemName);
caseSettings = getSystemSettings(systemName,caseName);
System = ch_system_init(mpc,caseSettings);
n = System.indices.nbus;

%% Setting the production in conventional gens, wind power, and load
mpc.gen(caseSettings.indControl,PG) = Pgen0;
mpc.wind(:,2) = Pwind0;
mpc.bus(caseSettings.indLoads,PD) = Pload0;
mpc.bus(caseSettings.indLoads,QD) = Pload0.*System.Q_P(caseSettings.indLoads);

%% Optimization problem 1: Look in a particular direction
% Parameters
idxVarPQ = caseSettings.indLoads;
gen_a = [];
gen_b = [];

% optimization problem
x = sdpvar(2*n,1);
v = sdpvar(2*n,1);
Ploads = sdpvar(n,1);
assign(v,ones(2*n,1));
assign(x,ones(2*n,1));
assign(Ploads,2*mpc.bus(:,PD));
eta = sdpvar(1,1);
Objective = -eta;
Constraints = fmlp_buildConstraints(mpc,x,Ploads,v,idxVarPQ,gen_a,gen_b);
Constraints = [Constraints,...
    Ploads == mpc.bus(:,PD)+eta*dirP];
options = sdpsettings('verbose',1,'showprogress',1,'fmincon.Algorithm', 'interior-point','usex0',1);
sol = optimize(Constraints,Objective,options);

%% Getting the results
Ploadsopt = value(Ploads);
xopt = value(x);