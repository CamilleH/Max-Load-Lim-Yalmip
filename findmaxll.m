function results = findmaxll(mpc,idxVarPQ,dir)
% FINDMAXLL finds the maximum loadability limit point in a particular
% direction by using the YALMIP implementation that is formulated with the
% voltages expressed in rectangular coordinates.
%    DIR is the direction of load increase
%% Open the system file and creates the necessary matrices
define_constants;
n = size(mpc.bus,1);
mu = mpc.bus(idxVarPQ,PD);

%% Scale to per-unit
mpc.bus(:,[PD,QD]) = mpc.bus(:,[PD,QD])/mpc.baseMVA;
mpc.gen(:,[PG,QG]) = mpc.gen(:,[PG,QG])/mpc.baseMVA;
mpc.gen(:,[QMAX,QMIN]) = mpc.gen(:,[QMAX,QMIN])/mpc.baseMVA;
if isfield(mpc,'wind') && ~isempty(mpc.wind)
    mpc.wind(:,[PG,QG]) = mpc.wind(:,[PG,QG])/mpc.baseMVA;
    mpc.wind(:,4) = mpc.wind(:,4)/mpc.baseMVA; % wind farm capacity
end

%% Optimization problem
% Parameters
% gen_a = [];
% gen_b = [];
dirP = zeros(n,1);
dirP(idxVarPQ) = dir;

% optimization problem
x = sdpvar(2*n,1);
% v = sdpvar(2*n-1,1);
Ploads = sdpvar(n,1);
% assign(v,ones(2*n-1,1));
assign(x,ones(2*n,1));
assign(Ploads,mpc.bus(:,PD));
lambda = sdpvar(1,1);
Objective = -lambda;
Constraints = fmlp_buildConstraints(mpc,x,Ploads,[],idxVarPQ,[],[]);
Constraints = [Constraints,...
    Ploads == mpc.bus(:,PD)+lambda*dirP];
options = sdpsettings('verbose',0,'showprogress',1,'fmincon.Algorithm', 'interior-point','usex0',1);
sol = optimize(Constraints,Objective,options);

% Extracting the result
Ploadsopt1 = value(Ploads);
lambda = norm(Ploadsopt1(idxVarPQ)-mu);
results.Ploads = value(Ploads);
results.x = value(x);
results.lambda = lambda;
end