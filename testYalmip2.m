% Test file for using Yalmip for solving the nonrelaxed problemsclear, clc;
clear, clc;
systemName = 'caseKarysNew';
caseName = 'load4';
systemName = 'case9static';
caseName = 'load5';

% ch_loadFigure();

%% Open the system file and creates the necessary matrices
define_constants;
mpc = loadcase(systemName);
n = size(mpc.bus,1);
idxpv = find(mpc.bus(:,BUS_TYPE) == PV);
npv = length(idxpv);
idxpq = find(mpc.bus(:,BUS_TYPE) == PQ);
npq = length(idxpq);
mu = mpc.bus([5;7;9],PD)/mpc.baseMVA;
figure(1);
hold on;
plot3(mu(1),mu(2),mu(3),'*','Color',[1 0 0],'MarkerSize',10,'MarkerFaceColor',[1 0 0]);

%% Optimization problem 1: Look in a particular direction
% Parameters
dirP = zeros(n,1);
% idxVarPQ = [5;6;8];
% dirP(idxVarPQ) = [4.088;0.9007;0.2279]-mpc.bus(idxVarPQ,PD);
idxVarPQ = 5;
dirP(idxVarPQ) = 1;
dirP = dirP/norm(dirP);
gen_a = [];
gen_b = [];

% optimization problem
idxVarPQ = [5;7;9];
x = sdpvar(2*n,1);
v = sdpvar(2*n-1,1);
Ploads = sdpvar(n,1);
assign(v,ones(2*n-1,1));
assign(x,ones(2*n,1));
assign(Ploads,mpc.bus(:,PD));
eta = sdpvar(1,1);
Objective = -eta;
Constraints = fmlp_buildConstraints(mpc,x,Ploads,v,idxVarPQ,gen_a,gen_b);
Constraints = [Constraints,...
    Ploads == mpc.bus(:,PD)+eta*dirP];
options = sdpsettings('verbose',2,'showprogress',1,'fmincon.Algorithm', 'interior-point','usex0',1);
sol = optimize(Constraints,Objective,options);

%% Displaying the results
Ploadsopt1 = value(Ploads);
xopt1 = value(x);
fmlp_printResults(mpc,xopt1,Ploadsopt1)
fprintf('Distance to mu: %.4g\n',norm(Ploadsopt1([5;6;8])-mu));
[gen_a,gen_b] = determineGenSetsAB(mpc,xopt1,Ploadsopt1);
plotLoadsOnFig(1,Ploadsopt1,[5 7 9])

Jac=buildJacobiansRec(systemName,caseName,mpc,xopt3,gen_a,gen_b);
[Veig,Deig] = eig(Jac);

keyboard;

%% Optimization problem 2: look for the closest point on the identified 
% surface

% Parameters
idxVarPQ = [5;7;9];

% optimization problem
x = sdpvar(2*n,1);
v = sdpvar(2*n-2,1);
Ploads = sdpvar(n,1);
assign(v,ones(2*n-2,1));
assign(x,xopt1);
assign(Ploads,2*mpc.bus(:,PD));
Objective = (Ploads(idxVarPQ)-mu)'*(Ploads(idxVarPQ)-mu);
Constraints = fmlp_buildConstraints(mpc,x,Ploads,v,idxVarPQ,gen_a,gen_b);
options = sdpsettings('verbose',2,'showprogress',1,'fmincon.Algorithm', 'interior-point','usex0',1);
sol = optimize(Constraints,Objective,options);

%% Displaying the results
Ploadsopt2 = value(Ploads);
xopt2 = value(x);
vopt2 = value(v);
fmlp_printResults(mpc,xopt2,Ploadsopt2)
fprintf('Distance to mu: %.4g\n',norm(Ploadsopt2(idxVarPQ)-mu));
[gen_a,gen_b] = determineGenSetsAB(mpc,xopt2,Ploadsopt2);
plotLoadsOnFig(1,Ploadsopt2,[5 7 9])

keyboard;

%% Opti 3
x = sdpvar(2*n,1);
v = sdpvar(2*n-2,1);
dX_pre = sdpvar(2*n,1);
dX_post = sdpvar(2*n,1);
critnode = 9;
Ploads = sdpvar(n,1);
assign(v,ones(2*n-2,1));
assign(x,xopt2);
assign(Ploads,Ploadsopt2);
assign(dX_pre,ones(2*n,1));
assign(dX_post,ones(2*n,1));
Objective = (Ploads(idxVarPQ)-mu)'*(Ploads(idxVarPQ)-mu);
constrArgs = {dX_pre,dX_post,critnode};
Constraints = fmlp_buildConstraints(mpc,x,Ploads,constrArgs,idxVarPQ,gen_a,gen_b);
options = sdpsettings('verbose',2,'showprogress',1,'fmincon.Algorithm','interior-point',...
    'fmincon.TolCon',2e-5,'usex0',1);
sol = optimize(Constraints,Objective,options);
%% Displaying the results
Ploadsopt3 = value(Ploads);
xopt3 = value(x);
fmlp_printResults(mpc,xopt3,Ploadsopt3)
fprintf('Distance to mu: %.4g\n',norm(Ploadsopt3(idxVarPQ)-mu));
[gen_a,gen_b] = determineGenSetsAB(mpc,xopt3,Ploadsopt3);
plotLoadsOnFig(1,Ploadsopt3,[5 7 9])

keyboard;

%% Opti 4
gen_b(3) = 0;
Jac=buildJacobiansRec(systemName,caseName,mpc,xopt3,gen_a,gen_b);
[Veig,Deig] = eig(Jac);

x = sdpvar(2*n,1);
v = sdpvar(2*n-2,1);
Ploads = sdpvar(n,1);
assign(v,ones(2*n-2,1));
assign(x,xopt3);
assign(Ploads,Ploadsopt3);
Objective = (Ploads(idxVarPQ)-mu)'*(Ploads(idxVarPQ)-mu);
Constraints = fmlp_buildConstraints(mpc,x,Ploads,v,idxVarPQ,gen_a,gen_b);
options = sdpsettings('verbose',2,'showprogress',1,'fmincon.Algorithm', 'interior-point','usex0',1);
sol = optimize(Constraints,Objective,options);
%% Displaying the results
Ploadsopt4 = value(Ploads);
xopt4 = value(x);
fmlp_printResults(mpc,xopt4,Ploadsopt4)
fprintf('Distance to mu: %.4g\n',norm(Ploadsopt4(idxVarPQ)-mu));
[gen_a,gen_b] = determineGenSetsAB(mpc,xopt4,Ploadsopt3);
plotLoadsOnFig(1,Ploadsopt4,[5 6 8])