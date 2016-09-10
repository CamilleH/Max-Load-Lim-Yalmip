% Test file for using Yalmip for solving the nonrelaxed problemsclear, clc;
clear, clc;
systemName = 'case9static';

%% Open the system file and creates the necessary matrices
define_constants;
mpc = openCase(systemName);
caseSettings = getSystemSettings(systemName,'load5');
System = ch_system_init(mpc,caseSettings);
n = System.indices.nbus;
idxpv = find(mpc.bus(:,BUS_TYPE) == PV);
npv = length(idxpv);
idxpq = find(mpc.bus(:,BUS_TYPE) == PQ);
npq = length(idxpq);
idxslack = find(mpc.bus(:,BUS_TYPE) == REF);
nslack = length(idxslack);
dirP = zeros(n,1);
dirP(caseSettings.indLoads) = caseSettings.rikt_vec_all;

% Add very small resistance to all lines
Y = System.Ybus_stat;
% for i = 1:n
%     Y(
% end
% System.Ybus_stat(System.Ybus_stat ~= 0) = System.Ybus_stat(System.Ybus_stat ~= 0)+1e-4;
% System.Ybus_stat(eye(n)==1) = System.Ybus_stat(eye(n)==1) - 2e-4;

alleis = eye(n);

% Create the Yi and Mi matrices for each node
Yis = cell(n,1);
Ytis = cell(n,1);
for i = 1:n
    ei = alleis(:,i);
    Yi = ei*ei.'*Y;
    Yis{i} = 0.5*[real(Yi+Yi.') imag(Yi.'-Yi);
        imag(Yi-Yi.') real(Yi+Yi.')];
    Ytis{i} = -0.5*[imag(Yi+Yi.') real(Yi-Yi.');
        real(Yi.'-Yi) imag(Yi+Yi.')];
end

Mis = cell(n,1);
for i = 1:n
    ei = alleis(:,i);
    Mis{i} = [ei*ei' zeros(n);
        zeros(n) ei*ei'];
end
es = alleis(:,idxslack);
Mislack_Vr = [es*es' zeros(n);
    zeros(n,2*n)];
Mislack_Vq = [zeros(n,2*n);
    zeros(n) es*es'];

if isempty(mpc.wind)
    mpc.wind = [0 0 0 0];
end

%% Yalmip code

x = sdpvar(2*n,1);
assign(x,ones(2*n,1));
eta = sdpvar(1,1);
Objective = -eta;
Constraints = [];
% active power flows at pv and pq
for i = 1:npv
    idxi = idxpv(i);
    Constraints = [Constraints,...
        x'*Yis{idxi}*x+mpc.bus(idxi,PD)+eta*dirP(idxi) == sum(mpc.gen(mpc.gen(:,GEN_BUS) == idxi,PG))+sum(mpc.wind(mpc.wind(:,1) == idxi,2))];
end
for i = 1:npq
    idxi = idxpq(i);
    Constraints = [Constraints,...
        x'*Yis{idxi}*x+mpc.bus(idxi,PD)+eta*dirP(idxi) == sum(mpc.gen(mpc.gen(:,GEN_BUS) == idxi,PG))+sum(mpc.wind(mpc.wind(:,1) == idxi,2))];
end
% reactive power flows at pq and gen_b in pv
for i = 1:npq
    idxi = idxpq(i);
    Constraints = [Constraints,...
        x'*Ytis{idxi}*x+mpc.bus(idxi,QD)+eta*System.Q_P(idxi)*dirP(idxi) == 0];
end
% reactive power limits on gens
for i = 1:npv
    idxi = idxpv(i);
    Constraints = [Constraints,...
        x'*Ytis{idxi}*x+mpc.bus(idxi,QD)+eta*System.Q_P(idxi).*dirP(idxi) <= mpc.gen(mpc.gen(:,GEN_BUS) == idxi,QMAX)];
end
% V = Vref at slack and Vq = 0
Constraints = [Constraints,...
    x'*Mislack_Vr*x == mpc.gen(mpc.gen(:,GEN_BUS) == idxslack,VG)^2];
Constraints = [Constraints,...
    x'*Mislack_Vq*x == 0];
% V limits at bus where gens are connected
for i = 1:npv
    idxi = idxpv(i);
    Constraints = [Constraints,...
        x'*Mis{idxi}*x <= mpc.gen(mpc.gen(:,GEN_BUS) == idxi,VG)^2];
end
options = sdpsettings('verbose',2,'showprogress',1,'fmincon.Algorithm', 'interior-point','usex0',1);
sol = optimize(Constraints,Objective,options);

%% Display results
xopt = value(x);
etaopt = value(eta);
Vropt = xopt(1:n);
Vqopt = xopt(n+1:end);
Vopt = Vropt+1i*Vqopt;
Qgen = zeros(npv,1);
Vgen = zeros(npv,1);

for i = 1:npv
    idxi = idxpv(i);
    Qgen(i) = xopt'*Ytis{idxi}*xopt+mpc.bus(idxi,QD)+etaopt*System.Q_P(idxi)*dirP(idxi);
    Vgen(i) = sqrt(xopt'*Mis{idxi}*xopt);
end

% Comparing with the limits
fprintf('Reactive power production at the generators and their limits\n');
fprintf('   Bus nb      Qgen       Qlim\n');
[mpc.gen(2:end,1) Qgen mpc.gen(2:end,QMAX)]

% Comparing the voltages at PV buses with their reference
fprintf('Voltages at the PV buses and their references')
fprintf('   Bus nb      V       Vref\n');
[mpc.gen(2:end,1) Vgen mpc.gen(2:end,VG)]

%% Now same thing with the constraints for SNB
gen_a = [0;0;1];
gen_b = [0;1;0];
idxGenA = find(gen_a);
idxGenB = find(gen_b);

% Jacobian from previous solution
JacPrev = buildJacobiansRec(mpc,npv,npq,idxpv,idxpq,idxGenA,idxGenB,Yis,Ytis,Mis,Mislack_Vq,Mislack_Vr,xopt);
[Veig,Deig] = eig(JacPrev);

% optimization problem
x = sdpvar(2*n,1);
v = sdpvar(2*n,1);
assign(x,xopt);
assign(v,ones(2*n,1));
eta = sdpvar(1,1);
Objective = -eta;
Constraints = [];
% active power flows at pv and pq
for i = 1:npv
    idxi = idxpv(i);
    Constraints = [Constraints,...
        x'*Yis{idxi}*x+mpc.bus(idxi,PD)+eta*dirP(idxi) == sum(mpc.gen(mpc.gen(:,GEN_BUS) == idxi,PG))+sum(mpc.wind(mpc.wind(:,1) == idxi,2))];
end
for i = 1:npq
    idxi = idxpq(i);
    Constraints = [Constraints,...
        x'*Yis{idxi}*x+mpc.bus(idxi,PD)+eta*dirP(idxi) == sum(mpc.gen(mpc.gen(:,GEN_BUS) == idxi,PG))+sum(mpc.wind(mpc.wind(:,1) == idxi,2))];
end
% reactive power flows at pq and gen_b in pv
for i = 1:npq
    idxi = idxpq(i);
    Constraints = [Constraints,...
        x'*Ytis{idxi}*x+mpc.bus(idxi,QD)+eta*System.Q_P(idxi)*dirP(idxi) == 0];
end
for i = 1:length(idxGenB)
    idxGeni = idxGenB(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Constraints = [Constraints,...
        x'*Ytis{idxi}*x+mpc.bus(idxi,QD)+eta*System.Q_P(idxi)*dirP(idxi) - mpc.gen(idxi,QMAX) == 0];
end
% reactive power limits on gens
for i = 1:length(idxGenA)
    idxGeni = idxGenA(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Constraints = [Constraints,...
        x'*Ytis{idxi}*x+mpc.bus(idxi,QD)+eta*System.Q_P(idxi).*dirP(idxi) <= mpc.gen(mpc.gen(:,GEN_BUS) == idxi,QMAX)];
end
% V = Vref at slack and Vq = 0
Constraints = [Constraints,...
    x'*Mislack_Vr*x == mpc.gen(mpc.gen(:,GEN_BUS) == idxslack,VG)^2];
Constraints = [Constraints,...
    x'*Mislack_Vq*x == 0];
% V limits at bus where gens are connected: upper limits for gen_b,
% enforced for gen a
for i = 1:length(idxGenA)
    idxGeni = idxGenA(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Constraints = [Constraints,...
        x'*Mis{idxi}*x == mpc.gen(mpc.gen(:,GEN_BUS) == idxi,VG)^2];
end
for i = 1:length(idxGenB)
    idxGeni = idxGenB(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Constraints = [Constraints,...
        x'*Mis{idxi}*x <= mpc.gen(mpc.gen(:,GEN_BUS) == idxi,VG)^2];
end
% Additional constraints for SNB = taking derivative of all
% equality constraints above to form the Jacobians
for i = 1:npv
    idxi = idxpv(i);
    Constraints = [Constraints,...
        x'*(Yis{idxi}+Yis{idxi}')*v == 0];
end
for i = 1:npq
    idxi = idxpq(i);
    Constraints = [Constraints,...
        x'*(Yis{idxi}+Yis{idxi}')*v== 0];
end
for i = 1:npq
    idxi = idxpq(i);
    Constraints = [Constraints,...
        x'*(Ytis{idxi}+Ytis{idxi}')*v == 0];
end
for i = 1:length(idxGenB)
    idxGeni = idxGenB(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Constraints = [Constraints,...
        x'*(Ytis{idxi}+Ytis{idxi}')*v == 0];
end
for i = 1:length(idxGenA);
    idxGeni = idxGenA(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Constraints = [Constraints,...
        x'*(Mis{idxi}+Mis{idxi}')*v == 0];
end
Constraints = [Constraints,...
        x'*(Mislack_Vq+Mislack_Vq')*v == 0];
Constraints = [Constraints,...
        x'*(Mislack_Vr+Mislack_Vr')*v == 0];
Constraints = [Constraints,...
        v'*v == 1];
%options = sdpsettings('verbose',2,'showprogress',1,'fmincon.Algorithm', 'interior-point','fmincon.TolCon',1e-5,'fmincon.TolX',1e-5,'usex0',1);
options = sdpsettings('verbose',2,'showprogress',1,'fmincon.Algorithm', 'interior-point','usex0',1);
sol = optimize(Constraints,Objective,options);

%% Display results
xopt = value(x);
etaopt = value(eta);
Vropt = xopt(1:n);
Vqopt = xopt(n+1:end);
Vopt = Vropt+1i*Vqopt;
Qgen = zeros(npv,1);
Vgen = zeros(npv,1);

for i = 1:npv
    idxi = idxpv(i);
    Qgen(i) = xopt'*Ytis{idxi}*xopt+mpc.bus(idxi,QD)+etaopt*System.Q_P(idxi)*dirP(idxi);
    Vgen(i) = sqrt(xopt'*Mis{idxi}*xopt);
end
fprintf('Eta:\n');
etaopt
fprintf('Loads:\n');
(mpc.bus(:,PD)+etaopt*dirP)'

% Comparing with the limits
fprintf('Reactive power production at the generators and their limits\n');
fprintf('   Bus nb      Qgen       Qlim\n');
[mpc.gen(2:end,1) Qgen mpc.gen(2:end,QMAX)]

% Comparing the voltages at PV buses with their reference
fprintf('Voltages at the PV buses and their references')
fprintf('   Bus nb      V       Vref\n');
[mpc.gen(2:end,1) Vgen mpc.gen(2:end,VG)]

%% Test with another objective function
gen_a = [0;0;1];
gen_b = [0;1;0];
idxGenA = find(gen_a);
idxGenB = find(gen_b);
idxVarPQ = [5;6;8];
idxFixedPQ = setdiff(1:n,idxVarPQ);
mu = mpc.bus(idxVarPQ,PD);

% Jacobians from the previous optimization problem
JacPrev = buildJacobiansRec(mpc,npv,npq,idxpv,idxpq,idxGenA,idxGenB,Yis,Ytis,Mis,Mislack_Vq,Mislack_Vr,xopt);
[Veig,Deig] = eig(JacPrev);

% optimization problem
x = sdpvar(2*n,1);
v = sdpvar(2*n,1);
assign(x,xopt);
assign(v,ones(2*n,1));
Ploads = sdpvar(n,1);
assign(x,ones(2*n,1));
assign(Ploads,2*mpc.bus(:,PD));
eta = sdpvar(1,1);
Objective = (Ploads(idxVarPQ) - mu)'*(Ploads(idxVarPQ) - mu);
Constraints = [];
% Loads constant equal to base case for some loads
for i = 1:length(idxFixedPQ)
    idxi = idxFixedPQ(i);
    Constraints = [Constraints,...
        Ploads(idxi) == mpc.bus(idxi,PD)];
end
% active power flows at pv and pq
for i = 1:npv
    idxi = idxpv(i);
    Constraints = [Constraints,...
        x'*Yis{idxi}*x+Ploads(idxi) == sum(mpc.gen(mpc.gen(:,GEN_BUS) == idxi,PG))+sum(mpc.wind(mpc.wind(:,1) == idxi,2))];
end
for i = 1:npq
    idxi = idxpq(i);
    Constraints = [Constraints,...
        x'*Yis{idxi}*x+Ploads(idxi) == sum(mpc.gen(mpc.gen(:,GEN_BUS) == idxi,PG))+sum(mpc.wind(mpc.wind(:,1) == idxi,2))];
end
% reactive power flows at pq and gen_b in pv
for i = 1:npq
    idxi = idxpq(i);
    Constraints = [Constraints,...
        x'*Ytis{idxi}*x+System.Q_P(idxi)*Ploads(idxi) == 0];
end
for i = 1:length(idxGenB)
    idxGeni = idxGenB(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Constraints = [Constraints,...
        x'*Ytis{idxi}*x+System.Q_P(idxi)*Ploads(idxi) - mpc.gen(idxi,QMAX) == 0];
end
% reactive power limits on gens
for i = 1:length(idxGenA)
    idxGeni = idxGenA(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Constraints = [Constraints,...
        x'*Ytis{idxi}*x+System.Q_P(idxi).*Ploads(idxi) <= mpc.gen(mpc.gen(:,GEN_BUS) == idxi,QMAX)];
end
% V = Vref at slack and Vq = 0
Constraints = [Constraints,...
    x'*Mislack_Vr*x == mpc.gen(mpc.gen(:,GEN_BUS) == idxslack,VG)^2];
Constraints = [Constraints,...
    x'*Mislack_Vq*x == 0];
% V limits at bus where gens are connected: upper limits for gen_b,
% enforced for gen a
for i = 1:length(idxGenA)
    idxGeni = idxGenA(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Constraints = [Constraints,...
        x'*Mis{idxi}*x == mpc.gen(mpc.gen(:,GEN_BUS) == idxi,VG)^2];
end
for i = 1:length(idxGenB)
    idxGeni = idxGenB(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Constraints = [Constraints,...
        x'*Mis{idxi}*x <= mpc.gen(mpc.gen(:,GEN_BUS) == idxi,VG)^2];
end
% Additional constraints for SNB = taking derivative of all
% equality constraints above to form the Jacobians
for i = 1:npv
    idxi = idxpv(i);
    Constraints = [Constraints,...
        x'*(Yis{idxi}+Yis{idxi}')*v == 0];
end
for i = 1:npq
    idxi = idxpq(i);
    Constraints = [Constraints,...
        x'*(Yis{idxi}+Yis{idxi}')*v== 0];
end
for i = 1:npq
    idxi = idxpq(i);
    Constraints = [Constraints,...
        x'*(Ytis{idxi}+Ytis{idxi}')*v == 0];
end
for i = 1:length(idxGenB)
    idxGeni = idxGenB(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Constraints = [Constraints,...
        x'*(Ytis{idxi}+Ytis{idxi}')*v == 0];
end
for i = 1:length(idxGenA);
    idxGeni = idxGenA(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Constraints = [Constraints,...
        x'*(Mis{idxi}+Mis{idxi}')*v == 0];
end
Constraints = [Constraints,...
        x'*(Mislack_Vq+Mislack_Vq')*v == 0];
Constraints = [Constraints,...
        x'*(Mislack_Vr+Mislack_Vr')*v == 0];
% Ploads positive
Constraints = [Constraints,...
        Ploads >= 0];
options = sdpsettings('verbose',2,'showprogress',1,'fmincon.Algorithm', 'interior-point','usex0',1);
sol = optimize(Constraints,Objective,options);

%% Displaying the results
Ploadsopt = value(Ploads);
xopt = value(x);
Vropt = xopt(1:n);
Vqopt = xopt(n+1:end);
Vopt = Vropt+1i*Vqopt;
Qgen = zeros(npv,1);
Vgen = zeros(npv,1);

for i = 1:npv
    idxi = idxpv(i);
    Qgen(i) = xopt'*Ytis{idxi}*xopt+System.Q_P(idxi)*Ploads(idxi);
    Vgen(i) = sqrt(xopt'*Mis{idxi}*xopt);
end

fprintf('Loads:\n');
Ploadsopt'

% Comparing with the limits
fprintf('Reactive power production at the generators and their limits\n');
fprintf('   Bus nb      Qgen       Qlim\n');
[mpc.gen(2:end,1) Qgen mpc.gen(2:end,QMAX)]

% Comparing the voltages at PV buses with their reference
fprintf('Voltages at the PV buses and their references\n')
fprintf('   Bus nb      V       Vref\n');
[mpc.gen(2:end,1) Vgen mpc.gen(2:end,VG)]