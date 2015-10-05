function Constraints = fmlp_buildConstraints(mpc,x,Ploads,addVar,systemName,idxVarPQ,gen_a,gen_b)
%% Getting the settings for this case
define_constants;
% caseSettings = getSystemSettings(systemName,caseName);
% System = ch_system_init(mpc,caseSettings);
n = size(mpc.bus,1);

% Power factor of the loads
indNZ = find(mpc.bus(:,PD));
Q_P = zeros(n,1);
Q_P(indNZ) = mpc.bus(indNZ,QD)./mpc.bus(indNZ,PD);

idxslack = find(mpc.bus(:,BUS_TYPE) == REF);
% if isempty(gen_a) && isempty(gen_b)
    idxpq = find(mpc.bus(:,BUS_TYPE) == PQ);
    idxpv = find(mpc.bus(:,BUS_TYPE) == PV);
% else
%     idxpv = mpc.gen(gen_a == 1,GEN_BUS);
%     idxpq = setdiff(1:n,[idxpv;idxslack]);
% end

npv = length(idxpv);
npq = length(idxpq);
idxGenA = find(gen_a);
idxGenB = find(gen_b);
idxFixedPQ = setdiff(1:n,idxVarPQ);

if isempty(idxGenB) && isempty(idxGenA)
    mode = 'pvpq';
elseif sum(gen_a & gen_b) > 0
    mode = 'sll';
    dX_pre = addVar{1};
    dX_post = addVar{2};
    critnode = addVar{3};
else
    mode = 'snb';
    v = addVar;
end

% if strcmp(mode,'sll')
%     idxSLL = find(gen_a&gen_b);
% end

fprintf('The current surface is %s \n',mode);

% Build matrices
[Yis,Ytis,Mis,Mislack_Vr,Mislack_Vq] = buildYMmats(mpc);

if ~isfield(mpc,'wind') || isempty(mpc.wind)
    mpc.wind = [0 0 0 0];
end

%% Building the constraints
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
        x'*Ytis{idxi}*x+Q_P(idxi)*Ploads(idxi) == 0];
end

if strcmp(mode,'snb') || strcmp(mode,'sll')
    % If we are looking specifically for SNB or SLL with fixed gen_a and gen_b
    % sets
    for i = 1:length(idxGenB)
        idxGeni = idxGenB(i);
        idxi = mpc.gen(idxGeni,GEN_BUS);
        Constraints = [Constraints,...
            x'*Ytis{idxi}*x+Q_P(idxi)*Ploads(idxi) - mpc.gen(idxi,QMAX) == 0];
    end
end

% reactive power limits on gens
if strcmp(mode,'snb') || strcmp(mode,'sll')
    listGen = mpc.gen(idxGenA,GEN_BUS);
elseif strcmp(mode,'pvpq')
    listGen = idxpv;
else
    error('Wrong mode');
end
for i = 1:length(listGen)
    idxi = listGen(i);
    Constraints = [Constraints,...
        x'*Ytis{idxi}*x+Q_P(idxi).*Ploads(idxi) <= mpc.gen(mpc.gen(:,GEN_BUS) == idxi,QMAX)];
end

% V = Vref at slack and Vq = 0
Constraints = [Constraints,...
    x'*Mislack_Vr*x == mpc.gen(mpc.gen(:,GEN_BUS) == idxslack,VG)^2];
Constraints = [Constraints,...
    x'*Mislack_Vq*x == 0];

% V limits at bus where gens are connected: upper limits for gen_b,
% enforced for gen a
if strcmp(mode,'snb') || strcmp(mode,'sll')
    for i = 1:length(idxGenA)
        idxGeni = idxGenA(i);
        idxi = mpc.gen(idxGeni,GEN_BUS);
        Constraints = [Constraints,...
            x'*Mis{idxi}*x == mpc.gen(mpc.gen(:,GEN_BUS) == idxi,VG)^2];
    end
end
if strcmp(mode,'snb') || strcmp(mode,'sll')
    listGen = mpc.gen(idxGenB,GEN_BUS);
elseif strcmp(mode,'pvpq')
    listGen = idxpv;
else
    error('Wrong mode');
end
for i = 1:length(listGen)
    idxi = listGen(i);
    Constraints = [Constraints,...
        x'*Mis{idxi}*x <= mpc.gen(mpc.gen(:,GEN_BUS) == idxi,VG)^2];
end

% Additional constraints for SNB = taking derivative of all
% equality constraints above to form the Jacobians
if strcmp(mode,'snb')
    Constraints = [Constraints,...
        v'*v == 1];
    Jac = buildJacobiansRec(mpc,x,gen_a,gen_b);
%     Jac(n+[n-1 n],:) = []; % Removing the equations for the slack
%     Jac(:,[1 n+1]) = []; % Removing the variables Vr and Vq of the slack
    Jac([2*n-1 2*n],:) = []; % Removing the equations for the slack
    Jac(:,[1 n+1]) = []; % Removing the variables Vr and Vq of the slack
    Constraints = [Constraints,...
       Jac*v == 0];
end

% Additional constraints for SLL = the variations of lambda (continuation
% parameter at the critical bus) with increased critical bus voltage is
% increasing on post-limit PV curve and decreasing on pre-limit PV curve
if strcmp(mode,'sll')
    ngensll = find(gen_a & gen_b);
    nbussll = mpc.gen(ngensll,GEN_BUS);
    gen_a_pre = gen_a;
    gen_b_pre = gen_b;
    gen_b_pre(nbussll) = 0;
    Jac_pre = buildJacobiansRec(mpc,x,gen_a_pre,gen_b_pre);
    gen_a_post = gen_a;
    gen_b_post = gen_b;
    gen_a_post(nbussll) = 0;
    Jac_post = buildJacobiansRec(mpc,x,gen_a_post,gen_b_post);
    
    % Constraint SLL 1: direction of increase of both Vr (critical) and
    % Vq (critical)
    Vr_inc = zeros(2*n,2*n);
    Vr_inc(critnode,critnode) = 1;
    Vq_inc = zeros(2*n,2*n);
    Vq_inc(critnode,critnode) = 1;
    Constraints = [Constraints,...
        x'*Vr_inc*dX_pre == 1,...
        x'*Vq_inc*dX_pre == 1,...
        x'*Vr_inc*dX_post == 1,...
        x'*Vq_inc*dX_post == 1];
    
    % Constraint SLL 2: no change in voltage for PV buses
    idxpv_pre = find(gen_a_pre);
    for i = 1:length(idxpv_pre)
        idxgenpvi = idxpv_pre(i);
        idxbuspvi = mpc.gen(idxgenpvi,GEN_BUS);
        Vr_pvi_pre = zeros(2*n,2*n);
        Vq_pvi_pre = zeros(2*n,2*n);
        Vr_pvi_pre(idxbuspvi,idxbuspvi) = 1;
        Vq_pvi_pre(idxbuspvi,idxbuspvi) = 1;
        Constraints = [Constraints,...
            x'*Vr_pvi_pre*dX_pre == 0,...
            x'*Vq_pvi_pre*dX_pre == 0];
    end
    idxpv_post = find(gen_a_post);
    for i = 1:length(idxpv_post)
        idxgenpvi = idxpv_post(i);
        idxbuspvi = mpc.gen(idxgenpvi,GEN_BUS);
        Vr_pvi_post = zeros(2*n,2*n);
        Vq_pvi_post = zeros(2*n,2*n);
        Vr_pvi_post(idxbuspvi,idxbuspvi) = 1;
        Vq_pvi_post(idxbuspvi,idxbuspvi) = 1;
        Constraints = [Constraints,...
            x'*Vr_pvi_post*dX_post == 0,...
            x'*Vq_pvi_post*dX_post == 0];
    end
    
    % Constraint SLL 3: V increases => PL-PG decreases prior to limit and
    % increases post limit in the critical node. Also, make sure that the
    % power factor is constant.
    Constraints = [Constraints,...
        Jac_pre(critnode,:)*dX_pre >= 0];
    Constraints = [Constraints,...
        (Jac_pre(critnode,:)-Q_P(critnode)*Jac_pre(n+critnode,:))*dX_pre == 0];
    Constraints = [Constraints,...
        Jac_post(critnode,:)*dX_post <= 0];
    Constraints = [Constraints,...
        (Jac_post(critnode,:)-Q_P(critnode)*Jac_post(n+critnode,:))*dX_pre == 0];
end

% Ploads positive
% Constraints = [Constraints,...
%     Ploads >= 0];
% Constraints = [Constraints,...
%     Ploads(idxVarPQ) >= 0];