function Jac = buildJacobiansRec(mpc,x,gen_a,gen_b)
define_constants;
n = size(mpc.bus,1);
idxpv = find(mpc.bus(:,BUS_TYPE) == PV);
npv = length(idxpv);
idxpq = find(mpc.bus(:,BUS_TYPE) == PQ);
npq = length(idxpq);
idxGenA = find(gen_a);
idxGenB = find(gen_b);

[Yis,Ytis,Mis,Mislack_Vr,Mislack_Vq] = buildYMmats(mpc);

n = npv+npq+1;
Jac = [];
% Active power flow equations at PV buses
for i = 1:npv
    idxi = idxpv(i);
    Jac = [Jac;...
        x'*(Yis{idxi}+Yis{idxi}')];
end
% Active power flow equations at PQ buses
for i = 1:npq
    idxi = idxpq(i);
    Jac = [Jac;...
        x'*(Yis{idxi}+Yis{idxi}')];
end
% Reactive power flow equations at PQ buses
for i = 1:npq
    idxi = idxpq(i);
    Jac = [Jac;...
        x'*(Ytis{idxi}+Ytis{idxi}')];
end
for i = 1:length(idxGenB)
    idxGeni = idxGenB(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Jac = [Jac;...
        x'*(Ytis{idxi}+Ytis{idxi}')];
end
for i = 1:length(idxGenA);
    idxGeni = idxGenA(i);
    idxi = mpc.gen(idxGeni,GEN_BUS);
    Jac = [Jac;...
        x'*(Mis{idxi}+Mis{idxi}')];
end
Jac = [Jac;...
    x'*(Mislack_Vr+Mislack_Vr')];
Jac = [Jac;...
    x'*(Mislack_Vq+Mislack_Vq')];
% Jac(:,n+1) = []; % Removing Vq for the slack from x
    