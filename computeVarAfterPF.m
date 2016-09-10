function [Qgen,Vgen] = computeVarAfterPF(mpc,x,Ploads)
define_constants;

% Power factor of the loads
n = size(mpc.bus,1);
indNZ = find(mpc.bus(:,PD));
Q_P = zeros(n,1);
Q_P(indNZ) = mpc.bus(indNZ,QD)./mpc.bus(indNZ,PD);

idxpv = find(mpc.bus(:,BUS_TYPE) == PV);
npv = length(idxpv);

Qgen = zeros(npv,1);
Vgen = zeros(npv,1);

[~,Ytis,Mis,~] = buildYMmats(mpc);
for i = 1:npv
    idxi = idxpv(i);
    Qgen(i) = x'*Ytis{idxi}*x+Q_P(idxi)*Ploads(idxi);
    Vgen(i) = sqrt(x'*Mis{idxi}*x);
end

