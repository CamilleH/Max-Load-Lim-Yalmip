function [Qgen,Vgen] = computeVarAfterPF(mpc)
define_constants;

% x is made of the real and imaginary parts of the voltage
V = mpc.bus(:,VM).*exp(1i*mpc.bus(:,VA)*pi/180);
x = [real(V);imag(V)];

idxpv = find(mpc.bus(:,BUS_TYPE) == PV);
npv = length(idxpv);
Qgen = zeros(npv,1);
Vgen = zeros(npv,1);

[~,Ytis,Mis,~] = buildYMmats(mpc);
for i = 1:npv
    idxi = idxpv(i);
    Qgen(i) = x'*Ytis{idxi}*x+mpc.bus(idxi,QD);
    Vgen(i) = sqrt(x'*Mis{idxi}*x);
end

