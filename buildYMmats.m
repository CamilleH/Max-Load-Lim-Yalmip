function [Yis,Ytis,Mis,Mislack_Vr,Mislack_Vq] = buildYMmats(mpc)
define_constants;

[Y,~] = makeYbus(mpc);
n = size(mpc.bus,1);
idxslack = find(mpc.bus(:,BUS_TYPE) == REF);

alleis = eye(n);
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