function [gen_a,gen_b] = determineGenSetsAB(mpc)
% This function determines the generators that are PV buses and the
% generators that have reached their Q limit and can be considered as PQ
% buses
define_constants;
% Scaling to pu
mpc.bus(:,[PD,QD]) = mpc.bus(:,[PD,QD])/mpc.baseMVA;
mpc.gen(:,[PG,QG]) = mpc.gen(:,[PG,QG])/mpc.baseMVA;
mpc.gen(:,[QMAX,QMIN]) = mpc.gen(:,[QMAX,QMIN])/mpc.baseMVA;
if isfield(mpc,'wind') && ~isempty(mpc.wind)
    mpc.wind(:,[PG,QG]) = mpc.wind(:,[PG,QG])/mpc.baseMVA;
    mpc.wind(:,4) = mpc.wind(:,4)/mpc.baseMVA; % wind farm capacity
end
% Computing Qg and Vg
idx_buspv = find(mpc.bus(:,BUS_TYPE) == PV);
idx_genpv = find(ismember(mpc.gen(:,GEN_BUS),idx_buspv));
ngen = size(mpc.gen,1);
[Qgen,Vgen] = computeVarAfterPF(mpc);
genQlim = abs(Qgen-mpc.gen(idx_genpv,QMAX));
genVlim = abs(Vgen-mpc.gen(idx_genpv,VG));
inGenB = genQlim<1e-5;
inGenA = genVlim<1e-5;
gen_a = zeros(ngen,1);
gen_b = zeros(ngen,1);
gen_a(idx_genpv(inGenA)) = 1;
gen_b(idx_genpv(inGenB)) = 1;