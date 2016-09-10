function [gen_a,gen_b] = determineGenSetsAB(mpc,x,Ploads)
% This function determines the generators that are PV buses and the
% generators that have reached their Q limit and can be considered as PQ
% buses
define_constants;
idx_buspv = find(mpc.bus(:,BUS_TYPE) == PV);
idx_genpv = ismember(mpc.gen(:,GEN_BUS),idx_buspv);
ngen = size(mpc.gen,1);
[Qgen,Vgen] = computeVarAfterPF(mpc,x,Ploads);
genQlim = abs(Qgen-mpc.gen(idx_genpv,QMAX));
genVlim = abs(Vgen-mpc.gen(idx_genpv,VG));
inGenB = genQlim<1e-5;
inGenA = genVlim<1e-5;
gen_a = zeros(ngen,1);
gen_b = zeros(ngen,1);
gen_a(idx_genpv(inGenA)) = 1;
gen_b(idx_genpv(inGenB)) = 1;