function fmlp_printResults(mpc,x,Ploads,systemName,caseName)
define_constants;
[Qgen,Vgen] = computeVarAfterPF(mpc,x,Ploads,systemName,caseName);

fprintf('Loads:\n');
Ploads'

% Comparing with the limits
fprintf('Reactive power production at the generators and their limits\n');
fprintf('   Bus nb      Qgen       Qlim\n');
[mpc.gen(2:end,1) Qgen mpc.gen(2:end,QMAX)]

% Comparing the voltages at PV buses with their reference
fprintf('Voltages at the PV buses and their references\n')
fprintf('   Bus nb      V       Vref\n');
[mpc.gen(2:end,1) Vgen mpc.gen(2:end,VG)]