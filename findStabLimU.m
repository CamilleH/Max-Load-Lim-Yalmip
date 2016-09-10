%% In this file, we perform many CPF with different values of the generations u
clear;clc;
% Defining parameters
% global Settings_ch StochMod
Settings_ch = ch_settings_init();
bifSettings = Settings_ch.bifSettings;
define_constants;
toSave = 1;
times = zeros(2,1);
indT = 1;
for conti = 0
    tic
    
    systemName = 'caseKarysNew';
    caseName = 'load4';
    systemName = 'case9static';
    caseName = 'loads568';
    wpNormal = 1;
    
    % System settings
    mpc = openCase(systemName);
    nbus = size(mpc.bus,1);
    caseSettings = getSystemSettings(systemName,caseName);
    StochMod = getSystemStochModel(systemName,caseName,wpNormal);
    rikt_vec_all = caseSettings.rikt_vec_all;
    
    indControlSCOPF = caseSettings.indControl;
    incP = caseSettings.indLoads;
    contingencies = caseSettings.contingencies;
    indStartU = caseSettings.indControl;
    %     startU_all = caseSettings.startU;
    startWP_all = caseSettings.startWP;
    startLoad = caseSettings.startLoad;
    indWP = caseSettings.indWP;
    
    % Load increase direction
    rikt = rikt_vec_all(1,:);
    dirP = zeros(nbus,1);
    dirP(caseSettings.indLoads) = rikt;
    
    % Defining the different cases for the CPF (different u for gen 2 and gen
    % 3)
    nbGenCases = 10;
    casesWP = linspace(0.01,2.25,nbGenCases);%zeros(nbGenCases,nbGenCases);
    casesGen3 = linspace(0,2.7,nbGenCases);%zeros(nbGenCases,nbGenCases);
    [casesWP,casesGen3] = meshgrid(casesWP,casesGen3);    
    
    % Running the CPFs and storing the results
    nbCasesTot = nbGenCases^2;
    indI = zeros(nbCasesTot,1);
    indJ = zeros(nbCasesTot,1);
    startPointUs = zeros(nbCasesTot,1);
    startPointWPs = zeros(nbCasesTot,1);
    failedCases_pf = zeros(nbCasesTot,1);
    Pload_pf = zeros(nbCasesTot,1);
    infoPts_pf = cell(nbCasesTot,1);
    
    % preparing for parfor loop
    for n = 1:nbCasesTot
        i = floor((n-1)/nbGenCases)+1;
        j = mod(n-1,nbGenCases)+1;
        indI(n) = i;
        indJ(n) = j;
        startPointUs(n) = casesGen3(i,1);
        startPointWPs(n) = casesWP(1,j);
    end
    
    for n = 1:nbCasesTot
            c = conti;
            startPointU = startPointUs(n);
            startPointWP = startPointWPs(n);

            startPoint = startLoad;
            chooseStartPoint = 2;
            flagCPF = 0;
            if flagCPF
                try
                    ch_CPF_Dyn;
                catch
                    keyboard;
                end
                try
                    ch_processAfterCPF;
                catch exception
                    keyboard;
                    fprintf('\n The process after the CPF failed.\n');
                    fprintf('The error was:\n%s\n\n',exception.getReport());
                    failedCases_pf(n) = 1;
                end
            end

            % using Yalmip
            Pgen0 = startPointU;
            if isempty(mpc.wind)
                mpc.wind = [0 0 0 0];
            end
            Pwind0 = startPointWP;
            Pload0 = startPoint;
            [sol,Ploadsopt,xopt] = findStabLim_opt(systemName,caseName,dirP,Pgen0,Pwind0,Pload0);
            if sol.problem ~= 0
                failedCases_pf(n) = 1;
            end
            [gen_a,gen_b] = determineGenSetsAB(mpc,xopt,Ploadsopt,systemName,caseName);
            if sum(gen_a&gen_b) > 0
                type = 2;
            else
                type = 1;
            end
            % Storing the values to be plotted
            infoPts_pf{n}.type = type;
            infoPts_pf{n}.gen_a = gen_a;
            infoPts_pf{n}.gen_b = gen_b;
            Pload_pf(n) = Ploadsopt(caseSettings.indLoads);
%             infoPts{n} = infoPt_cur;
            fprintf(1,'CPF %d out of %d completed.\n\n',n,nbCasesTot);
    end
    
    % Splitting in mtrices
    failedCases = zeros(nbGenCases,nbGenCases);
    Pload = zeros(nbGenCases,nbGenCases);
    infoPts = cell(nbGenCases,nbGenCases);
    for n = 1:nbCasesTot
        i = indI(n);
        j = indJ(n);
        failedCases(i,j) = failedCases_pf(n);
        Pload(i,j) = Pload_pf(n);
        infoPts{i,j} = infoPts_pf{n};
    end
    
    % Saving
    if toSave
        systemName = sprintf('%s_WPisPQ_cont%d',systemName,conti);
        saveName = sprintf('Results/KarysNew/Surfaces/%s-%s-PointsSurfaces',datestr(now,'yyyy-mmdd'),systemName);
        save(saveName,'Pload','casesWP','casesGen3','infoPts','failedCases');
    end
    times(indT) = toc;
    indT = indT+1;
end