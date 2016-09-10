% This script had two modes:
% 1. It searches for the stability limit in different directions in a
% parameter space, uniformly sampled.
% The search uses the function findPtStabLim which solves the non-linear
% optimization problem using the YALMIP interface.
% 2. It plots the results from a previous run with the search mode

clear;clc;close all;

systemName = 'caseKarysNew';
caseName = 'load4';
systemName = 'case9static';
caseName = 'loads568';

% System settings
mpc = openCase(systemName);
nbus = size(mpc.bus,1);
caseSettings = getSystemSettings(systemName,caseName);
define_constants;

% What to do, find or plot?
todo = 'plot';

if strcmp(todo,'find')
    if strcmp(systemName,'case9static')
        nb_pts = 10;
        theta = linspace(0,pi/2,nb_pts);%0:pi/60:pi/2;
        phi = linspace(0,pi/2,nb_pts);%0:pi/60:pi/2;
        len2 = length(theta);
        len1 = length(phi);
        
        infoPts = cell(len1,len2);
        failedCases = zeros(len1,len2);
        load5lim = zeros(len1,len2);
        load6lim = zeros(len1,len2);
        load8lim = zeros(len1,len2);
        iterNb = 1;
        for i = 1:len1
            for j = 1:len2
                rikt = [cos(phi(i))*sin(theta(j));sin(phi(i))*sin(theta(j));cos(theta(j))];
                % Load increase direction
                dirP = zeros(nbus,1);
                dirP(caseSettings.indLoads) = rikt;
                Pgen0 = mpc.gen(caseSettings.indControl,PG);
                if isempty(mpc.wind)
                    mpc.wind = [0 0 0 0];
                end
                Pwind0 = 0;
                Pload0 = mpc.bus(caseSettings.indLoads,PD);
                Pload0(1:end) = 0.001;
                [sol,Ploadsopt,xopt] = findPtStabLim(systemName,caseName,dirP,Pgen0,Pwind0,Pload0);
                if sol.problem ~= 0
                    failedCases(i,j) = 1;
                end
                [gen_a,gen_b] = determineGenSetsAB(mpc,xopt,Ploadsopt);
                if sum(gen_a&gen_b) > 0
                    type = 2;
                else
                    type = 1;
                end
                % Storing the values to be plotted
                infoPts{i,j}.type = type;
                infoPts{i,j}.gen_a = gen_a;
                infoPts{i,j}.gen_b = gen_b;
                load568_opt = Ploadsopt(caseSettings.indLoads);
                load5lim(i,j) = load568_opt(1);
                load6lim(i,j) = load568_opt(2);
                load8lim(i,j) = load568_opt(3);
                %             infoPts{n} = infoPt_cur;
                fprintf(1,'CPF %d out of %d completed.\n\n',iterNb,len1*len2);
                iterNb = iterNb+1;
            end
        end
        saveName = sprintf('Results/KarysNew/Surfaces/ieee9static-from0-essai.mat');
        save(saveName,'infoPts','load5lim','load6lim','load8lim','failedCases');
    end
end

%% Plot
if strcmp(todo,'plot')
    saveName = sprintf('Results/KarysNew/Surfaces/ieee9static-from0-essai.mat');
    resStabLim = load(saveName);
    len1 = size(resStabLim.infoPts,1);
    len2 = size(resStabLim.infoPts,2);
    colorMap = zeros(len1,len2,3);
    allTypes = zeros(len1,len2);
    isSLL = zeros(len1,len2);
    isSNB = zeros(len1,len2);
    
    for i = 1:len1
        for j = 1:len2
            infoPt_cur = resStabLim.infoPts{i,j};
            % We build a binary code for the points. First digit: type, then
            % one digit per gen for gen_a and for gen_b.
            typeBin = sprintf('%d',[infoPt_cur.type-1 infoPt_cur.gen_a.' infoPt_cur.gen_b.']);
            allTypes(i,j) = bin2dec(typeBin);
            isSLL(i,j) = sum(infoPt_cur.gen_a & infoPt_cur.gen_b);
            isSNB(i,j) = ~isSLL(i,j);
        end
    end
    allTypesSNB = allTypes(isSNB == 1);
    allTypesSLL = allTypes(isSLL == 1);
    
    % allTypes(end,end) = 9;
    diffTypes = unique(allTypes);
    diffTypesSNB = unique(allTypesSNB);
    diffTypesSLL = unique(allTypesSLL);
    colorIsSNB = ismember(diffTypes,diffTypesSNB);
    colorIsSLL = ismember(diffTypes,diffTypesSLL);
    nbcolors = length(diffTypes);
    colorTypes = jet(nbcolors);
    namecolors = cell(nbcolors,1);
    genas = zeros(nbcolors,3);
    genbs = zeros(nbcolors,3);
    idxsnb = 1;
    idxsll = 1;
    for i = 1:nbcolors
        [idxi,idxj] = find(allTypes == diffTypes(i));
        genas(i,:) = resStabLim.infoPts{idxi(1),idxj(1)}.gen_a;
        genbs(i,:) = resStabLim.infoPts{idxi(1),idxj(1)}.gen_b;
        if colorIsSNB(i)
            namecolors{i} = sprintf('SNB %d',idxsnb);
            idxsnb = idxsnb+1;
        else
            namecolors{i} = sprintf('SLL %d',idxsll);
            idxsll = idxsll+1;
        end
    end
    % figure
    for i = 1:len1
        for j = 1:len2
            infoPt_cur = resStabLim.infoPts{i,j};
            currColor = colorTypes(diffTypes == allTypes(i,j),:);
            colorMap(i,j,:) = currColor;
        end
    end
    
    figure
    hs = surf(resStabLim.load5lim,resStabLim.load6lim,resStabLim.load8lim,colorMap);
    xlabel('Load 5');
    ylabel('Load 6');
    zlabel('Load 8');
    colormap(colorTypes);
    c = colorbar;
    set(gca,'CLim',[0 1]);
    ytick_width = 1/nbcolors;
    ytick_pos = (ytick_width/2):ytick_width:(1-ytick_width/2);
    set(c,'YTick',ytick_pos,'YTickLabel',namecolors);
    
    % Print table with the gen sets
    txtSurface = cell(nbcolors,3);
    for i = 1:nbcolors
        txtSurface{i,1} = namecolors{i};
        txtSurface{i,2} = mat2str(genas(i,:));
        txtSurface{i,3} = mat2str(genbs(i,:));
    end
    htable = figure('Position',[200,100,250,250]);
    cnames = {'Surface','Set A','Set B'};
    tabCase = uitable('Parent',htable,'Position',[20 20 200 200],...
        'ColumnName',cnames,'Data',txtSurface,...
        'ColumnWidth',{60,60,60});
end