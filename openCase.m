function mpc = openCase(systemName,varargin)
define_constants;
Settings_ch = ch_settings_init();
try
    mpc = loadcase(systemName);
catch
    keyboard;
end

% Conversion to per-unit
mpc.bus(:,[PD,QD]) = mpc.bus(:,[PD,QD])/mpc.baseMVA;
mpc.gen(:,[PG,QG]) = mpc.gen(:,[PG,QG])/mpc.baseMVA;
if Settings_ch.bifSettings.includeSLL
    mpc.gen(:,[QMAX,QMIN]) = mpc.gen(:,[QMAX,QMIN])/mpc.baseMVA;
else
    mpc.gen(:,QMAX) = 9999;
    mpc.gen(:,QMIN) = -9999;
    if isfield(mpc,'ex')
        mpc.ex(:,9) = 9999;
    end
end

if isfield(mpc,'wind') && ~isempty(mpc.wind)
    mpc.wind(:,[PG,QG]) = mpc.wind(:,[PG,QG])/mpc.baseMVA;
    mpc.wind(:,4) = mpc.wind(:,4)/mpc.baseMVA; % wind farm capacity
end
% conversion to radian
mpc.bus(:,VA) = mpc.bus(:,VA)*pi/180;

if ~isempty(varargin)
    indNZ = find(mpc.bus(:,PD));
    Q_P = zeros(size(mpc.bus,1),1);
    Q_P(indNZ) = mpc.bus(indNZ,QD)./mpc.bus(indNZ,PD);
    mpc.bus(:,PD) = mpc.bus(:,PD)*1.2;
    mpc.bus(:,QD) = mpc.bus(:,PD).*Q_P;
end

