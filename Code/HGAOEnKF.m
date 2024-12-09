clear
addpath(genpath('../Data/Utils'))
% clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
earth_R = 6367;
period = [1880:2000];
simul_year = length(period);
ensemble_size = 500;
climo_size = 1000;
Mask_year = 1305;

Exp_sets = {'REAL','OSSE'};
iExp_set = 2;
Exp_set = Exp_sets{iExp_set};

Model_types = {'HGAOEnKF'};
iModel_type = 1;
Model_type = Model_types{iModel_type};

QCs = {'WiQC','NoQC'};
iQC = 2;
QC = QCs{iQC};

Prior_raw_names = {'CESM2'};
iPrior_raw_name = 1;
Prior_raw_name = Prior_raw_names{iPrior_raw_name};

Proxy_raw_names = {'PAGES2k' ,'Pseudoproxy'};
iProxy_raw_name = iExp_set;
Proxy_raw_name = Proxy_raw_names{iProxy_raw_name};

Prior_dir = ['../Data/prior/Prior_preprocess/' Prior_raw_name '/'];
Proxy_dir = ['../Data/proxy/Proxy_preprocess/' Proxy_raw_name '/' num2str(Mask_year) '/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the AOEnKF and AOEnKF_B
Locs = [1 2 3 4];
iLoc = 4; % optimal loc for AOEnKF,for OSSE 4, for REAL 3
Loc = Locs(iLoc);
Assim_dir = ['../Data/Assim_output/AOEnKF/' Exp_set '/'];
Output_name = [QC '_Loc' num2str(Loc) '_Mask_year' num2str(Mask_year) '.mat'];
load([Assim_dir Output_name])
filter_period = [period(1)-recon.sttime+1:period(1)-recon.sttime+simul_year];
IHCEnKF_prior_save = recon.priorvalue(:,filter_period);
IHCEnKF_poste_save = recon.postevalue(:,filter_period);
Inc_AOEnKF = IHCEnKF_poste_save-IHCEnKF_prior_save;

Assim_dir = ['../Data/Assim_output/AOEnKF_B/' Exp_set '/'];
Output_name = [QC '_Loc5' '_Mask_year' num2str(Mask_year) '.mat'];
load([Assim_dir Output_name])
filter_period = [period(1)-recon.sttime+1:period(1)-recon.sttime+simul_year];
IHCEnKF_prior_save = recon.priorvalue(:,filter_period);
IHCEnKF_poste_save = recon.postevalue(:,filter_period);
Inc_AOEnKF_B = IHCEnKF_poste_save-IHCEnKF_prior_save;

clear recon
Assim_dir = ['../Data/Assim_output/' Model_type '/' Exp_set '/'];
if ~exist(Assim_dir, 'dir')
    mkdir(Assim_dir)
end
% Linear combination of increments
ahps = [0:0.1:1];
iahp = 1;
ahp = ahps(iahp);
for iahp = 1:length(ahps)
    ahp = ahps(iahp);
    EnKF_prior_save = IHCEnKF_prior_save;
    EnKF_poste_save = IHCEnKF_prior_save+ahp.*Inc_AOEnKF+(1-ahp).*Inc_AOEnKF_B;

    Output_name = [QC '_Alpha' num2str(ahp, '%.1f') '_Mask_year' num2str(Mask_year) '.mat'];
    recon.priorvalue = EnKF_prior_save;
    recon.postevalue = EnKF_poste_save;
    recon.sttime = period(1);    
    save([Assim_dir Output_name],'recon')
end