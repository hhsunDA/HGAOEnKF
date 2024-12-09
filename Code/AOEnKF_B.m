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

Model_types = {'AOEnKF_B'};
iModel_type = 1;
Model_type = Model_types{iModel_type};

QCs = {'WiQC','NoQC'};
iQC = 2;
QC = QCs{iQC};

serial_update = false;

localize = false;
% localization_values = [earth_R*pi/4;earth_R*pi/2;earth_R*pi;earth_R*pi*2]; % unit:km
localization_values = [earth_R*pi/2]; % unit:km
num_localization = length(localization_values);

Prior_raw_names = {'CESM2'};
iPrior_raw_name = 1;
Prior_raw_name = Prior_raw_names{iPrior_raw_name};

Proxy_raw_names = {'PAGES2k','Pseudoproxy'};
iProxy_raw_name = iExp_set;
Proxy_raw_name = Proxy_raw_names{iProxy_raw_name};
 
Prior_dir = ['../Data/prior/Prior_preprocess/' Prior_raw_name '/'];
Proxy_dir = ['../Data/proxy/Proxy_preprocess/' Proxy_raw_name '/' num2str(Mask_year) '/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the prior
Prior_name = ['prior.mat'];
load([Prior_dir Prior_name])
lat = prior.lat;
lon = prior.lon;
lev = prior.lev;
model_size=length(lat)*length(lon);
Model_climo = prior.value(:,1:climo_size);
Hx_climo = prior.Ye(1:climo_size,:)';
Hx_size = size(Hx_climo,1);
Model_climo = cat(1, Model_climo, Hx_climo);

for imodel = 1:model_size
    mod_lat_all(imodel,1) = lat(fix((imodel-1)/length(lon))+1);
    mod_lon_all(imodel,1) = lon(mod(imodel-1,length(lon))+1);
end

for ilev = 1:length(lev)-1
    mod_lat_all = cat(1,mod_lat_all,mod_lat_all);
    mod_lon_all = cat(1,mod_lon_all,mod_lon_all);
end

% Load the proxy
Proxy_name = ['proxy.mat'];
load([Proxy_dir Proxy_name])

proxy_lat = proxy.Proxy_lat_all;
proxy_lon = proxy.Proxy_lon_all;

% Filter the time of proxy
filter_period = period(1)-proxy.sttime+1:period(1)-proxy.sttime+simul_year;
Proxy = proxy.Proxy(filter_period);
Proxy_R = proxy.Proxy_R(filter_period);
Proxy_SNR = proxy.Proxy_SNR(filter_period);
Proxy_idx = proxy.Proxy_idx(filter_period);
Proxy_lat = proxy.Proxy_lat(filter_period);
Proxy_lon = proxy.Proxy_lon(filter_period);

% Load the analog index
ANAind_name = 'ANA_ind.mat';
load([Proxy_dir ANAind_name])
filter_period = period(1)-ANA_ind.sttime+1:period(1)-ANA_ind.sttime+simul_year;
index_rms = ANA_ind.RMSEind(filter_period,:);

% cat the lat and lon of the prior and proxy, which is needed by GC localization
obs_lon_all = proxy_lon';
obs_lat_all = proxy_lat';
obs_grids_all = [obs_lon_all obs_lat_all];
update_lon_all = [mod_lon_all;obs_lon_all]; 
update_lat_all = [mod_lat_all;obs_lat_all];
update_grids_all = [update_lon_all update_lat_all];

% AOEnKF_B assimilation
EnKF_prior_save = zeros(model_size*length(lev),simul_year);
EnKF_poste_save = zeros(model_size*length(lev),simul_year);
Assim_dir = ['../Data/Assim_output/' Model_type '/' Exp_set '/'];
if ~exist(Assim_dir, 'dir')
    mkdir(Assim_dir)
end
for iloc = 1:num_localization
    localization_value = localization_values(iloc);
    CMat_all = GC(obs_grids_all,update_grids_all,localization_value); 
    xb_prime = Model_climo - mean(Model_climo,2); 
    for iassim = 1:simul_year
        obs_error_var = Proxy_R{iassim};
        zobs_total = Proxy{iassim};
        proxy_idx = Proxy_idx{iassim};
        CMat = CMat_all(proxy_idx,:); 
    
        if (strcmp(QC,'WiQC'))
            Hxens = Hx_climo(proxy_idx,:)';
            O = mean(Hxens);
            OMB = abs(O-zobs_total);
            stdR = sqrt(obs_error_var);
            assim_ind = find(OMB > stdR);
            obs_error_var(assim_ind) = obs_error_var(assim_ind).*(OMB(assim_ind)./stdR(assim_ind));  
        end
    
        ind_sample_year = index_rms(iassim,1:ensemble_size);
        Model_sample = Model_climo(:,ind_sample_year);
    
        xb_mean = mean(Model_sample,2);
        xb = xb_mean + xb_prime;
        EnKF_prior_save(:,iassim) = xb_mean(1:end-Hx_size);
        zobs = zobs_total';
        tic;
        xa = EnSRF(xb,obs_error_var,zobs,Hx_size,proxy_idx,CMat,localize,serial_update);
        elapsedTime(iassim) = toc;
        xa_mean = mean(xa,2);
        xa_ana_prime = xa(:,ind_sample_year) - mean(xa(:,ind_sample_year),2);
        xa = xa_mean + xa_ana_prime;
        EnKF_poste_save(:,iassim) = xa_mean(1:end-Hx_size);
    end % iassim
    Output_name = [QC '_Loc5_Mask_year' num2str(Mask_year) '.mat'];
    recon.priorvalue = EnKF_prior_save;
    recon.postevalue = EnKF_poste_save;
    recon.sttime = period(1);
    save([Assim_dir Output_name],'recon')
end

