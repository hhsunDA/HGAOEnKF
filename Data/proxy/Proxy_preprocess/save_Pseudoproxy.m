clear
addpath(genpath('../../Utils'))
% clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(1)
period = [851:2003];
simul_year = length(period);

Mask_year = 1945;

Prior_raw_names = {'CESM2'};
iPrior_raw_name = 1;
Prior_raw_name = Prior_raw_names{iPrior_raw_name};

Proxy_raw_names = {'PAGES2k','Pseudoproxy'};
iProxy_raw_name = 2;
Proxy_raw_name = Proxy_raw_names{iProxy_raw_name};

Obs_raw_names={'20CR','Truth'};
iObs_raw_name = 2;
Obs_raw_name = Obs_raw_names{iObs_raw_name};

Prior_dir = ['../../prior/Prior_preprocess/' Prior_raw_name '/'];
Proxy_dir = ['../../proxy/Proxy_preprocess/' Proxy_raw_name '/' num2str(Mask_year) '/'];
Obs_dir = ['../../obs/Obs_preprocess/' Obs_raw_name '/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess the truth
Obs_raw_dir = ['../../obs/Obs_raw/' Obs_raw_name '/'];
filename = ['regridded_truth_tas_tos.nc'];
lat = ncread([Obs_raw_dir,filename],'lat');
lon = ncread([Obs_raw_dir,filename],'lon');
lev = ncread([Obs_raw_dir,filename],'lev');
Norvar = ncread([Obs_raw_dir,filename],'tas_tos');
model_size = length(lon)*length(lat);
stds = ncread([Obs_raw_dir,filename],'stdtemp');
stds = squeeze(mean(stds,[1,2],'omitnan'));
stds_reshaped = repmat(reshape(stds,[1,1,length(lev),1]),[length(lon) length(lat) 1 size(Norvar,4)]);
Norvar = Norvar.*stds_reshaped;
Norvar = reshape(Norvar,length(lon)*length(lat)*length(lev),[]); % reshape as model_size*num_V x time

% Load the proxy
Proxy_name = ['proxy.mat'];
load(['../../proxy/Proxy_preprocess/PAGES2k/' num2str(Mask_year) '/' Proxy_name])
Proxy_lat_all = proxy.Proxy_lat_all;
Proxy_lon_all = proxy.Proxy_lon_all;
Proxy_error_var_all = proxy.Proxy_error_var_all;
Proxy_error_SNR_all = proxy.Proxy_error_SNR_all;
Proxy_ols_all = proxy.Proxy_ols_all;
Proxy_sea_all = proxy.Proxy_sea_all;
Ptype_all = proxy.Ptype_all;
filter_period = period(1)-proxy.sttime+1:period(1)-proxy.sttime+simul_year;
Proxy_lat = proxy.Proxy_lat(filter_period);
Proxy_lon = proxy.Proxy_lon(filter_period);
Proxy_R = proxy.Proxy_R(filter_period);
Proxy_idx = proxy.Proxy_idx(filter_period);
Proxy_SNR = proxy.Proxy_SNR(filter_period);
Proxy_ols = proxy.Proxy_ols(filter_period);
Proxy_sea = proxy.Proxy_sea(filter_period);
Psm_type = proxy.Psm_type(filter_period);
Ptype = proxy.Ptype(filter_period);

% find the index of the x that is used to calculate Hx
sea_lst = {[1:12],[6:8],[3:8],[6:11],[0:2],[-3:2],[0:5]};
dist = zeros(length(lon),length(lat),length(Proxy_lat_all));
for iobs = 1:length(Proxy_lat_all)
    for ilon = 1:length(lon)
        for ilat = 1:length(lat)
            dist(ilon,ilat,iobs) =  disHaversine(lon(ilon),lat(ilat),Proxy_lon_all(iobs),Proxy_lat_all(iobs));
        end
    end
end
dist = reshape(dist,length(lon)*length(lat),[]);
[~,index] = sort(dist,1);
nearest_index = index(1,:)';
nearest_index =  (Proxy_sea_all-1).*model_size + nearest_index;

Truth = Norvar(1:model_size,:);
Truth_seamean = zeros(model_size*length(sea_lst),simul_year);
for isea = 1:length(sea_lst)
    for iyear = 1:simul_year
        Truth_seamean((isea-1)*model_size+1:isea*model_size,iyear) = mean(Truth(:,iyear*12+sea_lst{isea}),2);
    end
end
for iyear = 1:simul_year
    Truth_annual(:,iyear) = mean(Norvar(:,iyear*12+sea_lst{1}),2);
end
Truth_eartime = 851;
filter_period = period(1)-Truth_eartime+1:period(1)-Truth_eartime+simul_year;
Truth_seamean = Truth_seamean(:,filter_period);
Truth_annual = Truth_annual(:,filter_period);
intec = Proxy_ols_all(:,1);
slope = Proxy_ols_all(:,2);
Proxy_true = Truth_seamean(nearest_index,:).*kron(slope,ones(1,simul_year))+kron(intec,ones(1,simul_year)); 

obs_error_var = Proxy_error_var_all;
noise = randn(size(Proxy_true,1),size(Proxy_true,2)).*sqrt(kron(obs_error_var,ones(1,size(Proxy_true,2))));
Proxy = Proxy_true + noise;
for iassim = 1:simul_year
    Proxy_tmp{iassim} = Proxy(Proxy_idx{iassim},iassim)';
end
Proxy = Proxy_tmp;

proxy.Proxy = Proxy;
proxy.Proxy_lat = Proxy_lat;
proxy.Proxy_lon = Proxy_lon;
proxy.Proxy_R = Proxy_R;
proxy.Proxy_idx = Proxy_idx;
proxy.Proxy_SNR = Proxy_SNR;
proxy.Proxy_ols = Proxy_ols;
proxy.Proxy_sea = Proxy_sea;
proxy.Psm_type = Psm_type;
proxy.Ptype = Ptype;
proxy.sttime = period(1);
proxy.Proxy_lon_all = Proxy_lon_all;
proxy.Proxy_lat_all = Proxy_lat_all;
proxy.Proxy_error_var_all = Proxy_error_var_all;
proxy.Proxy_error_SNR_all = Proxy_error_SNR_all;
proxy.Proxy_ols_all = Proxy_ols_all;
proxy.Proxy_sea_all = Proxy_sea_all;
proxy.Ptype_all = Ptype_all;

if ~exist(Proxy_dir,'dir')
    mkdir(Proxy_dir)
end
Proxy_name = ['proxy.mat'];
save([Proxy_dir Proxy_name],'proxy')

obs.value = Truth_annual;
obs.lat = lat;
obs.lon = lon;
obs.lev = lev;
obs.sttime = period(1);

if ~exist(Obs_dir,'dir')
    mkdir(Obs_dir)
end
Obs_name = ['obs.mat'];
save([[Obs_dir Obs_name]],'obs')