clear
addpath(genpath('../../Utils'))
% clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
period = [851:2003];
simul_year = length(period);

Mask_year = 1305;

Prior_raw_names = {'CESM2'};
iPrior_raw_name = 1;
Prior_raw_name = Prior_raw_names{iPrior_raw_name};

Proxy_raw_names = {'PAGES2k','Pseudoproxy'};
iProxy_raw_name = 1;
Proxy_raw_name = Proxy_raw_names{iProxy_raw_name};

Obs_raw_names={'20CR','Truth'};
iObs_raw_name = 1;
Obs_raw_name = Obs_raw_names{iObs_raw_name};

Prior_dir = ['../../prior/Prior_preprocess/' Prior_raw_name '/'];
Proxy_dir = ['../../proxy/Proxy_preprocess/' Proxy_raw_name '/' num2str(Mask_year) '/'];
Obs_dir = ['../../obs/Obs_preprocess/' Obs_raw_name '/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess the prior
Prior_raw_dir = ['../Prior_raw/' Prior_raw_name '/'];
filename = ['regridded_prior_tas_tos.nc'];
lat = ncread([Prior_raw_dir,filename],'lat');
lon = ncread([Prior_raw_dir,filename],'lon');
lev = ncread([Prior_raw_dir,filename],'lev');
Norvar = ncread([Prior_raw_dir,filename],'tas_tos');
model_size = length(lon)*length(lat);
stds = ncread([Prior_raw_dir,filename],'stdtemp');
stds = squeeze(mean(stds,[1,2],'omitnan'));
stds_reshaped = repmat(reshape(stds,[1,1,length(lev),1]),[length(lon) length(lat) 1 size(Norvar,4)]);
Norvar = Norvar.*stds_reshaped;
Norvar = reshape(Norvar,length(lon)*length(lat)*length(lev),[]); % reshape as model_size*num_V x time

% Load the proxy
Proxy_name = ['proxy.mat'];
load([Proxy_dir Proxy_name])
Proxy_lat_all = proxy.Proxy_lat_all;
Proxy_lon_all = proxy.Proxy_lon_all;
Proxy_ols_all = proxy.Proxy_ols_all;
Proxy_sea_all = proxy.Proxy_sea_all;

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

Prior = Norvar(1:model_size,:);
Prior_seamean = zeros(model_size*length(sea_lst),simul_year);
for isea = 1:length(sea_lst)
    for iyear = 1:simul_year
        Prior_seamean((isea-1)*model_size+1:isea*model_size,iyear) = mean(Prior(:,iyear*12+sea_lst{isea}),2);
    end
end
intec = Proxy_ols_all(:,1);
slope = Proxy_ols_all(:,2);
Ye = (Prior_seamean(nearest_index,:).*kron(slope,ones(1,simul_year))+kron(intec,ones(1,simul_year)))'; 

% save prior
for iyear = 1:simul_year
    Prior_annual(:,iyear) = mean(Norvar(:,iyear*12+sea_lst{1}),2);
end

Prior_eartime = 851;
filter_period = period(1)-Prior_eartime+1:period(1)-Prior_eartime+simul_year;
prior.value = Prior_annual(:,filter_period);
prior.lat = lat;
prior.lon = lon;
prior.lev = lev;
prior.Ye = Ye;
prior.sttime = period(1);

if ~exist(Prior_dir,'dir')
    mkdir(Prior_dir)
end
Prior_name = ['prior.mat'];
save([[Prior_dir Prior_name]],'prior')


