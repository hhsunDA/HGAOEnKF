clear
addpath(genpath('../../Utils'))
% clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% period = [851:2003]; % for Pseudoproxy, min of period = 851; for PAGES2k, minof period = 0
period = [0:2003]; % for Pseudoproxy, min of period = 851; for PAGES2k, minof period = 0
simul_year = length(period);
ensemble_size = 100;
climo_size = 1000;
Mask_year = 1305;
isosparse = true; % select the analogs based on isotropic sparse regions: Every model grid/box has at most 1 proxy for selecting 
co_length = 10; % select based on the continuous simulation years of the climate model

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
% Load the prior
Prior_name = ['prior.mat'];
load([Prior_dir Prior_name])
lat = prior.lat;
lon = prior.lon;
lev = prior.lev;
model_size = length(lat)*length(lon);
Model_climo = prior.value(1:length(lev)*model_size,1:climo_size);
Hx_climo = prior.Ye(1:climo_size,:)';
for imodel = 1:model_size
    mod_lat_all(imodel,1) = lat(fix((imodel-1)/length(lon))+1);
    mod_lon_all(imodel,1) = lon(mod(imodel-1,length(lon))+1);    
    wt(imodel,1) = cos(mod_lat_all(imodel,1)*pi/180);
end
% Load the proxy
Proxy_name = ['proxy.mat'];
load([Proxy_dir Proxy_name])
proxy_lat = proxy.Proxy_lat_all;
proxy_lon = proxy.Proxy_lon_all;

filter_period = period(1)-proxy.sttime+1:period(1)-proxy.sttime+simul_year;
Proxy = proxy.Proxy(filter_period);
Proxy_R = proxy.Proxy_R(filter_period);
Proxy_SNR = proxy.Proxy_SNR(filter_period);
Proxy_idx = proxy.Proxy_idx(filter_period);
Proxy_lat = proxy.Proxy_lat(filter_period);
Proxy_lon = proxy.Proxy_lon(filter_period);

% Isotropic box
if(isosparse) 
    % Define isotropic sparse grids
    lat_interval = 2;
    lon_interval = 2;
    lat_sparse = lat(lat_interval:lat_interval:length(lat));
    lon_sparse = lon(lon_interval:lon_interval:length(lon));
    for iproxy = 1:length(proxy_lat)
        [lat1,lat2,lon1,lon2] = return_boxindex(iproxy,lat_sparse,lon_sparse,proxy_lat,proxy_lon);  
        obs_flag(iproxy,1) = lat1;
        obs_flag(iproxy,2) = lat2;
        obs_flag(iproxy,3) = lon1;
        obs_flag(iproxy,4) = lon2;    
    end
    
    p_lon_unique(1) = proxy_lon(1);
    p_lat_unique(1) = proxy_lat(1);
    p_obs_flag(1,:) = obs_flag(1,:);
    Mask(1) = 1;
    number = 1;
    for iproxy = 1:length(proxy_lat)
        for j = 1:length(p_lon_unique)
            flagA = p_obs_flag(j,:);
            flagB = obs_flag(iproxy,:);
            if(isequal(flagA,flagB))
                logic = 1;        
            end    
        end
        if(logic ~= 1)
            number = number+1;
            p_lon_unique(number) = proxy_lon(iproxy);
            p_lat_unique(number) = proxy_lat(iproxy);
            p_obs_flag(number,:) = obs_flag(iproxy,:);
            Mask(number) = iproxy;
        end
        logic = 0;   
    end 
    
    % just filter the proxy network with new mask
    for iassim = 1:simul_year
        count_P = 1;
        filter = [];
        for jj = 1:length(Mask)
            Tmp = find(Proxy_idx{iassim} == Mask(jj));
            if(~isempty(Tmp))
                filter(count_P) = Tmp;
                count_P = count_P+1;
            end
        end
        proxy_filter{iassim} = filter;
    end
    
    for iassim = 1:simul_year
        Tmp = Proxy_idx{iassim};
        Proxy_idx{iassim} = Tmp(proxy_filter{iassim});
        Tmp = Proxy_R{iassim};
        Proxy_R{iassim} = Tmp(proxy_filter{iassim});
        Tmp = Proxy{iassim};
        Proxy{iassim} = Tmp(proxy_filter{iassim});
        Tmp = Proxy_lat{iassim};
        Proxy_lat{iassim} = Tmp(proxy_filter{iassim});
        Tmp = Proxy_lon{iassim};
        Proxy_lon{iassim} = Tmp(proxy_filter{iassim});
        Tmp = Proxy_SNR{iassim};
        Proxy_SNR{iassim} = Tmp(proxy_filter{iassim});
    end
end

% Select the analog ensemble
cycle = floor(simul_year/co_length);
for ii = 1:climo_size+1-co_length
    Hx_climo_timeseries(:,:,ii) = Hx_climo(:,ii:ii+co_length-1); % choose the continuous time of Ye as sample pools.
end

for i = 1:cycle    
    % AOEnKF_E
    for ii = 1:climo_size+1-co_length
        for jj = 1:co_length
            proxy = Proxy{(i-1)*co_length+jj};
            proxy_lat = Proxy_lat{(i-1)*co_length+jj};
            proxy_idx = Proxy_idx{(i-1)*co_length+jj};
            proxy_R = Proxy_R{(i-1)*co_length+jj};
            proxy_wt = cos(proxy_lat*pi/180);
            % lowlatind = find(proxy_wt > 0.9);
            % proxy_wt(lowlatind) = 1;           
            xf_diff_all = (proxy-squeeze(Hx_climo_timeseries(proxy_idx,jj,ii))')./sqrt(proxy_R).*proxy_wt;
            xf_diff(ii,jj) = sqrt(sum(xf_diff_all.^2));
        end
    end
    [~,index] = sort(sum(xf_diff,2));    
    for jj = 1:co_length
        index_rms((i-1)*co_length+jj,:) = index+jj-1;
    end
end

i = cycle+1;
    % AOEnKF_E
    for ii = 1:climo_size+1-co_length
        for jj = 1:simul_year-cycle*co_length
            proxy = Proxy{(i-1)*co_length+jj};
            proxy_lat = Proxy_lat{(i-1)*co_length+jj};
            proxy_idx = Proxy_idx{(i-1)*co_length+jj};
            proxy_R = Proxy_R{(i-1)*co_length+jj};
            proxy_wt = cos(proxy_lat*pi/180);
            % lowlatind = find(proxy_wt > 0.9);
            % proxy_wt(lowlatind) = 1;            
            xf_diff_all = (proxy-squeeze(Hx_climo_timeseries(proxy_idx,jj,ii))')./sqrt(proxy_R).*proxy_wt;
            xf_diff(ii,jj) = sqrt(sum(xf_diff_all.^2));
        end
    end
    [~,index] = sort(sum(xf_diff,2));    
    for jj = 1:1:simul_year-cycle*co_length
        index_rms((i-1)*co_length+jj,:) = index+jj-1;
    end
ANA_ind.RMSEind = index_rms;
ANA_ind.sttime = period(1);
ANAind_name = ['ANA_ind.mat'];
save([Proxy_dir ANAind_name],'ANA_ind')

% Validation in instrumental period
vf_period = [1880:2000];
vf_year = length(vf_period);

filter_period = vf_period(1)-ANA_ind.sttime+1:vf_period(1)-ANA_ind.sttime+vf_year;
index_rms = index_rms(filter_period,:);

Obs_name = ['obs.mat'];
load([Obs_dir Obs_name])
Obs = obs.value;
filter_period = vf_period(1)-obs.sttime+1:vf_period(1)-obs.sttime+vf_year;
Obs = Obs(1:length(lev)*model_size,filter_period);

rng(1)
ind = randperm(climo_size,ensemble_size);
for iassim = 1:vf_year  
    % OEnKF
    zens = Model_climo(:,ind)';
    OFF_prior(:,iassim) = mean(zens);   
    % AOEnKF_E
    zens = Model_climo(:,index_rms(iassim,1:ensemble_size))';
    ANA_prior(:,iassim) = mean(zens);  
end

U_mask = ~isnan(Obs)&~isnan(OFF_prior);

for ilev = 1:length(lev)
    for iassim = 1:vf_year
        mask = U_mask((ilev-1)*model_size+1:ilev*model_size,iassim);
        diff = squeeze(ANA_prior((ilev-1)*model_size+1:ilev*model_size,iassim))-Obs((ilev-1)*model_size+1:ilev*model_size,iassim);
        error_ANA_prior(iassim,ilev) = sqrt(sum((wt(mask).*(diff(mask).^2)))/sum(wt(mask)));

        diff = squeeze(OFF_prior((ilev-1)*model_size+1:ilev*model_size,iassim))-Obs((ilev-1)*model_size+1:ilev*model_size,iassim);
        error_OFF_prior(iassim,ilev) = sqrt(sum((wt(mask).*(diff(mask).^2)))/sum(wt(mask)));
    end 
end
plot(error_ANA_prior-error_OFF_prior)
mean(error_ANA_prior-error_OFF_prior)

plt_ANA = reshape(ANA_prior,length(lon),length(lat),length(lev),[]);
plt_Obs = reshape(Obs,length(lon),length(lat),length(lev),[]);

year = 1910;
plt_year = year-vf_period(1)+1;
proxy_year = year-period(1)+1;
proxy_lat = Proxy_lat{proxy_year};
proxy_lon = Proxy_lon{proxy_year};
proxy_SNR = Proxy_SNR{proxy_year};
for iproxy = 1:length(proxy_lon)
    if(proxy_lon(iproxy) < 0)
        proxy_lon(iproxy) = proxy_lon(iproxy)+360;
    else
        proxy_lon(iproxy) = proxy_lon(iproxy);    
    end 
end

Figure_dir = ['../../Figures/'];

ilev = 1;
ifig = 100;
figure(ifig)
lgd_color_txt = 'color_b_w_r.mat';plt_proxy = true;left_lgd = -2;right_lgd = 2;label_txt = 'K';
title_txt = [Obs_raw_name];field = plt_Obs(:,:,ilev,plt_year);
plt_field
set(gca,'fontsize',16);
Figure_name = [Obs_raw_name '_' num2str(year) '.png'];
print(['-f' num2str(ifig)],'-dpng',[Figure_dir Figure_name])

ifig = 101;
figure(ifig)
lgd_color_txt = 'color_b_w_r.mat';plt_proxy = true;left_lgd = -2;right_lgd = 2;label_txt = 'K';
title_txt = ['AOEnKF Prior'];field = plt_ANA(:,:,ilev,plt_year);
plt_field
set(gca,'fontsize',16);
Figure_name = [Prior_raw_name '_' num2str(year) '_' num2str(Mask_year) '.png'];
print(['-f' num2str(ifig)],'-dpng',[Figure_dir Figure_name])