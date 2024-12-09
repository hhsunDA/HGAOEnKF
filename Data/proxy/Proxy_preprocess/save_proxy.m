clear

% Define the earliest, latest time, and proxy network
filter_timemask = true;
period = [0:2003];
simul_year = length(period);

filter_spacemask = true;
Mask_year = 1945; % Define the proxy network according to the time of the year

% Read the raw proxy one by one
Proxy_raw_names = {'PAGES2k'};
iProxy_raw_name = 1;
Proxy_raw_name = Proxy_raw_names{iProxy_raw_name};
Proxy_raw_dir = ['../Proxy_raw/' Proxy_raw_name '/'];
filelist = dir([Proxy_raw_dir, '*.nc']);
sealist = {'1_2_3_4_5_6_7_8_9_10_11_12','6_7_8','3_4_5_6_7_8','6_7_8_9_10_11','-12_1_2','-9_-10_-11_-12_1_2','-12_1_2_3_4_5'};
for iproxy = 1:length(filelist)
    filename = [Proxy_raw_dir,filelist(iproxy).name];
    finfo = ncinfo(filename);
    Psm_type_all{iproxy} = finfo.Attributes(1).Value;
    [~, Proxy_sea_all(iproxy)] = ismember(finfo.Attributes(2).Value,sealist);
    Ptype_all{iproxy} = finfo.Attributes(3).Value;
    Proxy_all{iproxy} = ncread(filename,'value');
    Proxy_time{iproxy} = ncread(filename,'time');
    Proxy_lat_all(iproxy) = ncread(filename,'lat');
    Proxy_lon_all(iproxy) = ncread(filename,'lon');
    if(Proxy_lon_all(iproxy) > 180)
        Proxy_lon_all(iproxy) = Proxy_lon_all(iproxy)-360;  
    else
        Proxy_lon_all(iproxy) = Proxy_lon_all(iproxy);
    end
    Proxy_error_var_all(iproxy) = ncread(filename,'R');
    Proxy_error_SNR_all(iproxy) = ncread(filename,'SNR');
    Proxy_ols_all(iproxy,:) = ncread(filename,'ols');
end

% Find the earliest and latest time of proxy from the raw data
if( ~ filter_timemask)
    period(1) = inf;
    for iproxy = 1:length(Proxy_time)
        if(period(1) >= Proxy_time{iproxy}(1))
            period(1) = Proxy_time{iproxy}(1);
        end
    end
end

% Obtain each year's full proxy network
count_dup = 0; % DR99ABR01_d18O has 2 proxies in 1959 years;
for iassim = 1:simul_year
    count_P = 1;
    for iproxy = 1:length(filelist)
        Proxy_time_idx = find(Proxy_time{iproxy} == period(iassim));
        if( ~ isempty(Proxy_time_idx))
            if(length(Proxy_time_idx) > 1)
                count_dup = count_dup+1;
            end
            Proxy_filter_sy(count_P) = Proxy_all{iproxy}(Proxy_time_idx(1));
            Proxy_idx_sy(count_P) = iproxy;
            Proxy_lat_sy(count_P) = Proxy_lat_all(iproxy);
            Proxy_lon_sy(count_P) = Proxy_lon_all(iproxy);
            Proxy_R_sy(count_P) = Proxy_error_var_all(iproxy);
            Proxy_SNR_sy(count_P) = Proxy_error_SNR_all(iproxy);
            Psm_type_sy{count_P} = Psm_type_all{iproxy};        
            Proxy_ols_sy(:,count_P) = Proxy_ols_all(iproxy,:);
            Proxy_sea_sy(count_P) = Proxy_sea_all(iproxy); 
            Ptype_sy{count_P} = Ptype_all{iproxy};
            count_P = count_P+1; 
        end
    end
    Proxy{iassim} = Proxy_filter_sy;
    Proxy_idx{iassim} = Proxy_idx_sy;
    Proxy_lat{iassim} = Proxy_lat_sy;
    Proxy_lon{iassim} = Proxy_lon_sy;
    Proxy_R{iassim} = Proxy_R_sy;
    Proxy_SNR{iassim} = Proxy_SNR_sy;
    Psm_type{iassim} = Psm_type_sy;
    Proxy_ols{iassim} = Proxy_ols_sy;
    Proxy_sea{iassim} = Proxy_sea_sy;
    Ptype{iassim} = Ptype_sy;
    clear Proxy_filter_sy Proxy_idx_sy Proxy_lat_sy Proxy_lon_sy Proxy_R_sy Proxy_SNR_sy Psm_type_sy Proxy_ols_sy Proxy_sea_sy Ptype_sy% size of _sy changes every time
end

% Find the year that the proxy network is the densest
if( ~ filter_spacemask)
    Mask_year = period(1);
    Mask = Proxy_idx{Mask_year-period(1)+1};
    for iassim = 1:simul_year
        if(length(Mask) <= length(Proxy{iassim}))
            Mask = Proxy_idx{iassim};
            Mask_year = period(iassim);
        end
    end
end
  
% Filter the proxy network with new space mask, space mask is determined by
% the time of the year
Mask = Proxy_idx{Mask_year-period(1)+1};
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
    Tmp = Psm_type{iassim};
    Psm_type{iassim} = Tmp(proxy_filter{iassim});
    Tmp = Proxy_ols{iassim};
    Proxy_ols{iassim} = Tmp(proxy_filter{iassim});
    Tmp = Proxy_sea{iassim};
    Proxy_sea{iassim} = Tmp(proxy_filter{iassim});        
    Tmp = Ptype{iassim};
    Ptype{iassim} = Tmp(proxy_filter{iassim});
end    

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

Proxy_dir = ['./' Proxy_raw_name '/' num2str(Mask_year) '/'];
if ~exist(Proxy_dir,'dir')
    mkdir(Proxy_dir)
end
Proxy_name = ['proxy.mat'];
save([Proxy_dir Proxy_name],'proxy')