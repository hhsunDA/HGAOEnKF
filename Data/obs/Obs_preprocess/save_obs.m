clear
% clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
period = [1880:2000];
simul_year = length(period);

Obs_raw_names = {'BE','20CR','GISTEMP'};
iObs_raw_name = 2;
Obs_raw_name = Obs_raw_names{iObs_raw_name};

Obs_dir = ['../../obs/Obs_preprocess/' Obs_raw_name '/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess the obs
Obs_raw_dir = ['../Obs_raw/' Obs_raw_name '/'];
filename = ['regridded_test_tas_tos.nc'];

lat = ncread([Obs_raw_dir,filename],'lat');
lon = ncread([Obs_raw_dir,filename],'lon');
lev = ncread([Obs_raw_dir,filename],'lev');
Norvar = ncread([Obs_raw_dir,filename],'tas_tos');
model_size = length(lon)*length(lat);
Norvar = reshape(Norvar,length(lon)*length(lat)*length(lev),[]); % reshape as model_size*num_V x time

Obs_annual = squeeze(mean(reshape(Norvar,length(lon)*length(lat)*length(lev),12,[]),2));

Obs_eartime = 1880;
filter_period = period(1)-Obs_eartime+1:period(1)-Obs_eartime+simul_year;
obs.value = Obs_annual(:,filter_period);
obs.lat = lat;
obs.lon = lon;
obs.lev = lev;
obs.sttime = period(1);

if ~exist(Obs_dir,'dir')
    mkdir(Obs_dir)
end
Obs_name = ['obs.mat'];
save([[Obs_dir Obs_name]],'obs')