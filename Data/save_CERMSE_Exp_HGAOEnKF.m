clear
addpath(genpath('./Utils'))
% clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
period = [1880:2000];
simul_year = length(period);

Exp_sets = {'REAL','OSSE'};
iExp_set = 2;
Exp_set = Exp_sets{iExp_set};

Model_types = {'HGAOEnKF'};
iModel_type = 1;
Model_type = Model_types{iModel_type};

% Exp_names = {'WiQC','NoQC'};
Exp_names = {'NoQC'};
iExp_name = 1;
Exp_name = Exp_names{iExp_name};

Mask_years = [1305 1945];
% Mask_years = [1945];
iMask_year = 1;
Mask_year = Mask_years(iMask_year);

ahps = [0.0:0.1:1.0];
iahp = 1;
ahp = ahps(iahp);

Obs_raw_names={'20CR','Truth'};
iObs_raw_name = iExp_set;
Obs_raw_name = Obs_raw_names{iObs_raw_name};

Obs_dir = ['../Data/obs/Obs_preprocess/' Obs_raw_name '/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Obs_name = ['obs.mat'];
load([Obs_dir Obs_name])
lat = obs.lat;
lon = obs.lon;
lev = obs.lev;
model_size = length(lat)*length(lon);

for imodel = 1:model_size
    mod_lat_all(imodel,1) = lat(fix((imodel-1)/length(lon))+1);
    mod_lon_all(imodel,1) = lon(mod(imodel-1,length(lon))+1);    
    wt(imodel,1) = cos(mod_lat_all(imodel,1)*pi/180);
end

for iModel_type = 1:length(Model_types)
    Model_type = Model_types{iModel_type};
    for iExp_name = 1:length(Exp_names)
        Exp_name = Exp_names{iExp_name};     
        CE.(Model_type).(Exp_name).prior_gmt = zeros(length(Mask_years),length(ahps),length(lev));
        CE.(Model_type).(Exp_name).poste_gmt = zeros(length(Mask_years),length(ahps),length(lev));
        CE.(Model_type).(Exp_name).prior_fedmean = zeros(length(Mask_years),length(ahps),length(lev));
        CE.(Model_type).(Exp_name).poste_fedmean = zeros(length(Mask_years),length(ahps),length(lev));
        CE.(Model_type).(Exp_name).prior_fed = zeros(length(Mask_years),length(ahps),length(lev),length(lon),length(lat));
        CE.(Model_type).(Exp_name).poste_fed = zeros(length(Mask_years),length(ahps),length(lev),length(lon),length(lat));
        RMSE.(Model_type).(Exp_name).prior_gmt = zeros(length(Mask_years),length(ahps),length(lev));
        RMSE.(Model_type).(Exp_name).poste_gmt = zeros(length(Mask_years),length(ahps),length(lev));
        RMSE.(Model_type).(Exp_name).prior_fedmean = zeros(length(Mask_years),length(ahps),length(lev));
        RMSE.(Model_type).(Exp_name).poste_fedmean = zeros(length(Mask_years),length(ahps),length(lev));
        RMSE.(Model_type).(Exp_name).prior_fed = zeros(length(Mask_years),length(ahps),length(lev),length(lon),length(lat));
        RMSE.(Model_type).(Exp_name).poste_fed = zeros(length(Mask_years),length(ahps),length(lev),length(lon),length(lat));      
        for iMask_year = 1:length(Mask_years)
            Mask_year = Mask_years(iMask_year);
            for iahp = 1:length(ahps)
                ahp = ahps(iahp);
                Assim_dir = ['../Data/Assim_output/' Model_type '/' Exp_set '/'];
                Output_name = [Exp_name '_Alpha' num2str(ahp,'%.1f') '_Mask_year' num2str(Mask_year) '.mat'];
                load([Assim_dir Output_name])
                for ilev = 1:length(lev)
                    filter_period = [period(1)-recon.sttime+1:period(1)-recon.sttime+simul_year];
                    IHCEnKF_prior_save = recon.priorvalue((ilev-1)*model_size+1:ilev*model_size,filter_period);
                    IHCEnKF_poste_save = recon.postevalue((ilev-1)*model_size+1:ilev*model_size,filter_period); 
                    filter_period = [period(1)-obs.sttime+1:period(1)-obs.sttime+simul_year];
                    Truth = obs.value((ilev-1)*model_size+1:ilev*model_size,filter_period);
                    U_mask = std(Truth,0,2)==0;
                    IHCEnKF_prior_save(U_mask,:) = nan;
                    IHCEnKF_poste_save(U_mask,:) = nan;
                    Truth(U_mask,:) = nan;
                    Truth_gmt = global_mean_nan(Truth,wt);
                
                    Prior_gmt = global_mean_nan(IHCEnKF_prior_save,wt);
                    Poste_gmt = global_mean_nan(IHCEnKF_poste_save,wt);
                    CE_Prior_fed = CEfield(IHCEnKF_prior_save,Truth);
                    CE_Poste_fed = CEfield(IHCEnKF_poste_save,Truth);
                    RMSE_Prior_fed = RMSEfield(IHCEnKF_prior_save,Truth);
                    RMSE_Poste_fed = RMSEfield(IHCEnKF_poste_save,Truth);

                    mask = ~isnan(CE_Prior_fed);
                    CE.(Model_type).(Exp_name).prior_gmt(iMask_year,iahp,ilev) = 1-sum((Prior_gmt-Truth_gmt).^2)/(length(Truth_gmt)*var(Truth_gmt));
                    CE.(Model_type).(Exp_name).poste_gmt(iMask_year,iahp,ilev) = 1-sum((Poste_gmt-Truth_gmt).^2)/(length(Truth_gmt)*var(Truth_gmt));                    
                    CE.(Model_type).(Exp_name).prior_fedmean(iMask_year,iahp,ilev) = sum(wt(mask).*CE_Prior_fed(mask))/sum(wt(mask));
                    CE.(Model_type).(Exp_name).prior_fed(iMask_year,iahp,ilev,:,:) = reshape(CE_Prior_fed,length(lon),length(lat));
                    CE.(Model_type).(Exp_name).poste_fedmean(iMask_year,iahp,ilev) = sum(wt(mask).*CE_Poste_fed(mask))/sum(wt(mask));
                    CE.(Model_type).(Exp_name).poste_fed(iMask_year,iahp,ilev,:,:) = reshape(CE_Poste_fed,length(lon),length(lat));

                    RMSE.(Model_type).(Exp_name).prior_gmt(iMask_year,iahp,ilev) = sqrt(sum((Prior_gmt-Truth_gmt).^2)/length(Truth_gmt));
                    RMSE.(Model_type).(Exp_name).poste_gmt(iMask_year,iahp,ilev) = sqrt(sum((Poste_gmt-Truth_gmt).^2)/length(Truth_gmt));            
                    RMSE.(Model_type).(Exp_name).prior_fedmean(iMask_year,iahp,ilev) = sum(wt(mask).*RMSE_Prior_fed(mask))/sum(wt(mask));
                    RMSE.(Model_type).(Exp_name).prior_fed(iMask_year,iahp,ilev,:,:) = reshape(RMSE_Prior_fed,length(lon),length(lat));                   
                    RMSE.(Model_type).(Exp_name).poste_fedmean(iMask_year,iahp,ilev) = sum(wt(mask).*RMSE_Poste_fed(mask))/sum(wt(mask));
                    RMSE.(Model_type).(Exp_name).poste_fed(iMask_year,iahp,ilev,:,:) = reshape(RMSE_Poste_fed,length(lon),length(lat));
                end % ilev
            end % iLoc
        end % imask
    end % iname
end % itype

save(['CERMSE_Exp_' Model_type '_' Exp_set '.mat'],'CE','RMSE')