clear
addpath(genpath('../Data/Utils'))
% clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
year = 1999;
earth_R = 6367;
period = [year:year];
simul_year = length(period);
ensemble_size = 40;
climo_size = 1000;
Mask_year = 1945;

Exp_sets = {'REAL','OSSE'};
iExp_set = 1;
Exp_set = Exp_sets{iExp_set};

Model_types = {'AOEnKF'};
iModel_type = 1;
Model_type = Model_types{iModel_type};

QCs = {'WiQC','NoQC'};
iQC = 2;
QC = QCs{iQC};

serial_update = true;

localize = false;
localization_values = [earth_R*pi]; % unit:km
num_localization = length(localization_values);

Prior_raw_names = {'CESM2'};
iPrior_raw_name = 1;
Prior_raw_name = Prior_raw_names{iPrior_raw_name};

Proxy_raw_names = {'PAGES2k','Pseudoproxy'};
iProxy_raw_name = iExp_set;
Proxy_raw_name = Proxy_raw_names{iProxy_raw_name};

Obs_raw_names={'20CR','Truth'};
iObs_raw_name = 1;
Obs_raw_name = Obs_raw_names{iObs_raw_name};

Prior_dir = ['../Data/prior/Prior_preprocess/' Prior_raw_name '/'];
Proxy_dir = ['../Data/proxy/Proxy_preprocess/' Proxy_raw_name '/' num2str(Mask_year) '/'];
Obs_dir = ['../Data/obs/Obs_preprocess/' Obs_raw_name '/'];
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

Proxy_lat_all = proxy.Proxy_lat_all;
Proxy_lon_all = proxy.Proxy_lon_all;
Proxy_ols_all = proxy.Proxy_ols_all;
Ptype_all = proxy.Ptype_all;

% Filter the time of proxy
filter_period = period(1)-proxy.sttime+1:period(1)-proxy.sttime+simul_year;
Proxy = proxy.Proxy(filter_period);
Proxy_R = proxy.Proxy_R(filter_period);
Proxy_SNR = proxy.Proxy_SNR(filter_period);
Proxy_idx = proxy.Proxy_idx(filter_period);
Proxy_lat = proxy.Proxy_lat(filter_period);
Proxy_lon = proxy.Proxy_lon(filter_period);
Proxy_sea = proxy.Proxy_sea(filter_period);
Ptype = proxy.Ptype(filter_period);

filter_space = 50;
% filter_space = find(cellfun(@(x) contains(x, 'coral'), Ptype{1}));
% filter_space = filter_space(6);

% Load the obs
Obs_name = ['obs.mat'];
load([Obs_dir Obs_name])
Obs = obs.value;
filter_period = period(1)-obs.sttime+1:period(1)-obs.sttime+simul_year;
Obs = Obs(:,filter_period);

% Load the analog index
ANAind_name = 'ANA_ind.mat';
load([Proxy_dir ANAind_name])
filter_period = period(1)-ANA_ind.sttime+1:period(1)-ANA_ind.sttime+simul_year;
index_rms = ANA_ind.RMSEind(filter_period,:);

% cat the lat and lon of the prior and proxy, which is needed by GC localization
update_lon_all = [mod_lon_all; Proxy_lon_all]; 
update_lat_all = [mod_lat_all; Proxy_lat_all];
update_grids_all = [update_lon_all update_lat_all];
Proxy_grids_all = [Proxy_lon_all Proxy_lat_all];

% AOEnKF assimilation
EnKF_prior_save = zeros(model_size*length(lev),simul_year);
EnKF_poste_save = zeros(model_size*length(lev),simul_year);
Assim_dir = ['../Data/Assim_output/' Model_type '/' Exp_set '/'];
if ~exist(Assim_dir, 'dir')
    mkdir(Assim_dir)
end
for iloc = 1:num_localization
    localization_value = localization_values(iloc);
    CMat_all = GC(Proxy_grids_all,update_grids_all,localization_value); 
    for iassim = 1:simul_year
        obs_error_var = Proxy_R{iassim}(filter_space);
        zobs_total = Proxy{iassim}(filter_space);
        proxy_idx = Proxy_idx{iassim}(filter_space);
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

        xb = Model_sample;
        EnKF_prior_save(:,iassim) = mean(xb(1:end-Hx_size,:),2);
        zobs = zobs_total';
        tic;
        xa = EnSRF(xb,obs_error_var,zobs,Hx_size,proxy_idx,CMat,localize,serial_update);
        elapsedTime(iassim) = toc;
        EnKF_poste_save(:,iassim) = mean(xa(1:end-Hx_size,:),2);
    end % iassim
    Output_name = [QC '_Loc' num2str(iloc) '_Mask_year' num2str(Mask_year) '_stepbystep.mat'];
    recon.priorvalue = EnKF_prior_save;
    recon.postevalue = EnKF_poste_save;
    recon.sttime = period(1);
    % save([Assim_dir Output_name],'recon')
end % iloc

%% Plot step by step
Figure_dir = ['../Data/Figures/'];


plt_year = year-period(1)+1;
proxy_year = year-period(1)+1;
proxy_lat = Proxy_lat{proxy_year}(filter_space);
proxy_lon = Proxy_lon{proxy_year}(filter_space);
proxy_SNR = Proxy_SNR{proxy_year}(filter_space);
for iproxy = 1:length(proxy_lon)
    if(proxy_lon(iproxy) < 0)
        proxy_lon(iproxy) = proxy_lon(iproxy)+360;
    else
        proxy_lon(iproxy) = proxy_lon(iproxy);    
    end 
end

plt_Inc = reshape(EnKF_poste_save-EnKF_prior_save,length(lon),length(lat),length(lev),[]);
plt_ANA = reshape(EnKF_prior_save,length(lon),length(lat),length(lev),[]);
plt_Obs = reshape(Obs,length(lon),length(lat),length(lev),[]);

% linear regression between Hx and x
for ilon = 1:length(lon)
    for ilat = 1:length(lat)
        dist(ilon,ilat) =  disHaversine(lon(ilon),lat(ilat),proxy_lon,proxy_lat);
    end
end
dist = reshape(dist,length(lon)*length(lat),[]);
[~,index] = sort(dist,1);

x_idx = 7000;
x_lat = mod_lat_all(x_idx);
x_lon = mod_lon_all(x_idx);

Hxb_climo = Model_climo(end-Hx_size+1:end,:);
Hxb_sample = xb(end-Hx_size+1:end,:);
Hxa_sample = xa(end-Hx_size+1:end,:);
xb_sample1 = xb(index(1),:);
xb_sample2 = xb(x_idx,:);
xa_sample1 = xa(index(1),:);
xa_sample2 = xa(x_idx,:);

handles.ens_size = ensemble_size;
handles.climo_size = climo_size;
handles.obs = zobs;
handles.obs_error_sd = sqrt(obs_error_var);
handles.Hxb_climo = Hxb_climo(proxy_idx,:);
handles.Hxb_ens = Hxb_sample(proxy_idx,:);
handles.Hxa_ens = Hxa_sample(proxy_idx,:);
handles.xb_ens = cat(1,xb_sample1, xb_sample2);
handles.xa_ens = cat(1,xa_sample1, xa_sample2);
handles.ptype = Ptype_all(proxy_idx);

location = 2; % 1: locate at proxy site, 2: locate at moderate distance
prior_ens(1, :) = handles.Hxb_ens;
prior_ens(2, :) = handles.xb_ens(location,:);
poste_ens(1, :) = handles.Hxa_ens;
poste_ens(2, :) = handles.xa_ens(location,:);
ptype = handles.ptype;
% The best fit line on the ensemble
prior_mean = mean(prior_ens, 2);
prior_cov  = cov(prior_ens(1, :), prior_ens(2, :));
slope      = prior_cov(1, 2) / var(prior_ens(1, :));
intercept  = prior_mean(2) - slope * prior_mean(1);

best_x = [-1000 1000];
best_y = slope * best_x + intercept;

xlower = handles.obs - 4*handles.obs_error_sd;
xupper = handles.obs + 4*handles.obs_error_sd;
ylower = 0;
yupper = 7.0;
plt_num = 6;
plt_point = [floor(handles.ens_size/plt_num):floor(handles.ens_size/plt_num):handles.ens_size];
prior_ens = prior_ens(:,plt_point);
poste_ens = poste_ens(:,plt_point);

ilev = 1;
ifig = 100;
figure(ifig)
lgd_color_txt = 'color_b_w_r.mat';plt_proxy = false;left_lgd = -2;right_lgd = 2;label_txt = 'K';
title_txt = [Obs_raw_name];field = plt_Obs(:,:,ilev,plt_year);
plt_field_Hx_x
set(gca,'fontsize',16);
Figure_name = [Obs_raw_name '_' num2str(year) '.png'];
print(['-f' num2str(ifig)],'-dpng',[Figure_dir Figure_name])
ifig = ifig +1;

figure(ifig)
lgd_color_txt = 'color_b_w_r.mat';plt_proxy = true;left_lgd = -2;right_lgd = 2;label_txt = 'K';
title_txt = ['AOEnKF prior'];field = plt_ANA(:,:,ilev,plt_year);
plt_field_Hx_x
set(gca,'fontsize',16);
Figure_name = ['ANA' '_' num2str(year) '_' num2str(Mask_year) '.png'];
print(['-f' num2str(ifig)],'-dpng',[Figure_dir Figure_name])
ifig = ifig +1;

figure(ifig)
lgd_color_txt = 'color_b_w_r.mat';plt_proxy = true;left_lgd = -0.5;right_lgd = 0.5;label_txt = 'K';
title_txt = ['Increment'];field = plt_Inc(:,:,ilev,plt_year);
plt_field_Hx_x
set(gca,'fontsize',16);
Figure_name = ['Inc' '_' num2str(year) '_' num2str(Mask_year) '.png'];
print(['-f' num2str(ifig)],'-dpng',[Figure_dir Figure_name])
ifig = ifig +1;

figure(ifig)
lgd_color_txt = 'color_b_w_r.mat';plt_proxy = true;left_lgd = -0.1;right_lgd = 0.1;label_txt = 'K';
title_txt = [''];field = zeros(length(lon),length(lat));
plt_field_Hx_x
set(gca,'fontsize',16);
Figure_name = ['Proxy' '_' num2str(year) '_' num2str(Mask_year) '.png'];
print(['-f' num2str(ifig)],'-dpng',[Figure_dir Figure_name])
ifig = ifig +1;

atts = stylesheet; % get the default fonts and colors
figXmin   = 450; 
figYmin   = 250; 
figWidth  = 670; 
figHeight = 500;
%% Go ahead and plot the initial observational error distribution
handles.figure1 = figure('Position', [figXmin figYmin figWidth figHeight], ...
    'Units', 'Pixels', ...
    'Name', 'proxy space update');

handles.h_obs_plot = plt_gaussian(handles.obs, handles.obs_error_sd, 1);
set(handles.h_obs_plot, 'Color', atts.red, 'Linewidth', 2.0);
hold on

% Plot an asterisk
handles.h_obs_ast = plot(handles.obs, 0, '*', 'MarkerSize', 10, 'Color', atts.red, 'LineWidth',2.0);
axis([xlower xupper ylower yupper]);
xlabel([ptype]);
ylabel('PDF');
set(gca, 'FontSize', atts.fontsize)
set(gca, 'XGrid', 'on')

% plot([xlower xupper], [0 0], 'k', 'Linewidth', 1.7);
% set(gca, 'FontUnits', 'Normalized');

handles.h_Hxb_climo = plot(handles.Hxb_climo, 1, '*', 'MarkerSize', 10, 'Color', atts.grey, 'LineWidth', 0.5);

handles.h_climo_pdf = plt_gaussian(mean(handles.Hxb_climo), std(handles.Hxb_climo), 1);
set(handles.h_climo_pdf, 'linewidth', 0.5, 'color', atts.grey);

handles.h_Hxb_ens = plot(prior_ens(1, :), 1, '*', 'MarkerSize', 10, 'Color', atts.green, 'LineWidth', 2.0);

handles.h_prior_pdf = plt_gaussian(mean(handles.Hxb_ens), std(handles.Hxb_ens), 1);
set(handles.h_prior_pdf, 'linewidth', 2, 'color', atts.green);

handles.h_Hxa_ens = plot(poste_ens(1, :), 0.5, '*', 'MarkerSize', 10, 'Color', atts.blue, 'LineWidth', 2.0);

handles.h_poste_pdf = plt_gaussian(mean(handles.Hxa_ens), std(handles.Hxa_ens), 1);
set(handles.h_poste_pdf, 'linewidth', 2, 'color', atts.blue);

% Plot lines connecting the prior and posterior ensemble members
for i = 1:plt_num
    x_line = [prior_ens(1, i), poste_ens(1, i)];
    y_line = [1, 0.5];
    handles.h_poste_lines(i) = plot(x_line, y_line, 'Color', atts.grey, 'LineWidth', 1.0);
end
hold off
Figure_name = ['Probability.png'];
print(handles.figure1,'-dpng',[Figure_dir Figure_name])


% Set a basic plotting domain range that includes mean +/- 3 obs SDs
obslower = min([handles.Hxb_ens handles.Hxa_ens]);
obsupper = max([handles.Hxb_ens handles.Hxa_ens]);
unobslower = min([handles.xb_ens(location,:) handles.xa_ens(location,:)]);
unobsupper = max([handles.xb_ens(location,:) handles.xa_ens(location,:)]);
xlower = 0;
xupper = 1;
ylower = 0;
yupper = 1;

atts = stylesheet;  % get the default fonts and color
figureXmin   = 100;
figureYmin   = 250;
figureWidth  = 550;
figureHeight = 500;
%% Create a figure layout
handles.figure1 = figure('Name', 'twod_ensemble', ...
    'Position', [figureXmin figureYmin figureWidth figureHeight]);
plotposition = [ 50/figureWidth 120/figureHeight  60/figureWidth 350/figureHeight; ...
    115/figureWidth 120/figureHeight 400/figureWidth 350/figureHeight; ...
    115/figureWidth  50/figureHeight 400/figureWidth  55/figureHeight];

%% Create the axis for the unobserved marginal graphic
handles.h_unobMarginal = axes('Position', plotposition(1,:), ...
    'FontName', atts.fontname, ...
    'FontSize', atts.fontsize);
axes(handles.h_unobMarginal);
handles.h_xb_ens = ...
    plot(0, prior_ens(2, :), '*', 'MarkerSize', 10, 'Color', atts.green, 'LineWidth',2.0);
axis([xlower xupper unobslower unobsupper]);
set(gca, 'FontUnits', 'Normalized');
ylabel('SAT', 'Fontsize', atts.fontsize, 'FontUnits', 'Normalized');
hold on

% Now need to sort the state variable ensemble to get nice ordering
[~, sort_ind] = sort(prior_ens(2, :));
for i = 1:plt_num
    y(i) = plt_point(i) / (handles.ens_size + 1);
    handles.h_marg_state(i) = ...
        plot(y(i), poste_ens(2, sort_ind(i)), '*', 'MarkerSize', 10, 'Color', atts.blue, 'LineWidth', 2.0);
    % Also plot a segment in grey
    handles.h_state_inc(i) = plot([y(i), y(i)], ...
        [prior_ens(2, sort_ind(i)), poste_ens(2, sort_ind(i))], 'Color', atts.grey, 'LineWidth', 1.0);
end

%% Create a subplot for the observed variable marginal
handles.h_obMarginal = axes('Position', plotposition(3,:), ...
    'FontName', atts.fontname, ...
    'FontSize', atts.fontsize);
axes(handles.h_obMarginal);
handles.h_Hxb_ens = ...
    plot(prior_ens(1, :), 0, '*', 'MarkerSize', 10, 'Color', atts.green, 'LineWidth', 2.0);
axis([obslower obsupper ylower yupper]);
xlabel(ptype, 'FontSize', atts.fontsize, 'FontUnits', 'Normalized');
set(gca, 'FontUnits', 'Normalized');
hold on
handles.h_obs_marg = plot(handles.obs, 0, '*', 'MarkerSize', 12, 'LineWidth', 2.0);
set(handles.h_obs_marg,'Color',atts.red)

% Need to sort ensemble to get nice ordering for increments
[~, sort_obs_ind] = sort(prior_ens(1, :));
for i = 1:plt_num
    handles.h_marg_update(i) = ...
        plot(poste_ens(1, sort_obs_ind(i)), y(i), '*', 'MarkerSize', 10, 'Color', atts.blue, 'LineWidth', 2.0);
    % Also plot a segment in grey
    handles.h_marg_inc(i) = ...
        plot([prior_ens(1, sort_obs_ind(i)), poste_ens(1, sort_obs_ind(i))], ...
        [y(i), y(i)], 'Color', atts.grey, 'LineWidth', 1.0);
end

%% Get a subplot for the joint
handles.h_joint = axes('Position', plotposition(2,:), ...
    'FontName', atts.fontname, ...
    'FontSize', atts.fontsize);
axes(handles.h_joint);
handles.h_Hxb_xb_ens = ...
    plot(prior_ens(1, :), prior_ens(2, :), '*', 'MarkerSize', 10, 'Color', atts.green, 'LineWidth', 2.0);
axis([obslower obsupper unobslower unobsupper]);
grid on
title('Joint Distribution', 'Fontsize', atts.fontsize, 'FontUnits', 'Normalized');
set(gca, 'FontUnits', 'Normalized')
hold on
handles.h_best_fit = plot(best_x, best_y, 'LineWidth', 2.0);
set(handles.h_best_fit, 'Color', atts.green);

for i = 1:plt_num
    handles.h_joint_update(i) = plot(poste_ens(1, i), poste_ens(2, i), ...
        '*', 'MarkerSize', 10, 'Color', atts.blue, 'LineWidth', 2.0);
    handles.h_joint_inc(i) = plot([prior_ens(1, i), poste_ens(1, i)], ...
        [prior_ens(2, i), poste_ens(2, i)], 'Color', atts.grey, 'LineWidth', 1.0);
end

% Turn off the unwanted tick labels
set(handles.h_unobMarginal, 'Xticklabel', []);
set(handles.h_joint       , 'Xticklabel', [], 'Yticklabel', []);
set(handles.h_obMarginal  , 'Yticklabel', []);
hold off
Figure_name = ['linearregression.png'];
print(handles.figure1,'-dpng',[Figure_dir Figure_name])