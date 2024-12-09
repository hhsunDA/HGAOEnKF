function RMSE = RMSEfield(recon,verify)
    model_size = size(recon,1);
    for imodel = 1:model_size
        rec_time = squeeze(recon(imodel,:));
        ver_time = squeeze(verify(imodel,:));
        mask = ~(isnan(ver_time)|isnan(rec_time));
        if(sum(mask) >= 50)
            diff = rec_time-ver_time;
            if(var(ver_time,'omitnan')>1e-4)
                RMSE(imodel,1) = sqrt(sum(diff(mask).^2)/sum(mask));
            else
                RMSE(imodel,1) = nan;
            end            
        else
                RMSE(imodel,1) = nan;
        end
    end