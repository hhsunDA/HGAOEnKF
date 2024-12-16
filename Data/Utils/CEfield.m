function CE = CEfield(recon,verify)
    model_size = size(recon,1);
    for imodel = 1:model_size
        rec_time = squeeze(recon(imodel,:));
        ver_time = squeeze(verify(imodel,:));
        mask = ~(isnan(ver_time)|isnan(rec_time));
        if(sum(mask) >= 50) % set min length of validation years
            diff = rec_time-ver_time;
            if(var(ver_time(mask)) > 1e-3) % set min variance of SST
                CE(imodel,1) = 1-sum(diff(mask).^2)/(sum(mask)*var(ver_time(mask)));
            else
                CE(imodel,1) = nan;
            end
        else
            CE(imodel,1) = nan;
        end
    end