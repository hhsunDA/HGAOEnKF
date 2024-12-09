function xa = EnSRF(xb,obs_error_var,obs_total,Hx_size,proxy_idx,CMat,localize,serial_update)
    ensemble_size = size(xb,2);
    nobsgrid = length(obs_total);
    rn     = 1.0 / (ensemble_size - 1);
    if(serial_update)
        % Serial update
        xb_mean = mean(xb,2);
        xb_prime = xb - xb_mean;
        % initialize xa_mean, xa_prime
        xa_mean = xb_mean;
        xa_prime = xb_prime;
        for iobs = 1:nobsgrid
            xb_mean        = xa_mean;                              % (n+p) x 1
            xb_prime       = xa_prime;                            % (n+p) x ens
            Hxb_mean = xb_mean(end-Hx_size+1:end);
            Hxb_mean_sub = Hxb_mean(proxy_idx);

            Hx_prime   = xb_prime(end-Hx_size+1:end,:);
            Hx_prime_sub = Hx_prime(proxy_idx,:);
            hxb_mean       = Hxb_mean_sub(iobs);                             % 1 x 1
            hxb_prime      = Hx_prime_sub(iobs,:);                          % 1 x ens
            hpbht        = hxb_prime*hxb_prime'*rn;                     % 1 x 1
            gainfact     = (hpbht+obs_error_var(iobs))/hpbht*(1.0-sqrt(obs_error_var(iobs)/(hpbht+obs_error_var(iobs))));
            pbht         = (xb_prime*hxb_prime')*rn;                    % (n+p) x 1
            if(localize)
                Cvect    = CMat(iobs,:);
                kfgain   = Cvect'.*pbht/(hpbht+obs_error_var(iobs));% (n+p) x 1
            else
                kfgain   = pbht/(hpbht+obs_error_var(iobs));        % (n+p) x 1
            end
            xa_mean        = xb_mean + kfgain*(obs_total(iobs)-hxb_mean);     % (n+p) x 1 
            xa_prime       = xb_prime - gainfact*kfgain*hxb_prime;    % (n+p) x ens
            xa = xa_mean + xa_prime;
        end
    else
        % Matrix update
        xb_mean = mean(xb,2);
        xb_prime = xb - xb_mean;
        Hxb_mean = xb_mean(end-Hx_size+1:end);
        Hxb_mean_sub = Hxb_mean(proxy_idx); 
        R = diag(obs_error_var);

        Hx_prime = xb_prime(end-Hx_size+1:end,:);
        Hx_prime_sub = Hx_prime(proxy_idx,:);
        BHT = xb_prime*Hx_prime_sub'*rn;
        HBHT = Hx_prime_sub*Hx_prime_sub'*rn;
        KFgain = BHT/(HBHT+R);
        reduceKFgain = BHT/sqrtm(HBHT+R)/(sqrtm(HBHT+R)+sqrtm(R));    
        xa_mean = xb_mean + KFgain*(obs_total-Hxb_mean_sub);     
        xa_prime = xb_prime -reduceKFgain*Hx_prime_sub;
        xa = xa_mean + xa_prime;             
    end

end