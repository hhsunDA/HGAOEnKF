function global_mean = global_mean_nan(field,wt)
    simul_year = size(field,2);
    for iyear = 1:simul_year
        mask = ~isnan(field(:,iyear));
        global_mean(iyear,1) = sum((wt(mask).*field(mask,iyear)),1)./sum(wt(mask));
    end