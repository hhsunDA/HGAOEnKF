function Loc = GC(Proxy_grids_all,update_grids_all,cut)
    for iobs = 1:size(Proxy_grids_all,1)     
        for iCut = 1:size(update_grids_all,1)
            dist = disHaversine(Proxy_grids_all(iobs,1),Proxy_grids_all(iobs,2),update_grids_all(iCut,1),update_grids_all(iCut,2));
            r = dist ./ (0.5*cut);
            if ( dist >= cut )
                V(iobs,iCut) = 0.0;
            elseif ( dist >= 0.5*cut & dist < cut )
                V(iobs,iCut) = (r^5 / 12.0 - r^4 / 2.0 + r^3 * 5.0 / 8.0 + r^2 * 5.0 / 3.0 - 5.0*r + 4.0 - 2.0 / (3.0 * r));
            else
                V(iobs,iCut) = (r^5 * (-0.25 ) + r^4 / 2.0 + r^3 * 5.0/8.0 - r^2 * 5.0/3.0 + 1.0);
            end
        end
    end
    Loc = V;
end