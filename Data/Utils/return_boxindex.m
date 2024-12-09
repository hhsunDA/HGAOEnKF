function [lat1,lat2,lon1,lon2] = return_boxindex(j,lat_sparse,lon_sparse,obs_grids_lat,obs_grids_lon)

% lon is a circle
if(obs_grids_lon(j)<0)
    if(obs_grids_lon(j)>lon_sparse(end)-360 && obs_grids_lon(j)<=lon_sparse(end)-360+(lon_sparse(1)-lon_sparse(2))/2)
        lon1=length(lon_sparse);
        lon2=1;
    elseif(obs_grids_lon(j)>lon_sparse(end)-360+(lon_sparse(1)-lon_sparse(2))/2 && obs_grids_lon(j)<=lon_sparse(end)-360+(lon_sparse(1)-lon_sparse(2)))
        lon1=1;
        lon2=length(lon_sparse); 
    else
    [~,index]=sort((abs(lon_sparse-(obs_grids_lon(j)+360))));
        lon1=index(1);
        lon2=index(2);
    end
elseif(obs_grids_lon(j)>=0)
    if(obs_grids_lon(j)>lon_sparse(end)-360 && obs_grids_lon(j)<=lon_sparse(end)-360+(lon_sparse(1)-lon_sparse(2))/2)
        lon1=length(lon_sparse);
        lon2=1;
    elseif(obs_grids_lon(j)>lon_sparse(end)-360+(lon_sparse(1)-lon_sparse(2))/2 && obs_grids_lon(j)<=lon_sparse(end)-360+(lon_sparse(1)-lon_sparse(2)))
        lon1=1;
        lon2=length(lon_sparse);         
    else
    [~,index]=sort((abs(lon_sparse-obs_grids_lon(j))));
        lon1=index(1);
        lon2=index(2);
    end
end

[~,index]=sort(abs(lat_sparse-obs_grids_lat(j)));
lat1=index(1);
lat2=index(2);