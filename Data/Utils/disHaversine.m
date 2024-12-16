function dis = disHaversine(lon1,lat1,lon2,lat2)
    R = 6367;
    if(lon1 > 180)
        lon1 = lon1-360;   
    end
    if(lon2 > 180)
        lon2 = lon2-360;   
    end
    dlon = lon2-lon1;
    dlat = lat2-lat1;
    h = sin(dlat*pi/180.0/2)^2+cos(lat1*pi/180)*cos(lat2*pi/180.0)*sin(dlon*pi/180.0/2)^2;
    dis = 2*R*asin(sqrt(h));
end