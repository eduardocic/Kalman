function H = matrizH(xt, yt, xr, yr)
    
    r = sqrt((xt - xr)^2 + (yt - yr)^2);
    
    H(1,1) = -(yt - yr)/(r*r);
    H(1,2) = 0;
    H(1,3) = (xt - xr)/(r*r);
    H(1,4) = 0;
    
    H(2,1) = (xt - xr)/r;
    H(2,2) = 0;
    H(2,3) = (yt - yr)/r;
    H(2,4) = 0;
end
