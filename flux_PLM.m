function [flux] = flux_PLM(dir,dt,ds,vel,imu)
%a duplication of the function in ATHENA FullRT_flux.c

%the upwind slope
    delq1 = (imu(3) - imu(2))/ds;
    delq2 = (imu(2) - imu(1))/ds;
    if(delq1*delq2 >0.0)
        dqi0 = 2.0*delq1*delq2/(delq1+delq2);
    else
        dqi0 = 0.0;
    end
    
    distance = ds; 
    distance = distance - ((vel*dt)); 
    
    %The upwind flux
    flux = imu(2) - distance*(dqi0/2.0); 
end