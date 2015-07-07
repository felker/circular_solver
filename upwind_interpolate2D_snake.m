function i_flux = upwind_interpolate2D_snake(intensity, dt, ds, v,is,ie,js,je,dir)
%UPWIND_INTERPOLATE1D_decoupled For finite volume, pick one of several methods for
%advection step into cell. Does both directions at once
% Does only one ray at a time
    
%First index is LHS boundary
%Second index is top boundary 
[nx,ny] = size(intensity);
i_flux = zeros(nx,ny);
Jarray = zeros(5,1); 

if dir==1
    for j=js:je
        for i=is:ie
            vel = 0.5*(v(i-1,j) + v(i,j)); 
            if(vel > 0.0)
                for k=1:3
                    %NOTE!! ATHENA starts at is=3
                    Jarray(k) = intensity(i-3+k,j); 
                end
                i_flux(i,j) = flux_PLM(dir,dt,ds,vel,Jarray);
            elseif(vel <0.0)
                for k=1:3
                    Jarray(k) = intensity(i+2-k,j);
                end
                i_flux(i,j) = flux_PLM(dir,dt,ds,-vel,Jarray); 
            end
        end
    end
elseif dir==2
    for i=is:ie
        for j=js:je
            vel = 0.5*(v(i,j-1) + v(i,j)); 
            if(vel > 0.0)
                for k=1:3
                    %NOTE!! ATHENA starts at is=3
                    Jarray(k) = intensity(i,j-3+k); 
                end
                i_flux(i,j) = flux_PLM(dir,dt,ds,vel,Jarray);
            elseif(vel <0.0)
                for k=1:3
                    Jarray(k) = intensity(i,j+2-k);
                end
                i_flux(i,j) = flux_PLM(dir,dt,ds,-vel,Jarray); 
            end
        end
    end
end
    
    


