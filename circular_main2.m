%circular_main2.m
%Sovles 2D vacuum RTE in cylindrical coordinates

%this only tests advection step, no source terms allowed right now

%------------------------ PARAMETERS ------------------ %
clear all;
N = 12; 
lx1 = 1.0; 
lx2 = 2*pi; 
nx1 = 128;
nx2 = 128;

c = 1.0;

%Using properly cell centered positions makes automatically slighlty more
%challenging. 

%In both dimensions, we don't go all the way to one end of the interval

%specify final radius lr and number of active radial cells, nr. code fixes
%first point at dr/2= lr/2*nr
%then, the actual length of the domain is lr+dr
dx1 = lx1/(nx1 - 1/2);
dx2 = lx2/(nx2);

dt = 0.00001;
nt = 20000;
normalization_tol = 1e-6;

%------------------------ BOUNDARY CONDITIONS  ------------------ %
%%% van Leer interpolation requires two bc 
num_ghost = 2; 
is = 1+num_ghost; ie = nx1+num_ghost;
js= 1+num_ghost; je= nx2+num_ghost;

%number of physical cells
nx1_r = nx1;
nx2_r = nx2;
nx1 = nx1+2*num_ghost;
nx2 = nx2+2*num_ghost;

%------------------------ BUILD SPATIAL MESH ------------------ %
%Radiation Spatial Discretization
%Radiation samples are volume averaged quantites
%these should be cell centers
rr=linspace(dx1/2,lx1,nx1_r)';
tt=linspace(0,lx2-dx2,nx2_r)'; %should modify to make mod(2pi) arithmetic automatic

%make new array with one extra point to handle boundaries 
rr_b=linspace(0,lx1+dx1/2,nx1_r+1)';
%dont repeat boundaries in circular array!
tt_b=linspace(dx2/2,lx2-dx2/2,nx2_r)';

%plot them
% hold on;
% for i=1:nr_r+1
%     %the r boundaries are straight lines, not arcs here
%     %this should actually use a sampled parameterized theta to get arcs
%     polar([tt_b' tt_b(1)]',rr_b(i).*ones(ntheta_r+1,1),'-k'); 
%     %concatenate last boundary line
% end
% for i=1:ntheta_r
%     polar(tt_b(i).*ones(nr_r+1,1),rr_b,'-k');
% end
% hold off;

%------------------------ BUILD MOMENTUM ANGULAR MESH ------------------ %
%Basic uniform Cartesian basis angular discretization
%This produces a local discretization when the metric is in the canonical
%form. I.e. the angular parameterization refers to a local
%orthonormal tetrad
[ncells,nxa,mu,mu_b,pw] = uniform_angles2D(N,pi/2);
phi_bin = 1; %used to select phi level for plotting, injection of IC

%There is a notational confusion whereby nxa(2) includes the two poles.
%Therefore, there are actually nxa(2)-1 cells in the phi direction
%whereas in theta, we dont duplicate the last boundary ray (it is implied that 
% it is cyclic). Maybe want to be consistent in future
if (nxa(2) -1) >= 2
    num_phi_cells = nxa(2)-1; 
else
    num_phi_cells = 1;
end
dxA1 = 2*pi./(nxa(1)); 
dxA2 = pi./(nxa(2));
     
%Precompute the inverse tetrad transformation of the discretized rays
%Coordinate Basis components
mu_coordinate = zeros(nx1,nx2,nxa(1),num_phi_cells,3);
mu_b_coordinate = zeros(nx1,nx2,nxa(1),nxa(2),3);
for j=1:nx2
    %basis bectors only change with phi coordinate
    if j >= js && j<=je %deal with ghost cells 
        cos_tt = cos(tt(j-num_ghost));
        sin_tt = sin(tt(j-num_ghost));
    elseif j < js
        cos_tt = cos(tt(nx2_r-2+j));
        sin_tt = sin(tt(nx2_r-2+j));
    elseif j > je
        cos_tt = cos(tt(j-je));
        sin_tt = sin(tt(j-je));
    end
    for i=1:nx1
        if i >= is && i<=ie %deal with radial ghost cells 
            rr2 = rr(i-num_ghost); 
        elseif i < is
            rr2 = rr(is-num_ghost); 
        elseif i > ie
            rr2 = rr(ie-num_ghost);
        end
        for k=1:nxa(1)
            for l=1:num_phi_cells;   
                mu_coordinate(i,j,k,l,1) = mu(k,l,1)*cos_tt + sin_tt*mu(k,l,2); 
                mu_coordinate(i,j,k,l,2) = -mu(k,l,1)*sin_tt/rr2 + cos_tt*mu(k,l,2)/rr2;
                mu_b_coordinate(i,j,k,l,1) = mu_b(k,l,1)*cos_tt + sin_tt*mu_b(k,l,2); 
                mu_b_coordinate(i,j,k,l,2) = -mu_b(k,l,1)*sin_tt/rr2 + cos_tt*mu_b(k,l,2)/rr2;
                %add \hat{k}^3?
                assert(abs(polar_norm(squeeze(mu_b_coordinate(i,j,k,l,:)),rr2)- 1.0) < normalization_tol);
                assert(abs(polar_norm(squeeze(mu_coordinate(i,j,k,l,:)),rr2)- 1.0) < normalization_tol);
            end
            l=nxa(2);
            mu_b_coordinate(i,j,k,l,1) = mu_b(k,l,1)*cos_tt + sin_tt*mu_b(k,l,2); 
            mu_b_coordinate(i,j,k,l,2) = -mu_b(k,l,1)*sin_tt/rr2 + cos_tt*mu_b(k,l,2)/rr2;
            assert(abs(polar_norm(squeeze(mu_b_coordinate(i,j,k,l,:)),rr2)- 1.0) < normalization_tol);
        end
    end
end
%--------------------- CHARACTERIZE GEODESIC CURVATURE ------------------ %
% C^i and dx^A/dk can be precomputed since they do not depend on coord time
% in this coordinate system

% Compute C^i dtheta/dk^i
 vtheta = zeros(nx1_r,nx2_r,nxa(1),num_phi_cells);
 coordinate_source = zeros(nx1_r,nx2_r,nxa(1),num_phi_cells);
 for n=1:nx1_r  
     for p=1:nx2_r
         for j=1:nxa(1)
             for l=1:num_phi_cells
                 %atan2 returns negative angles. Only in those cases do i
                 %want to add 2pi. Fix this in snake_angles!!
                 xAtemp = atan2(mu(j,l,2),mu(j,l,1));
                 xAtemp = xAtemp + 2*pi*(xAtemp <0); 
                 vtheta(n,p,j,l) = -(sin(xAtemp - tt(p))^3 +(sin(3*xAtemp ...
                     -tt(p)) + sin(xAtemp - 3*tt(p)))/2)/rr(n); 
                 coordinate_source(n,p,j,l) = -(mu(j,l,1)*cos(tt(p)) +mu(j,l,2)*sin(tt(p)) + ...
                     3*sin(xAtemp-tt(p))^2*cos(xAtemp-tt(p)) +(3*cos(3*xAtemp - ...
                     tt(p)) + cos(xAtemp - 3*tt(p)))/2)/rr(n); 
             end
         end
     end
 end
%------------------------ GR PHASE VOLUME ------------------ %
cell_vol = zeros(nx1_r,nx2_r,nxa(1),nxa(2));
for i=1:nxa(1)
    for j=1:nxa(2)
    cell_vol(:,:,i,j) = dx1*dx2*dxA1*(rr*ones(nx2_r,1)');
    end
end

%------------------------ INTENSITY AND FLUID VELOCITY ------------------ %
%Monochromatic specific intensity, boundary conditions at 2,nz-1 
intensity = zeros(nx1,nx2,nxa(1),num_phi_cells); 

%--------------- JIANG14 variables, dimensionless constants --------------%
a_r =  5.670373e-8; %stefan boltzman constant
adiabatic = 5/3; %gamma=adiabatic index
R = 1; %8.311; %ideal gas constant

%characteristic values (arbitrary?)
a_0 = 0.1; %characteristic velocity
T_0 = 1.0; %characteristic temperature
P_0 = a_r;%1.0; %characteristic gas pressure

C = c/a_0; %dimensionlesss speed of light
P = a_r*T_0^4/P_0;  %measure of relative importance of rad and gas pressure

%------------------------ PROBLEM SETUP ------------------ %
%Absorption opacities
rho_a = zeros(nx1,nx2);
%Scattering opacities
rho_s = zeros(nx1,nx2);
%Fluid density, temperature
density = ones(nx1,nx2);
temp = ones(nx1,nx2);

%------------------------ PRE-TIMESTEPPING SETUP ------------------ %
%Calculate Radiation CFL numbers
% CFL CONDITION IS NOT TRIVIAL-- need to do stability analysis

%------------------------ OUTPUT VARIABLES------------------------------ %
output_interval = 1000; 
num_output = 8; %number of data to output
num_pts = nt/output_interval; 
time_out = dt*linspace(0,nt+output_interval,num_pts+1); %extra pt for final step
y_out = zeros(num_pts+1,num_output);

%Explicit-Implicit operator splitting scheme
%-------------------------------------------
for i=0:nt
    time = dt*i; %#ok<NOPTS> 
    %Substep 0: Time series output, boundary condition, moments update
    %Update moments
    if ~mod(i,output_interval)
    end
    
    %Boundary conditions
    for j=1:num_ghost
        %absorbing bcs at max radius
        intensity(ie+j,:,:,:) = 0.0;
        %absorbing bcs at min radius 
        intensity(is-j,:,:,:) = 0.0; %shouldnt need to do this
        %periodic bcs at max and min theta
        intensity(:,js-j,:,:) = intensity(:,je-j+1,:,:);
        intensity(:,je+j,:,:) = intensity(:,js+j-1,:,:);
    end  
    
    %DEBUG
    %Make the innermost radial ring be absorbing to avoid volume issues
    intensity(is,:,:,:) = 0.0;

    %Inject thick ray from top right edge of cylinder downward and to the
    %left
    injection_ir = ie;
    %first quadrant of theta
    injection_jtheta = num_ghost+nx2_r/8; % center of beam
    %width of beam, 1/4 of quadrant
    beam_width = nx2_r/16;
    %we want the injected beam to be close to the -\partial_r basis vector 
    injection_theta = 7; 
    theta_width = 4; 
    injection_phi = phi_bin;
    for j=1:num_ghost
        intensity(injection_ir+j,(injection_jtheta-beam_width/2+1):(injection_jtheta+1)...%beam_width/2)...
        ,injection_theta... %,(injection_theta-theta_width/2+1):(injection_theta+theta_width/2)...
    ,injection_phi) = 1.0;
    end 
    
    %Substep #1: Explicitly advance transport term
    net_flux = zeros(nx1,nx2,nxa(1),num_phi_cells);
    %r-flux
    for j=1:nxa(1) %do all nx, ny at once
        for l=1:num_phi_cells 
        %cannot pull mu out from partial, since it changes with x-position 
        i_flux = upwind_interpolate2D_snake(mu_coordinate(:,:,j,l,1).*(intensity(:,:,j,l)),dt,dx1,C.*sign(mu_coordinate(:,:,j,l,1)),is,ie+1,js,je+1,1);
        %net_flux(is:ie,js:je,j,l) = dt*C/(dx2*dx1)*(dx2*i_flux(is+1:ie+1,js:je) - dx2*[zeros(nx1_r, 1) ones(nx1_r, nx2_r -1)].*i_flux(is:ie,js:je));
        net_flux(is:ie,js:je,j,l) = dt*C*(dxA1*dx2*(rr_b(2:nx1_r+1,1)*ones(nx2_r,1)').*i_flux(is+1:ie+1,js:je) - (dxA1*dx2*rr_b(1:nx1_r,1)*ones(nx2_r,1)').*i_flux(is:ie,js:je))./(cell_vol(:,:,j,l));
        end
    end %end of ray trace loop, x-direction
    
%     %phi-flux
    for j=1:nxa(1) %do all nx, ny at once
        for l=1:num_phi_cells
            %THIS IS WHERE MONOTONICITY IS BROKEN. MU_COORDINATE(:,2) is
            %huge, unlike mu_coordinate(:,1)
        i_flux = upwind_interpolate2D_snake(mu_coordinate(:,:,j,l,2).*(intensity(:,:,j,l)),dt,dx2,C.*sign(mu_coordinate(:,:,j,l,2)),is,ie+1,js,je+1,2);
     %   net_flux(is:ie,js:je,j,l) = net_flux(is:ie,js:je,j,l)+ dt*C/dx2*(i_flux(is:ie,js+1:je+1) - i_flux(is:ie,js:je));
        net_flux(is:ie,js:je,j,l) = net_flux(is:ie,js:je,j,l)+ dt*C*(dx1*dxA1*i_flux(is:ie,js+1:je+1) - dx1*dxA1*i_flux(is:ie,js:je))./(cell_vol(:,:,j,l));
        end
    end    %end of ray trace loop, y-direction
    
     %Substep #1.1: Compute solid angular fluxes        
     %Manually hardcode the donor cell method since we only have vtheta
     angular_flux = zeros(nx1_r,nx2_r,nxa(1),num_phi_cells);
     i_flux = zeros(nx1_r,nx2_r,nxa(1),num_phi_cells,2);

     for n=1:nx1_r
        for p=1:nx2_r
            for j=1:nxa(1) %edit out for easy periodic circular shifts
                for l=1:num_phi_cells
                    %velocity at interface u_{j+1/2} via linear
                    %interpolation
                    if j-1 ==0
                        ave_vel = 0.5*(vtheta(n,p,j,l) + vtheta(n,p,nxa(1),l));
                    else
                        ave_vel = 0.5*(vtheta(n,p,j,l) + vtheta(n,p,j-1,l));
                    end
                    if ave_vel > 0
                    %if vtheta(n,p,j,l) > 0 %counterclockwise speed
                        %debugging nonconservation of angular flux:
                        %-- maybe this shouldnt be periodic in \theta
                        if j-1 == 0
                          i_flux(n,p,j,l,1) = ave_vel*intensity(n+num_ghost,p+num_ghost,nxa(1),l,1); %flux entering "left" boundary                          
                        else
                          i_flux(n,p,j,l,1) = ave_vel*intensity(n+num_ghost,p+num_ghost,j-1,l,1); %flux entering "left" boundary                          
                        end
                    %elseif vtheta(n,p,j,l) < 0 %clockwise
                    elseif ave_vel < 0
                        i_flux(n,p,j,l,1) = ave_vel*intensity(n+num_ghost,p+num_ghost,j,l,1);  %flux leaving left boundary          
                    end
                end
            end %end of angular loops
            for j=1:nxa(1)
                for l=1:num_phi_cells
                    if j == nxa(1)                 
                       % angular_flux(n,p,j,l) = dt*C/dxA1.*(i_flux(n,p,1,l,1) - i_flux(n,p,j,l,1));                   
                        angular_flux(n,p,j,l) = dt*C./cell_vol(n,p,j,l).*(dx2*dx1*rr(n).*i_flux(n,p,1,l,1) - dx2*dx1*rr(n).*i_flux(n,p,j,l,1));                   
                    else
                       % angular_flux(n,p,j,l) = dt*C/dxA1.*(i_flux(n,p,j+1,l,1) - i_flux(n,p,j,l,1));
                        angular_flux(n,p,j,l) = dt*C./cell_vol(n,p,j,l).*(dx2*dx1*rr(n).*i_flux(n,p,j+1,l,1) - dx2*dx1*rr(n).*i_flux(n,p,j,l,1));
                    end

                end
            end
            
        end
    end %end of spatial loops
    %check that the angular flux formulation is truly conservative
%     if(abs(max(max(max(sum(angular_flux,3))))) > normalization_tol)
%        error('conservation in local angular bins violated!');
%     end
    %Substep #2: Add coordinate source terms
    intensity_new = intensity(is:ie,js:je,:,:) ...
                        -net_flux(is:ie,js:je,:,:) - angular_flux;
    intensity_new = intensity_new + dt*C*intensity_new.*(coordinate_source); 
%     if(min(min(min(intensity_new))) < 0.0)
%        [Y I] = min(intensity_new); 
%         I = squeeze(I);
%         Y = squeeze(Y);
%         [Y2 I2] = max(Y);
%         [Y3 I3] = max(Y2); %these are scalars
%         %climb back up the tree
%         indices = zeros(1,3);
%         indices(3) = I3;
%         indices(2) = I2(indices(3));
%         indices(1) = I(indices(2),indices(3));
%         indices
%         intensity_new(indices(1),indices(2),indices(3))
%                error('negative intensity!');
%     end    
%     
%     %DEBUG THIS  
%      if max(max(max(max(abs(intensity_new))))) > 1.0
%         %tiny function for finding a certain element in my 4D intensity array
%         [Y I] = max(intensity_new); 
%         I = squeeze(I);
%         Y = squeeze(Y);
%         [Y2 I2] = max(Y);
%         [Y3 I3] = max(Y2); %these are scalars
%         %climb back up the tree
%         indices = zeros(1,3);
%         indices(3) = I3;
%         indices(2) = I2(indices(3));
%         indices(1) = I(indices(2),indices(3));
%         indices
%         intensity_new(indices(1),indices(2),indices(3))
%         %truncate the cached data at current timestep
%         return;
%      end
     intensity(is:ie,js:je,:,:) = intensity_new;
    %------------------------ NON-TIME SERIES OUTPUT ------------------ %
    if ~mod(i,output_interval)
        time
        h1 = figure(2);
        clf(h1); 
        set(h1,'name','Covariant polar solution','numbertitle','off');
        time_title = sprintf('t = %.3f (s)',time); %dont know where to add this above all subplots
        for j=1:nxa(1)
            %Ray intensity plots
            l=phi_bin; %select phi bin
            h(j) = subplot(3,4,j);%,'Parent',figure(2),'Clim',[0 1]); 
            
            %have to augment the intensity data with zeros corresponding to
            % the zero radius position. pcolor help page unclear which
            % vertex corresponds to shaded cell
            %use padarray to add a row and column of zeros
            intensity_plot=padarray(intensity(is:ie,js:je,j,l),[1,1],'post');
            %since we pass matrices for the coordinates, do not transpose
            %intensity matrix
            h(j) = pcolor(rr_b*cos([6.1850 tt_b' ]),rr_b*sin([6.1850 tt_b']),intensity_plot);
            %turn off grid lines
            set(h(j), 'EdgeColor', 'none');
            axis equal tight
            caxis manual
            caxis([0 1]);
            %this string formatter doesnt work on laptop for some reason
            %subtitle = sprintf('$$\hat{k}^i_{Cartesian}$$ =(%0.3f, %0.3f,%0.3f)',mu(j,l,1),mu(j,l,2),mu(j,l,3));
            %subtitle = sprintf('mu =(%0.3f, %0.3f,%0.3f)',mu(j,l,1),mu(j,l,2),mu(j,l,3));
            %title(subtitle);
            subtitle = ['$$\hat{k}^i_{Cartesian} = $$ (',num2str(mu(j,l,1),'%.3f'),',',num2str(mu(j,l,2),'%.3f'),',',...
                num2str(mu(j,l,3),'%.3f'),')'];
            title(subtitle,'Interpreter','latex');      
            %xlabel('r');
            %ylabel('y + A sin(kx)');
            colorbar
       end
        pause(0.1)
    end
end

%Plot mean intensity
J = zeros(nx1_r,nx2_r);
for i=1:nx1_r
    for j=1:nx2_r
        for k=1:nxa(1)
            J(i,j) = J(i,j) + intensity(i+num_ghost,j+num_ghost,k,1)*pw(k);
        end
    end
end
J_plot=padarray(J,[1,1],'post');
figure
h=pcolor(rr_b*cos([6.1850 tt_b' ]),rr_b*sin([6.1850 tt_b']),J_plot);
set(h, 'EdgeColor', 'none');
