%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based off code from clyinder.m%                                                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Setup: STEP 0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GENERAL FLOW CONSTANTS
lx     = 400;      % number of cells in x-direction
ly     = 100;      % number of cells in y-direction
uMax   = 0.1;      % velocity ofinflow
Re     = 500;      % Reynolds number
nu     = uMax * ly / Re;  % kinematic viscosity
omega  = 1. / (3*nu+1./2.);      % relaxation parameter
maxT   = 40000;  % total number of iterations
tPlot  = 50;      % cycles between plots


% D2Q9 LATTICE CONSTANTS
w  = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
cx = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
cy = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];
%Note that the directions are shifted over one
opp = [ 1,   4,  5,  2,  3,    8,   9,   6,   7]; 
col = [2:(ly-1)]; %outside walls
in  = 1;   % position of inlet
out = lx;  % position of outlet

[y,x] = meshgrid(1:ly,1:lx); % get coordinate of matrix indices
  
obst_x = lx/5+1;   % position of the cylinder
obst_y = ly/2+3;   
obst_r = ly/10+1;  % radius of the cylinder
obst =(x-obst_x).^2 + (y-obst_y).^2 <= 1.5*obst_r.^2;
obst(:,[1,ly]) = 1 ;   % Location of top/bottom boundary
obst(290:310,40:60) = 1 ; % back square
obst(140:150,25:75) = 1 ; % rectangle 
bbRegion = find(obst); % Boolean mask for bounce-back cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTIAL CONDICTIONS: STEP 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIAL CONDITION: 
L = ly-2; y_phys = y-1.5;
rho0 = 1;
F=zeros(9,400,100);

for i=1:9
    F(i,:,:) = rho0 .* w(i);
end

% MAIN LOOP (TIME CYCLES)
for cycle = 1:maxT

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MACROSCOPIC VALUES: STEP 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rho = sum(F,1);
    ux = F(2,:,:) + F(6,:,:) + F(9,:,:) - F(4,:,:) - F(7,:,:) - F(8,:,:);
    ux=ux./rho;
    
    uy = F(3,:,:) + F(6,:,:) + F(7,:,:) - F(5,:,:) - F(8,:,:) - F(9,:,:);
    uy=uy./rho;
	  
    % MACROSCOPIC BOUNDARY CONDITIONS
    % Inlet: linear constant flow profile
    y_phys = col-1.5;
    ux(:,in,col) = uMax; 
    uy(:,in,col) = 0; % no input y velocity
    rho(:,in,col) = 1; 
    
    %Outlet: Constant density and 0 velocity
    rho(:,out,col) = 1;
    ux(:,out,col) = 0;
    uy(:,out,col)  = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ZOU/HE Boundry Condictions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % BOUNDARY CONDITIONS: INLET (Zou/He BC)
    F(2,in,col) = F(4,in,col) + 2/3*rho(1,in,col).*ux(1,in,col); 
    F(6,in,col) = F(8,in,col) + 1/2*(F(5,in,col)-F(3,in,col)) ... 
                                    + 1/2*rho(1,in,col).*uy(:,in,col) ...
                                    + 1/6*rho(1,in,col).*ux(:,in,col); 
    F(9,in,col) = F(7,in,col) + 1/2*(F(3,in,col)-F(5,in,col)) ... 
                                    - 1/2*rho(:,in,col).*uy(:,in,col) ...
                                    + 1/6*rho(:,in,col).*ux(:,in,col); 

    % BOUNDARY CONDITIONS: OUTLET (Zou/He BC)
    F(4,out,col) = F(2,out,col) - 2/3*rho(:,out,col).*ux(:,out,col); 
    F(8,out,col) = F(6,out,col) + 1/2*(F(3,out,col)-F(5,out,col)) ... 
                                      - 1/2*rho(:,out,col).*uy(:,out,col) ...
                                      - 1/6*rho(:,out,col).*ux(:,out,col); 
    F(7,out,col) = F(9,out,col) + 1/2*(F(5,out,col)-F(3,out,col)) ... 
                                      + 1/2*rho(:,out,col).*uy(:,out,col) ...
                                      - 1/6*rho(:,out,col).*ux(:,out,col); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EQUILBRIUM DISTRIBUTION: STEP 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for i=1:9
      cu = 3*(cx(i)*ux+cy(i)*uy);
      term1 =(cx(i)*ux+cy(i)*uy);
      term2 = term1.*term1;
      term3 = ux.^2+uy.^2;
      fEq(i,:,:)  = rho .* w(i) .*( 1 + 3* term1 + 4.5*term2  - 3/2*term3);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMULATE COLLISION: STEP 4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    for i=1:9
      fprop(i,:,:) = F(i,:,:) - omega .* (F(i,:,:)-fEq(i,:,:));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STREAMING: STEP 5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    for i=1:9
       F(i,:,:) = circshift(fprop(i,:,:), [0,cx(i),cy(i)]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BOUNDRY CONDICTIONS OBSTACLE BOUNCE-BACK: STEP 6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    for i=1:9
         F(i,bbRegion) = F(opp(i),bbRegion);
    end
    
    
    % VISUALIZATION
    if (mod(cycle,tPlot)==1) % displace every 50th plot
        u = reshape(sqrt(ux.^2+uy.^2),lx,ly); %absolute velocity at each point
        u(bbRegion) = 0; % solid regions given zero velocity
        imagesc(u'); %plots image
        axis equal off; drawnow
    end
end