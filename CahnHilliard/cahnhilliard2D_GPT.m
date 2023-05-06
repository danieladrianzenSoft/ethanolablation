clear
clc

%% initialize params

K = 3e-10;
r0 = 0.03; % Netti 2003, 0.035cm
%r0 = round(0.311/(2*10),3); %nominal inner diameter of 24 gauge needle is 0.2032mm. Here we change to cm. 
Rtot = 1; % Netti 2003, cm
%Q = (0.1/1e3)/60; % 0.1uL/min = 0.1mm^3/min Netti 2003
Q = 1/3600;
Vol = 0.50; %total injection volume -> CHANGE TO 0.5ML - 2.5ML for humans
lambda = 13.16*1e4; % Netti 2003, 13.16kPa = 13.16*1e4 Barye
%mu = [6.58*1e4,1*1e5]; % Netti 2003, 6.58kPa = 6.58*1e4 Barye
mu = 6.58*1e4;
phi = 0.2; % Netti 2003, Porosity, dimensionless
numr = 400;
t_end = 2*60*60;
t_injection_end = Vol/Q;
t_injection = linspace(0.01,t_injection_end,600);
t_relaxation = [t_injection_end+1,t_end];

%t = 0.01:2:20000;
%t = linspace(0.01,12*60*60,400);
%rho_etoh = 0.789; % density of ethanol.
c0 = 0.5;
D_S = 5*10^(-6); %(cm^2/s)
D_C = D_S;
P = 1;
%gamma = (0.01/((2*pi)^2));
%gamma = 1;
gamma = 0.01;

rc = 0.3;
rs = 0.2;

% values
Nx=100;
Ny=100;
x = linspace(-1,1,Nx);
y = linspace(-1,1,Ny);
dx = x(2)-x(1);
dy = y(2)-y(1);
dr = sqrt(dx^2+dy^2);
u0=zeros(Nx,Ny);
i0=0;
j0=0;


%% CAHN HILLIARD


[xx,yy]=meshgrid(x,y);
c_mask = (xx-i0).^2+(yy-j0).^2<=rc^2;

rc_perimeter = bwperim(u0);
u0(~rc_perimeter)=-1;
u0(c_mask)=1;
u0(rc_perimeter)=0;

[perim_x,perim_y] = find(rc_perimeter>0);
perim_inds = [perim_x,perim_y];

perim = sub2ind([Ny,Nx],perim_y,perim_x);

u0(~rc_perimeter)=0;
u0(rc_perimeter)=1;
u0(iPerim_inds)=1;

D_CH = 1e-6;

IC = reshape(u0,[],1);


% Variables for incremental capture
count = 1;
frameStep = 1;
tSpan = [0,options.EndTime];
S = getSparsity(length(x),length(y));
S = sparse(S);
%opts1 = odeset('Vectorized','on');
opts1 = odeset('JPattern',S,'Vectorized','on','RelTol',1e-3,'AbsTol',1e-4);
domain_x = 2:length(x)-1;
domain_y = 2:length(y)-1;
domain_x_rep = repmat(domain_x,1,length(domain_y));
domain_y_rep = repelem(domain_y,length(domain_x));
%domain_y_rep = repmat(domain_y,1,length(domain_x));
domain_inds = sub2ind([length(y),length(x)],domain_y_rep,domain_x_rep);
[perim_x,perim_y] = ind2sub([length(y),length(x)],perimeter);
%opts1 = odeset('Vectorized','on','RelTol',1e-3,'AbsTol',1e-4);

tic
[t,ut] = ode15s(@(t,u) cahn_hilliard(t,u,x,y,perimeter,domain_inds,D,gamma), tSpan, IC, opts1);
toc

u = zeros(length(x),length(y),length(t));
for i = 1:length(t)
    u(:,:,i) = reshape(ut(i,1:length(x)*length(y)),length(x),length(y));
end

writer = VideoWriter(options.FileName);
open(writer);
f = figure();
if options.ImgStyle == "binary"
    colormap("pink")
else
    colormap(options.Colormap)
end
ax = axes("Parent",f);
%u0(perimeter) = -1;
perimeter_mask = nan(length(y),length(x));
perimeter_mask(perimeter) = 1;
heatmap = pcolor(ax,x,y,u0);
set(heatmap,'EdgeColor','none')
colorbar
hold on
%maskPlot=pcolor(ax,x,y,perimeter_mask);
%set(maskPlot,'EdgeColor','none');

%set(maskPlot,'Color','w')
%circle=u0;
%circle(perimeter)=-1;
%pcolor(ax,x,y,circle);

%clim([-1,1])
clim([0,1])
frame = getframe(1);
writeVideo(writer,frame);

for k = 1:length(t)
    if (strcmp(options.CaptureMode,'incremental'))
       if (count == frameStep)
           if (strcmp(options.ImgStyle,'true'))
               %u(perim_x,perim_y,k) = -1;
               heatmap = pcolor(ax,x,y,u(:,:,k));
               %maskPlot = pcolor(ax,x,y,perimeter_mask);
               %set(maskPlot,'EdgeColor','none');
               set(heatmap,'EdgeColor','none')
               %circle=u(:,:,k);
               %circle(perimeter)=-1;
               %pcolor(ax,x,y,circle);
               colorbar
               %clim([-1,1])
               clim([0,1])
           else
               uMod = round((u+1)/2); % Binarizes image
               %uMod(perim_x,perim_y,k) = -1;
               heatmap = pcolor(ax,x,y,uMod(:,:,k));
               %maskPlot = pcolor(ax,x,y,perimeter_mask);
               %set(maskPlot,'EdgeColor','none');
               set(heatmap,'EdgeColor','none')
               %circle=u(:,:,k);
               %circle(perimeter)=-1;
               %pcolor(ax,x,y,circle);
               colorbar
               %clim([-1,1])
               clim([0,1])
           end
           frame = getframe(1);
           writeVideo(writer,frame);
           count = 0;
           frameStep = frameStep + 1;
       end
       count = count+1;
       % Standard video mode
       else
           if (mod(k,options.FrameSpacing) == 0)
               if (strcmp(options.ImgStyle,'true'))
                   heatmap = pcolor(ax,x,y,u(:,:,k));
                   set(heatmap,'EdgeColor','none')
                   colorbar
                   clim([0,1])
               else
                   uMod = round((u+1)/2);
                   heatmap = pcolor(ax,x,y,uMod(:,:,k));
                   set(heatmap,'EdgeColor','none')
                   colorbar
                   clim([0,1])
               end
               frame = getframe(1);
               writeVideo(writer,frame);
           end
    end
end

close(writer);


function out = cahn_hilliard(t,u,x,y,perimeter,domain_inds,D,gamma) 
% Inputs:
    % t - time point
    % u - concentration
    % x, y - 2D spatial vectors x and y in cartesian coordinates
    % perimeter - column indices of elements in x by y matrix that form an
    %   internal boundary, and where we will apply a no flux
    %   boundary condition
    % domain_inds - column indices of all elements in the x by y matrix
    %   excluding elements in the edges of the domain
    % D - diffusion coefficient (equivalent to mobility function)
    % gamma - surface tension parameter
% Outputs:
    % out
% Function:

        uLaplace = zeros(numel(x)*numel(y),size(u,2));
        uMu = zeros(numel(x)*numel(y),size(u,2));
        uMuLaplace = zeros(numel(x)*numel(y),size(u,2));
        dx = x(2)-x(1);
        dy = y(2)-y(1);
        Nx = length(x);
        Ny = length(y);

        u_m = u;

        %%%%%% BOUNDARY CONDITIONS %%%%%%%

        left = sub2ind([Ny,Nx],1:Ny,ones(1,Ny)); % left + Ny = i+1,j
        top = sub2ind([Ny,Nx],ones(1,Nx),1:Nx); % top + 1 = i,j+1
        right = sub2ind([Ny,Nx],1:Ny,Nx*ones(1,Ny)); % right - Ny = i-1,j
        bottom = sub2ind([Ny,Nx],Ny*ones(1,Nx),1:Nx); % bottom - 1 = i,j-1

        % top left corner
        uLaplace(left(1),:) = 2*(u_m(left(1)+Ny,:)-u_m(left(1),:))/dx^2 + 2*(u_m(left(1)+1,:)-u_m(left(1),:))/dy^2;
        uMu(left(1),:) = u_m(left(1),:).^3 - u_m(left(1),:) - gamma*uLaplace(left(1),:);
        uMuLaplace(left(1),:) = 2*(uMu(left(1)+Ny,:)-uMu(left(1),:))/dx^2 + 2*(uMu(left(1)+1,:)-uMu(left(1),:))/dy^2;
        % top right corner
        uLaplace(right(1),:) = 2*(u_m(right(1)-Ny,:)-u_m(right(1),:))/dx^2 + 2*(u_m(right(1)+1,:)-u_m(right(1),:))/dy^2;
        uMu(right(1),:) = u_m(right(1),:).^3 - u_m(right(1),:) - gamma*uLaplace(right(1),:);
        uMuLaplace(right(1),:) = 2*(uMu(right(1)-Ny,:)-uMu(right(1),:))/dx^2 + 2*(uMu(right(1)+1,:)-uMu(right(1),:))/dy^2;
        % bottom left corner
        uLaplace(left(end),:) = 2*(u_m(left(end)+Ny,:)-u_m(left(end),:))/dx^2 + 2*(u_m(left(end)-1,:)-u_m(left(end),:))/dy^2;
        uMu(left(end),:) = u_m(left(end),:).^3 - u_m(left(end),:) - gamma*uLaplace(left(end),:);
        uMuLaplace(left(end),:) = 2*(uMu(left(end)+Ny,:)-uMu(left(end),:))/dx^2 + 2*(uMu(left(end)-1,:)-uMu(left(end),:))/dy^2;
        % bottom right corner
        uLaplace(right(end),:) = 2*(u_m(right(end)-Ny,:)-u_m(right(end),:))/dx^2 + 2*(u_m(right(end)-1,:)-u_m(right(end),:))/dy^2;
        uMu(right(end),:) = u_m(right(end)).^3 - u_m(right(end),:) - gamma*uLaplace(right(end),:);
        uMuLaplace(right(end),:) = 2*(uMu(right(end)-Ny,:)-uMu(right(1),:))/dx^2 + 2*(uMu(right(end)-1,:)-uMu(right(end),:))/dy^2;
        % top
        uLaplace(top(2:end-1),:) = (u_m(top(2:end-1)-Ny,:)-2*u_m(top(2:end-1),:)+u_m(top(2:end-1)+Ny,:))/dx^2+2*(u_m(top(2:end-1)+1,:)-u_m(top(2:end-1),:))/dy^2;
        uMu(top(2:end-1),:) = u_m(top(2:end-1),:).^3 - u_m(top(2:end-1),:) - gamma*uLaplace(top(2:end-1),:);
        uMuLaplace(top(2:end-1),:) = (uMu(top(2:end-1)-Ny,:)-2*uMu(top(2:end-1),:)+uMu(top(2:end-1)+Ny,:))/dx^2+2*(uMu(top(2:end-1)+1,:)-uMu(top(2:end-1),:))/dy^2;
        % left
        uLaplace(left(2:end-1),:) = 2*(u_m(left(2:end-1)+Ny,:)-u_m(left(2:end-1),:))/dx^2+(u_m(left(2:end-1)-1,:)-2*u_m(left(2:end-1),:)+u_m(left(2:end-1)+1,:))/dy^2;
        uMu(left(2:end-1),:) = u_m(left(2:end-1),:).^3 - u_m(left(2:end-1),:) - gamma*uLaplace(left(2:end-1),:);
        uMuLaplace(left(2:end-1),:) = 2*(uMu(left(2:end-1)+Ny,:)-uMu(left(2:end-1),:))/dx^2+(uMu(left(2:end-1)-1,:)-2*uMu(left(2:end-1),:)+uMu(left(2:end-1)+1,:))/dy^2;
        % right
        uLaplace(right(2:end-1),:) = 2*(u_m(right(2:end-1)-Ny,:)-u_m(right(2:end-1),:))/dx^2 + (u_m(right(2:end-1)-1,:)-2*u_m(right(2:end-1),:)+u_m(right(2:end-1)+1,:))/dy^2;
        uMu(right(2:end-1),:) = u_m(right(2:end-1),:).^3 - u_m(right(2:end-1),:) - gamma*uLaplace(right(2:end-1),:);
        uMuLaplace(right(2:end-1),:) = 2*(uMu(right(2:end-1)-Ny,:)-uMu(right(2:end-1),:))/dx^2 + (uMu(right(2:end-1)-1,:)-2*uMu(right(2:end-1),:)+uMu(right(2:end-1)+1,:))/dy^2;
        % bottom
        uLaplace(bottom(2:end-1),:) = (u_m(bottom(2:end-1)-Ny,:)-2*u_m(bottom(2:end-1),:)+u_m(bottom(2:end-1)+Ny,:))/dx^2 + 2*(u_m(bottom(2:end-1)-1,:)-u_m(bottom(2:end-1),:))/dy^2;
        uMu(bottom(2:end-1),:) = u_m(bottom(2:end-1),:).^3 - u_m(bottom(2:end-1),:) - gamma*uLaplace(bottom(2:end-1),:);
        uMuLaplace(bottom(2:end-1),:) = (uMu(bottom(2:end-1)-Ny,:)-2*uMu(bottom(2:end-1),:)+uMu(bottom(2:end-1)+Ny,:))/dx^2 + 2*(uMu(bottom(2:end-1)-1,:)-uMu(bottom(2:end-1),:))/dy^2;

        %%%%%% EQUATIONS %%%%%%%

        uLaplace(domain_inds,:) = (u_m(domain_inds-Ny,:)-2*u_m(domain_inds,:)+u_m(domain_inds+Ny,:))/dx^2 + (u_m(domain_inds-1,:)-2*u_m(domain_inds,:)+u_m(domain_inds+1,:))/dy^2;
        uMu(domain_inds,:) = u_m(domain_inds,:).^3 - u_m(domain_inds,:) - gamma*uLaplace(domain_inds,:);
        uMuLaplace(domain_inds,:) = (uMu(domain_inds-Ny,:)-2*uMu(domain_inds,:)+uMu(domain_inds+Ny,:))/dx^2 + (uMu(domain_inds-1,:)-2*uMu(domain_inds,:)+uMu(domain_inds+1,:))/dy^2;

        %%%%%% BOUNDARY CONDITION AT INTERIOR CIRCULAR BOUNDARY %%%%%%

        [yPos_ind,xPos_ind] = ind2sub([Nx,Ny],perimeter);
        yPos_ind = Ny - yPos_ind + 1;
        xPos = x(xPos_ind)'./(sqrt(x(xPos_ind).^2 + y(yPos_ind).^2))';
        yPos = y(yPos_ind)'./(sqrt(x(xPos_ind).^2 + y(yPos_ind).^2))';

        u_m_ip1j = (u_m(perimeter-1,:).*yPos - u_m(perimeter+1,:).*yPos + u_m(perimeter-Ny,:).*xPos)./(xPos);
        u_m_in1j = (u_m(perimeter+Ny,:).*xPos + u_m(perimeter+1,:).*yPos - u_m(perimeter-1,:).*yPos)./(xPos);
        u_m_ijp1 = (u_m(perimeter-1,:).*yPos + u_m(perimeter-Ny,:).*xPos - u_m(perimeter+Ny,:).*xPos)./(yPos);
        u_m_ijn1 = (u_m(perimeter+Ny,:).*xPos - u_m(perimeter-Ny,:).*xPos + u_m(perimeter+1,:).*yPos)./(yPos);

        uLaplace(perimeter,:) = (u_m_in1j-2*u_m(perimeter,:)+u_m_ip1j)/dx^2 + (u_m_ijn1-2*u_m(perimeter,:)+u_m_ijp1)/dy^2;
        uMu(perimeter,:) = u_m(perimeter,:).^3 - u_m(perimeter,:) - gamma*uLaplace(perimeter,:);
        uMuLaplace(perimeter,:) = (uMu(perimeter-Ny,:)-2*uMu(perimeter,:)+uMu(perimeter+Ny,:))/dx^2 + (uMu(perimeter-1,:)-2*uMu(perimeter,:)+uMu(perimeter+1,:))/dy^2;
        
        duT = D*uMuLaplace;

        %%%%% SETTING OUTPUT %%%%%%%

        out = duT;

    end 