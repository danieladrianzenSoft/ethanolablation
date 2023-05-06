%2D cartesian coordinate spinodal decomposition
    function out = cahn_hilliard_2D(t,u,x,y,perim_inds,perimeter,domain_inds,D,gamma) 
        %out = zeros(length(x)*length(y),1);
        uLaplace = zeros(numel(x),numel(y));
        %uLaplaceTest = zeros(numel(x),numel(y));
        uMu = zeros(numel(x),numel(y));
        uMuLaplace = zeros(numel(x),numel(y));
        %uLaplace = zeros(numel(x)*numel(y),size(u,2));
        %uMu = zeros(numel(x)*numel(y),size(u,2));
        %uMuLaplace = zeros(numel(x)*numel(y),size(u,2));
        dx = x(2)-x(1);
        dy = y(2)-y(1);
        %dr = sqrt(dx^2+dy^2);
        Nx = length(x);
        Ny = length(y);

        %%%%%% RESHAPING INPUT TO 2X2 MATRIX %%%%%%%
        u_m = reshape(u(1:length(x)*length(y)),length(x),length(y));
        %u_m = u;

        %%%%%% BOUNDARY CONDITIONS %%%%%%%

        %BC: x=0, y=0, (x_n=1, y_n=1), dnU=0, J.n=0
        %bound = sub2ind(size(u),1,1);
        %boundip1 = sub
        left = sub2ind([Ny,Nx],1:Ny,ones(1,Ny)); % left + Ny = i+1,j
        top = sub2ind([Ny,Nx],ones(1,Nx),1:Nx); % top + 1 = i,j+1
        right = sub2ind([Ny,Nx],1:Ny,Nx*ones(1,Ny)); % right - Ny = i-1,j
        bottom = sub2ind([Ny,Nx],Ny*ones(1,Nx),1:Nx); % bottom - 1 = i,j-1
        %left = sub2ind([Ny,Nx],2:Ny-1,ones(1,Ny-2)); % left + Ny = i+1,j
        %top = sub2ind([Ny,Nx],ones(1,Nx-2),2:Nx-1); % top + 1 = i,j+1
        %right = sub2ind([Ny,Nx],2:1:Ny-1,Nx*ones(1,Ny-2)); % right - Ny = i-1,j
        %bottom = sub2ind([Ny,Nx],Ny*ones(1,Nx-2),2:Nx-1); % bottom - 1 = i,j-1

        % top left corner
        %uLaplaceTest(1,1) = 2*(u_m(2,1)-u_m(1,1))/dx^2 + 2*(u_m(1,2)-u_m(1,1))/dy^2;
        uLaplace(left(1)) = 2*(u_m(left(1)+Ny)-u_m(left(1)))/dx^2 + 2*(u_m(left(1)+1)-u_m(left(1)))/dy^2;
        uMu(left(1)) = u_m(left(1)).^3 - u_m(left(1)) - gamma*uLaplace(left(1));
        uMuLaplace(left(1)) = 2*(uMu(left(1)+Ny)-uMu(left(1)))/dx^2 + 2*(uMu(left(1)+1)-uMu(left(1)))/dy^2;
        % top right corner
        %uLaplaceTest(Nx,1) = 2*(u_m(Nx-1,1)-u_m(Nx,1))/dx^2 + 2*(u_m(Nx,2)-u_m(Nx,1))/dy^2;
        uLaplace(right(1)) = 2*(u_m(right(1)-Ny)-u_m(right(1)))/dx^2 + 2*(u_m(right(1)+1)-u_m(right(1)))/dy^2;
        uMu(right(1)) = u_m(right(1)).^3 - u_m(right(1)) - gamma*uLaplace(right(1));
        uMuLaplace(right(1)) = 2*(uMu(right(1)-Ny)-uMu(right(1)))/dx^2 + 2*(uMu(right(1)+1)-uMu(right(1)))/dy^2;
        % bottom left corner
        %uLaplaceTest(1,Ny) = 2*(u_m(2,Ny)-u_m(1,Ny))/dx^2 + 2*(u_m(1,Ny-1)-u_m(1,Ny))/dy^2;
        uLaplace(left(end)) = 2*(u_m(left(end)+Ny)-u_m(left(end)))/dx^2 + 2*(u_m(left(end)-1)-u_m(left(end)))/dy^2;
        uMu(left(end)) = u_m(left(end)).^3 - u_m(left(end)) - gamma*uLaplace(left(end));
        uMuLaplace(left(end)) = 2*(uMu(left(end)+Ny)-uMu(left(end)))/dx^2 + 2*(uMu(left(end)-1)-uMu(left(end)))/dy^2;
        % bottom right corner
        %uLaplaceTest(Nx,Ny) = 2*(u_m(Nx-1,Ny)-u_m(Nx,Ny))/dx^2 + 2*(u_m(Nx,Ny-1)-u_m(Nx,Ny))/dy^2;
        uLaplace(right(end)) = 2*(u_m(right(end)-Ny)-u_m(right(end)))/dx^2 + 2*(u_m(right(end)-1)-u_m(right(end)))/dy^2;
        uMu(right(end)) = u_m(right(end)).^3 - u_m(right(end)) - gamma*uLaplace(right(end));
        uMuLaplace(right(end)) = 2*(uMu(right(end)-Ny)-uMu(right(1)))/dx^2 + 2*(uMu(right(end)-1)-uMu(right(end)))/dy^2;
        % top
        %uLaplaceTest(2:Nx-1,1) = (u_m(1:Nx-2,1)-2*u_m(2:Nx-1,1)+u_m(3:Nx,1))/dx^2 + 2*(u_m(2:Nx-1,2)-u_m(2:Nx-1,1))/dy^2;
        uLaplace(top(2:end-1)) = (u_m(top(2:end-1)-Ny)-2*u_m(top(2:end-1))+u_m(top(2:end-1)+Ny))/dx^2+2*(u_m(top(2:end-1)+1)-u_m(top(2:end-1)))/dy^2;
        uMu(top(2:end-1)) = u_m(top(2:end-1)).^3 - u_m(top(2:end-1)) - gamma*uLaplace(top(2:end-1));
        uMuLaplace(top(2:end-1)) = (uMu(top(2:end-1)-Ny)-2*uMu(top(2:end-1))+uMu(top(2:end-1)+Ny))/dx^2+2*(uMu(top(2:end-1)+1)-uMu(top(2:end-1)))/dy^2;
        % left
        %uLaplaceTest(1,2:Ny-1) = 2*(u_m(2,2:Ny-1)-u_m(1,2:Ny-1))/dx^2 + (u_m(1,1:Ny-2)-2*u_m(1,2:Ny-1)+u_m(1,3:Ny))/dy^2;
        uLaplace(left(2:end-1)) = 2*(u_m(left(2:end-1)+Ny)-u_m(left(2:end-1)))/dx^2+(u_m(left(2:end-1)-1)-2*u_m(left(2:end-1))+u_m(left(2:end-1)+1))/dy^2;
        uMu(left(2:end-1)) = u_m(left(2:end-1)).^3 - u_m(left(2:end-1)) - gamma*uLaplace(left(2:end-1));
        uMuLaplace(left(2:end-1)) = 2*(uMu(left(2:end-1)+Ny)-uMu(left(2:end-1)))/dx^2+(uMu(left(2:end-1)-1)-2*uMu(left(2:end-1))+uMu(left(2:end-1)+1))/dy^2;
        % right
        %uLaplaceTest(Nx,2:Ny-1) = 2*(u_m(Nx-1,2:Ny-1)-u_m(Nx,2:Ny-1))/dx^2 + (u_m(Nx,1:Ny-2)-2*u_m(Nx,2:Ny-1)+u_m(Nx,3:Ny))/dy^2;
        uLaplace(right(2:end-1)) = 2*(u_m(right(2:end-1)-Ny)-u_m(right(2:end-1)))/dx^2 + (u_m(right(2:end-1)-1)-2*u_m(right(2:end-1))+u_m(right(2:end-1)+1))/dy^2;
        uMu(right(2:end-1)) = u_m(right(2:end-1)).^3 - u_m(right(2:end-1)) - gamma*uLaplace(right(2:end-1));
        uMuLaplace(right(2:end-1)) = 2*(uMu(right(2:end-1)-Ny)-uMu(right(2:end-1)))/dx^2 + (uMu(right(2:end-1)-1)-2*uMu(right(2:end-1))+uMu(right(2:end-1)+1))/dy^2;
        % bottom
        %uLaplaceTest(2:Nx-1,Ny) = (u_m(1:Nx-2,Ny)-2*u_m(2:Nx-1,Ny)+u_m(3:Nx,Ny))/dx^2 + 2*(u_m(2:Nx-1,Ny-1)-u_m(2:Nx-1,Ny))/dy^2;
        uLaplace(bottom(2:end-1)) = (u_m(bottom(2:end-1)-Ny)-2*u_m(bottom(2:end-1))+u_m(bottom(2:end-1)+Ny))/dx^2 + 2*(u_m(bottom(2:end-1)-1)-u_m(bottom(2:end-1)))/dy^2;
        uMu(bottom(2:end-1)) = u_m(bottom(2:end-1)).^3 - u_m(bottom(2:end-1)) - gamma*uLaplace(bottom(2:end-1));
        uMuLaplace(bottom(2:end-1)) = (uMu(bottom(2:end-1)-Ny)-2*uMu(bottom(2:end-1))+uMu(bottom(2:end-1)+Ny))/dx^2 + 2*(uMu(bottom(2:end-1)-1)-uMu(bottom(2:end-1)))/dy^2;

%         uLaplace(1,1) = 2*(u_m(2,1)-u_m(1,1))/dx^2 + 2*(u_m(1,2)-u_m(1,1))/dy^2;
%         uMu(1,1) = u_m(1,1).^3 - u_m(1,1) - gamma*uLaplace(1,1);
%         uMuLaplace(1,1) = 2*(uMu(2,1)-uMu(1,1))/dx^2 + 2*(uMu(1,2)-uMu(1,1))/dy^2;
%         %BC: x=0, 0<y<Y, (x_n=1, y_n=2:Ny-1), dxU=0, J.n_x=0
%         uLaplace(1,2:Ny-1) = 2*(u_m(2,2:Ny-1)-u_m(1,2:Ny-1))/dx^2 + (u_m(1,1:Ny-2)-2*u_m(1,2:Ny-1)+u_m(1,3:Ny))/dy^2;
%         uMu(1,2:Ny-1) = u_m(1,2:Ny-1).^3 - u_m(1,2:Ny-1) - gamma*uLaplace(1,2:Ny-1);
%         uMuLaplace(1,2:Ny-1) = 2*(uMu(2,2:Ny-1)-uMu(1,2:Ny-1))/dx^2 + (uMu(1,1:Ny-2)-2*uMu(1,2:Ny-1)+uMu(1,3:Ny))/dy^2;
%         %BC: x=0, y=Y, (x_n=1, y_n=Ny), dnU=0, J.n=0
%         uLaplace(1,Ny) = 2*(u_m(2,Ny)-u_m(1,Ny))/dx^2 + 2*(u_m(1,Ny-1)-u_m(1,Ny))/dy^2;
%         uMu(1,Ny) = u_m(1,Ny).^3 - u_m(1,Ny) - gamma*uLaplace(1,Ny);
%         uMuLaplace(1,Ny) = 2*(uMu(2,Ny)-uMu(1,Ny))/dx^2 + 2*(uMu(1,Ny-1)-uMu(1,Ny))/dy^2;
%         %BC: 0<x<X, y=0, (x_n=2:Nx-1, y_n=1), dyU=0, J.n_y=0
%         uLaplace(2:Nx-1,1) = (u_m(1:Nx-2,1)-2*u_m(2:Nx-1,1)+u_m(3:Nx,1))/dx^2 + 2*(u_m(2:Nx-1,2)-u_m(2:Nx-1,1))/dy^2;
%         uMu(2:Nx-1,1) = u_m(2:Nx-1,1).^3 - u_m(2:Nx-1,1) - gamma*uLaplace(2:Nx-1,1);
%         uMuLaplace(2:Nx-1,1) = (uMu(1:Nx-2,1)-2*uMu(2:Nx-1,1)+uMu(3:Nx,1))/dx^2 + 2*(uMu(2:Nx-1,2)-uMu(2:Nx-1,1))/dy^2;
%         %BC: x=X, y=0, (x_n=Nx, y_n=1), dnU=0, J.n=0
%         uLaplace(Nx,1) = 2*(u_m(Nx-1,1)-u_m(Nx,1))/dx^2 + 2*(u_m(Nx,2)-u_m(Nx,1))/dy^2;
%         uMu(Nx,1) = u_m(Nx,1).^3 - u_m(Nx,1) - gamma*uLaplace(Nx,1);
%         uMuLaplace(Nx,1) = 2*(uMu(Nx-1,1)-uMu(Nx,1))/dx^2 + 2*(uMu(Nx,2)-uMu(Nx,1))/dy^2;
%         %BC: 0<x<X, y=Y, (x_n=2:Nx-1, y_n=Ny), dyU=0, J.n_y=0
%         uLaplace(2:Nx-1,Ny) = (u_m(1:Nx-2,Ny)-2*u_m(2:Nx-1,Ny)+u_m(3:Nx,Ny))/dx^2 + 2*(u_m(2:Nx-1,Ny-1)-u_m(2:Nx-1,Ny))/dy^2;
%         uMu(2:Nx-1,Ny) = u_m(2:Nx-1,Ny).^3 - u_m(2:Nx-1,Ny) - gamma*uLaplace(2:Nx-1,Ny);
%         uMuLaplace(2:Nx-1,Ny) = (uMu(1:Nx-2,Ny)-2*uMu(2:Nx-1,Ny)+uMu(3:Nx,Ny))/dx^2 + 2*(uMu(2:Nx-1,Ny-1)-uMu(2:Nx-1,Ny))/dy^2;
%         %BC: x=X, 0<y<Y, (x_n=Nx, y_n=2:Ny-1), dxU=0, J.x=0
%         uLaplace(Nx,2:Ny-1) = 2*(u_m(Nx-1,2:Ny-1)-u_m(Nx,2:Ny-1))/dx^2 + (u_m(Nx,1:Ny-2)-2*u_m(Nx,2:Ny-1)+u_m(Nx,3:Ny))/dy^2;
%         uMu(Nx,2:Ny-1) = u_m(Nx,2:Ny-1).^3 - u_m(Nx,2:Ny-1) - gamma*uLaplace(Nx,2:Ny-1);
%         uMuLaplace(Nx,2:Ny-1) = 2*(uMu(Nx-1,2:Ny-1)-uMu(Nx,2:Ny-1))/dx^2 + (uMu(Nx,1:Ny-2)-2*uMu(Nx,2:Ny-1)+uMu(Nx,3:Ny))/dy^2;
%         %BC: x=X, y=Y, (x_n=Nx, y_n=Ny), dnU=0, J.n=0
%         uLaplace(Nx,Ny) = 2*(u_m(Nx-1,Ny)-u_m(Nx,Ny))/dx^2 + 2*(u_m(Nx,Ny-1)-u_m(Nx,Ny))/dy^2;
%         uMu(Nx,Ny) = u_m(Nx,Ny).^3 - u_m(Nx,Ny) - gamma*uLaplace(Nx,Ny);
%         uMuLaplace(Nx,Ny) = 2*(uMu(Nx-1,Ny)-uMu(Nx,Ny))/dx^2 + 2*(uMu(Nx,Ny-1)-uMu(Nx,Ny))/dy^2;

        %%%%%% EQUATIONS %%%%%%%

        uLaplace(domain_inds) = (u_m(domain_inds-Ny)-2*u_m(domain_inds)+u_m(domain_inds+Ny))/dx^2 + (u_m(domain_inds-1)-2*u_m(domain_inds)+u_m(domain_inds+1))/dy^2;
        uMu(domain_inds) = u_m(domain_inds).^3 - u_m(domain_inds) - gamma*uLaplace(domain_inds);
        uMuLaplace(domain_inds) = (uMu(domain_inds-Ny)-2*uMu(domain_inds)+uMu(domain_inds+Ny))/dx^2 + (uMu(domain_inds-1)-2*uMu(domain_inds)+uMu(domain_inds+1))/dy^2;

        %uLaplace(2:Nx-1,2:Ny-1) = (u_m(1:Nx-2,2:Ny-1)-2*u_m(2:Nx-1,2:Ny-1)+u_m(3:Nx,2:Ny-1))/dx^2 + (u_m(2:Nx-1,1:Ny-2)-2*u_m(2:Nx-1,2:Ny-1)+u_m(2:Nx-1,3:Ny))/dy^2;
        %uMu(2:Nx-1,2:Ny-1) = u_m(2:Nx-1,2:Ny-1).^3 - u_m(2:Nx-1,2:Ny-1) - gamma*uLaplace(2:Nx-1,2:Ny-1);
        %uMuLaplace(2:Nx-1,2:Ny-1) = (uMu(1:Nx-2,2:Ny-1)-2*uMu(2:Nx-1,2:Ny-1)+uMu(3:Nx,2:Ny-1))/dx^2 + (uMu(2:Nx-1,1:Ny-2)-2*uMu(2:Nx-1,2:Ny-1)+uMu(2:Nx-1,3:Ny))/dy^2;
        
        %uLaplace(domain_inds,:) = (u_m(domain_inds-Ny,:)-2*u_m(domain_inds,:)+u_m(domain_inds+Ny,:))/dx^2 + (u_m(domain_inds-1,:)-2*u_m(domain_inds,:)+u_m(domain_inds+1,:))/dy^2;
        %uMu(domain_inds,:) = u_m(domain_inds,:).^3 - u_m(domain_inds,:) - gamma*uLaplace(domain_inds,:);
        %uMuLaplace(domain_inds,:) = (uMu(domain_inds-Ny,:)-2*uMu(domain_inds,:)+uMu(domain_inds+Ny,:))/dx^2 + (uMu(domain_inds-1,:)-2*uMu(domain_inds,:)+uMu(domain_inds+1,:))/dy^2;

%         uLaplace(cavityDomain(:,1),cavityDomain(:,2)) = (u_m(cavityDomain(:,1)-1,cavityDomain(:,2))-2*u_m(cavityDomain(:,1),cavityDomain(:,2))+u_m(cavityDomain(:,1)+1,cavityDomain(:,2)))/dx^2 + ...
%                                                         (u_m(cavityDomain(:,1),cavityDomain(:,2)-1)-2*u_m(cavityDomain(:,1),cavityDomain(:,2))+u_m(cavityDomain(:,1),cavityDomain(:,2)+1))/dy^2;
%         uMu(cavityDomain(:,1),cavityDomain(:,2)) = u_m(cavityDomain(:,1),cavityDomain(:,2)).^3 - u_m(cavityDomain(:,1),cavityDomain(:,2)) - gamma*uLaplace(cavityDomain(:,1),cavityDomain(:,2));
%         uMuLaplace(cavityDomain(:,1),cavityDomain(:,2)) = (uMu(cavityDomain(:,1)-1,cavityDomain(:,2))-2*uMu(cavityDomain(:,1),cavityDomain(:,2))+uMu(cavityDomain(:,1)+1,cavityDomain(:,2)))/dx^2 + ...
%                                                           (uMu(cavityDomain(:,1),cavityDomain(:,2)-1)-2*uMu(cavityDomain(:,1),cavityDomain(:,2))+uMu(cavityDomain(:,1),cavityDomain(:,2)+1))/dy^2;

        %uLaplace(perimeter(:,1),perimeter(:,2)) = 2*(u_m(perimeter(:,1)-1,perimeter(:,2))-u_m(perimeter(:,1),perimeter(:,2)))/dx^2 + 2*(u_m(perimeter(:,1),perimeter(:,2)-1)-u_m(perimeter(:,1),perimeter(:,2)))/dy^2;
        %uMu(perimeter(:,1),perimeter(:,2)) = u_m(perimeter(:,1),perimeter(:,2)).^3 - u_m(perimeter(:,1),perimeter(:,2)) - gamma*uLaplace(perimeter(:,1),perimeter(:,2));
        %uMuLaplace(perimeter(:,1),perimeter(:,2)) = 2*(uMu(perimeter(:,1)-1,perimeter(:,2))-uMu(perimeter(:,1),perimeter(:,2)))/dx^2 + 2*(uMu(perimeter(:,1),perimeter(:,2)-1)-uMu(perimeter(:,1),perimeter(:,2)))/dy^2;
        
        %pos = [x(perim_inds(:,1)),y(perim_inds(:,2))];
        %xPos = pos(:,1)./sqrt(pos(:,1).^2 + pos(:,2).^2);
        %yPos = pos(:,2)./sqrt(pos(:,1).^2 + pos(:,2).^2);
        %u_m_ip1j = (u_m(perimeter-1).*yPos - u_m(perimeter+1).*yPos + u_m(perimeter-Ny).*xPos)./(xPos);
        %u_m_in1j = (u_m(perimeter+Ny).*xPos + u_m(perimeter+1).*yPos - u_m(perimeter-1).*yPos)./(xPos);
        %u_m_ijp1 = (u_m(perimeter-1).*yPos + u_m(perimeter-Ny).*xPos - u_m(perimeter+Ny).*xPos)./(yPos);
        %u_m_ijn1 = (u_m(perimeter+Ny).*xPos - u_m(perimeter-Ny).*xPos + u_m(perimeter+1).*yPos)./(yPos);


        %uLaplace(perimeter) = (u_m(perimeter-Ny)-2*u_m(perimeter)+u_m(perimeter+Ny))/dx^2 + (u_m(perimeter-1)-2*u_m(perimeter)+u_m(perimeter+1))/dy^2;
        %u_m_ijp1
        %uLaplace(perimeter) = (u_m_in1j-2*u_m(perimeter)+u_m_ip1j)/dx^2 + (u_m_ijn1-2*u_m(perimeter)+u_m_ijp1)/dy^2;
        %uMu(perimeter) = u_m(perimeter).^3 - u_m(perimeter) - gamma*uLaplace(perimeter);
        %uMuLaplace(perimeter) = (uMu(perimeter-Ny)-2*uMu(perimeter)+uMu(perimeter+Ny))/dx^2 + (uMu(perimeter-1)-2*uMu(perimeter)+uMu(perimeter+1))/dy^2;
        %pos = [x(perimeter(:,1)),y(perimeter(:,2))];
        %xPos = pos(:,1)./sqrt(pos(:,1).^2 + pos(:,2).^2);
        %yPos = pos(:,2)./sqrt(pos(:,1).^2 + pos(:,2).^2);
        %u_m_ip1j = (u_m(perimeter(:,1),perimeter(:,2)-1).*yPos - u_m(perimeter(:,1),perimeter(:,2)+1).*yPos + u_m(perimeter(:,1)-1,perimeter(:,2)).*xPos)./(xPos);
        %u_m_in1j = (u_m(perimeter(:,1)+1,perimeter(:,2)).*xPos + u_m(perimeter(:,1),perimeter(:,2)+1).*yPos - u_m(perimeter(:,1),perimeter(:,2)-1).*yPos)./(xPos);
        %u_m_ijp1 = (u_m(perimeter(:,1),perimeter(:,2)-1).*yPos + u_m(perimeter(:,1)-1,perimeter(:,2)).*xPos - u_m(perimeter(:,1)+1,perimeter(:,2)).*xPos)./(yPos);
        %u_m_ijn1 = (u_m(perimeter(:,1)+1,perimeter(:,2)).*xPos - u_m(perimeter(:,1)-1,perimeter(:,2)).*xPos + u_m(perimeter(:,1),perimeter(:,2)+1).*yPos)./(yPos);

        %NONE
        % uLaplace(perimeter(:,1),perimeter(:,2)) = (u_m(perimeter(:,1)-1,perimeter(:,2))-2*u_m(perimeter(:,1),perimeter(:,2))+u_m(perimeter(:,1)+1,perimeter(:,2)))/dx^2 + (u_m(perimeter(:,1),perimeter(:,2)-1)-2*u_m(perimeter(:,1),perimeter(:,2))+u_m(perimeter(:,1),perimeter(:,2)+1))/dy^2;
        % u_m_ijp1
        % uLaplace(perimeter(:,1),perimeter(:,2)) = (u_m(perimeter(:,1)-1,perimeter(:,2))-2*u_m(perimeter(:,1),perimeter(:,2))+u_m(perimeter(:,1)+1,perimeter(:,2)))/dx^2 + (u_m(perimeter(:,1),perimeter(:,2)-1)-2*u_m(perimeter(:,1),perimeter(:,2))+u_m_ijp1)/dy^2;
        % u_m_ijn1 % DOES NOT FINISH 
        % uLaplace(perimeter(:,1),perimeter(:,2)) = (u_m(perimeter(:,1)-1,perimeter(:,2))-2*u_m(perimeter(:,1),perimeter(:,2))+u_m(perimeter(:,1)+1,perimeter(:,2)))/dx^2 + (u_m_ijn1-2*u_m(perimeter(:,1),perimeter(:,2))+u_m(perimeter(:,1),perimeter(:,2)+1))/dy^2;
        % u_m_ip1j %
        % uLaplace(perimeter(:,1),perimeter(:,2)) = (u_m(perimeter(:,1)-1,perimeter(:,2))-2*u_m(perimeter(:,1),perimeter(:,2))+u_m_ip1j)/dx^2 + (u_m(perimeter(:,1),perimeter(:,2)-1)-2*u_m(perimeter(:,1),perimeter(:,2))+u_m(perimeter(:,1),perimeter(:,2)+1))/dy^2;
        % u_m_ip1j, u_m_ijn1
        % uLaplace(perimeter(:,1),perimeter(:,2)) = (u_m_ip1j-2*u_m(perimeter(:,1),perimeter(:,2))+u_m(perimeter(:,1)-1,perimeter(:,2)))/dx^2 + (u_m(perimeter(:,1),perimeter(:,2)+1)-2*u_m(perimeter(:,1),perimeter(:,2))+u_m_ijn1)/dy^2;

        %uMu(perimeter(:,1),perimeter(:,2)) = u_m(perimeter(:,1),perimeter(:,2)).^3 - u_m(perimeter(:,1),perimeter(:,2)) - gamma*uLaplace(perimeter(:,1),perimeter(:,2));
        %uMuLaplace(perimeter(:,1),perimeter(:,2)) = (uMu(perimeter(:,1)+1,perimeter(:,2))-2*uMu(perimeter(:,1),perimeter(:,2))+uMu(perimeter(:,1)-1,perimeter(:,2)))/dx^2 + (uMu(perimeter(:,1),perimeter(:,2)+1)-2*uMu(perimeter(:,1),perimeter(:,2))+uMu(perimeter(:,1),perimeter(:,2)-1))/dy^2;

        %uLaplace(perimeter(:,1),perimeter(:,2)) = 2*(u_m(perimeter(:,1)-1,perimeter(:,2))-u_m(perimeter(:,1),perimeter(:,2)))/dx^2 + 2*(u_m(perimeter(:,1),perimeter(:,2)-1)-u_m(perimeter(:,1),perimeter(:,2)))/dy^2;
        %uMu(perimeter(:,1),perimeter(:,2)) = u_m(perimeter(:,1),perimeter(:,2)).^3 - u_m(perimeter(:,1),perimeter(:,2)) - gamma*uLaplace(perimeter(:,1),perimeter(:,2));
        %uMuLaplace(perimeter(:,1),perimeter(:,2)) = 2*(uMu(perimeter(:,1)-1,perimeter(:,2))-uMu(perimeter(:,1),perimeter(:,2)))/dx^2 + 2*(uMu(perimeter(:,1),perimeter(:,2)-1)-uMu(perimeter(:,1),perimeter(:,2)))/dy^2;
        %uLaplace(perimeter) = (2*u_m(iPerimeter)-2*u_m(perimeter))/dr^2;

        %uMu(perimeter) = u_m(perimeter).^3 - u_m(perimeter) - gamma*uLaplace(perimeter);
        %uMuLaplace(perimeter) = (2*uMu(iPerimeter)-2*uMu(perimeter))/dr^2;
        
        duT = D*uMuLaplace;

        %%%%%% RESHAPING OUTPUT TO COLUMN VECTOR %%%%%%%
        %out(1:length(x)*length(y),1) = reshape(duT,[],1);
        out(1:length(x)*length(y),1) = reshape(duT,[],1);
        %out = duT;

    end 