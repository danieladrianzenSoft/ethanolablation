%2D cartesian coordinate spinodal decomposition
function out = cahn_hilliard_2D_vec(t,u,x,y,interior_bounds,domain_inds,D,gamma) 
    %function out = cahn_hilliard(t,u,x,y,perimeter,oPerimeter,iPerimeter,cavityDomain,domain_inds,D,gamma) 
        uLaplace = zeros(numel(x)*numel(y),size(u,2));
        uMu = zeros(numel(x)*numel(y),size(u,2));
        uMuLaplace = zeros(numel(x)*numel(y),size(u,2));
        dx = x(2)-x(1);
        dy = y(2)-y(1);
        Nx = length(x);
        Ny = length(y);

        %%%%%% RESHAPING INPUT TO 2X2 MATRIX %%%%%%%
        u_m = u;

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
        uLaplace(left(1),:) = 2*(u_m(left(1)+Ny,:)-u_m(left(1),:))/dx^2 + 2*(u_m(left(1)+1,:)-u_m(left(1),:))/dy^2;
        uMu(left(1),:) = u_m(left(1),:).^3 - u_m(left(1),:) - gamma*uLaplace(left(1),:);
        uMuLaplace(left(1),:) = 2*(uMu(left(1)+Ny,:)-uMu(left(1),:))/dx^2 + 2*(uMu(left(1)+1,:)-uMu(left(1),:))/dy^2;
        % top right corner
        %uLaplaceTest(Nx,1) = 2*(u_m(Nx-1,1)-u_m(Nx,1))/dx^2 + 2*(u_m(Nx,2)-u_m(Nx,1))/dy^2;
        uLaplace(right(1),:) = 2*(u_m(right(1)-Ny,:)-u_m(right(1),:))/dx^2 + 2*(u_m(right(1)+1,:)-u_m(right(1),:))/dy^2;
        uMu(right(1),:) = u_m(right(1),:).^3 - u_m(right(1),:) - gamma*uLaplace(right(1),:);
        uMuLaplace(right(1),:) = 2*(uMu(right(1)-Ny,:)-uMu(right(1),:))/dx^2 + 2*(uMu(right(1)+1,:)-uMu(right(1),:))/dy^2;
        % bottom left corner
        %uLaplaceTest(1,Ny) = 2*(u_m(2,Ny)-u_m(1,Ny))/dx^2 + 2*(u_m(1,Ny-1)-u_m(1,Ny))/dy^2;
        uLaplace(left(end),:) = 2*(u_m(left(end)+Ny,:)-u_m(left(end),:))/dx^2 + 2*(u_m(left(end)-1,:)-u_m(left(end),:))/dy^2;
        uMu(left(end),:) = u_m(left(end),:).^3 - u_m(left(end),:) - gamma*uLaplace(left(end),:);
        uMuLaplace(left(end),:) = 2*(uMu(left(end)+Ny,:)-uMu(left(end),:))/dx^2 + 2*(uMu(left(end)-1,:)-uMu(left(end),:))/dy^2;
        % bottom right corner
        %uLaplaceTest(Nx,Ny) = 2*(u_m(Nx-1,Ny)-u_m(Nx,Ny))/dx^2 + 2*(u_m(Nx,Ny-1)-u_m(Nx,Ny))/dy^2;
        uLaplace(right(end),:) = 2*(u_m(right(end)-Ny,:)-u_m(right(end),:))/dx^2 + 2*(u_m(right(end)-1,:)-u_m(right(end),:))/dy^2;
        uMu(right(end),:) = u_m(right(end)).^3 - u_m(right(end),:) - gamma*uLaplace(right(end),:);
        uMuLaplace(right(end),:) = 2*(uMu(right(end)-Ny,:)-uMu(right(1),:))/dx^2 + 2*(uMu(right(end)-1,:)-uMu(right(end),:))/dy^2;
        % top
        %uLaplaceTest(2:Nx-1,1) = (u_m(1:Nx-2,1)-2*u_m(2:Nx-1,1)+u_m(3:Nx,1))/dx^2 + 2*(u_m(2:Nx-1,2)-u_m(2:Nx-1,1))/dy^2;
        uLaplace(top(2:end-1),:) = (u_m(top(2:end-1)-Ny,:)-2*u_m(top(2:end-1),:)+u_m(top(2:end-1)+Ny,:))/dx^2+2*(u_m(top(2:end-1)+1,:)-u_m(top(2:end-1),:))/dy^2;
        uMu(top(2:end-1),:) = u_m(top(2:end-1),:).^3 - u_m(top(2:end-1),:) - gamma*uLaplace(top(2:end-1),:);
        uMuLaplace(top(2:end-1),:) = (uMu(top(2:end-1)-Ny,:)-2*uMu(top(2:end-1),:)+uMu(top(2:end-1)+Ny,:))/dx^2+2*(uMu(top(2:end-1)+1,:)-uMu(top(2:end-1),:))/dy^2;
        % left
        %uLaplaceTest(1,2:Ny-1) = 2*(u_m(2,2:Ny-1)-u_m(1,2:Ny-1))/dx^2 + (u_m(1,1:Ny-2)-2*u_m(1,2:Ny-1)+u_m(1,3:Ny))/dy^2;
        uLaplace(left(2:end-1),:) = 2*(u_m(left(2:end-1)+Ny,:)-u_m(left(2:end-1),:))/dx^2+(u_m(left(2:end-1)-1,:)-2*u_m(left(2:end-1),:)+u_m(left(2:end-1)+1,:))/dy^2;
        uMu(left(2:end-1),:) = u_m(left(2:end-1),:).^3 - u_m(left(2:end-1),:) - gamma*uLaplace(left(2:end-1),:);
        uMuLaplace(left(2:end-1),:) = 2*(uMu(left(2:end-1)+Ny,:)-uMu(left(2:end-1),:))/dx^2+(uMu(left(2:end-1)-1,:)-2*uMu(left(2:end-1),:)+uMu(left(2:end-1)+1,:))/dy^2;
        % right
        %uLaplaceTest(Nx,2:Ny-1) = 2*(u_m(Nx-1,2:Ny-1)-u_m(Nx,2:Ny-1))/dx^2 + (u_m(Nx,1:Ny-2)-2*u_m(Nx,2:Ny-1)+u_m(Nx,3:Ny))/dy^2;
        uLaplace(right(2:end-1),:) = 2*(u_m(right(2:end-1)-Ny,:)-u_m(right(2:end-1),:))/dx^2 + (u_m(right(2:end-1)-1,:)-2*u_m(right(2:end-1),:)+u_m(right(2:end-1)+1,:))/dy^2;
        uMu(right(2:end-1),:) = u_m(right(2:end-1),:).^3 - u_m(right(2:end-1),:) - gamma*uLaplace(right(2:end-1),:);
        uMuLaplace(right(2:end-1),:) = 2*(uMu(right(2:end-1)-Ny,:)-uMu(right(2:end-1),:))/dx^2 + (uMu(right(2:end-1)-1,:)-2*uMu(right(2:end-1),:)+uMu(right(2:end-1)+1,:))/dy^2;
        % bottom
        %uLaplaceTest(2:Nx-1,Ny) = (u_m(1:Nx-2,Ny)-2*u_m(2:Nx-1,Ny)+u_m(3:Nx,Ny))/dx^2 + 2*(u_m(2:Nx-1,Ny-1)-u_m(2:Nx-1,Ny))/dy^2;
        uLaplace(bottom(2:end-1),:) = (u_m(bottom(2:end-1)-Ny,:)-2*u_m(bottom(2:end-1),:)+u_m(bottom(2:end-1)+Ny,:))/dx^2 + 2*(u_m(bottom(2:end-1)-1,:)-u_m(bottom(2:end-1),:))/dy^2;
        uMu(bottom(2:end-1),:) = u_m(bottom(2:end-1),:).^3 - u_m(bottom(2:end-1),:) - gamma*uLaplace(bottom(2:end-1),:);
        uMuLaplace(bottom(2:end-1),:) = (uMu(bottom(2:end-1)-Ny,:)-2*uMu(bottom(2:end-1),:)+uMu(bottom(2:end-1)+Ny,:))/dx^2 + 2*(uMu(bottom(2:end-1)-1,:)-uMu(bottom(2:end-1),:))/dy^2;

        %%%%%% EQUATIONS %%%%%%%

        uLaplace(domain_inds,:) = (u_m(domain_inds-Ny,:)-2*u_m(domain_inds,:)+u_m(domain_inds+Ny,:))/dx^2 + (u_m(domain_inds-1,:)-2*u_m(domain_inds,:)+u_m(domain_inds+1,:))/dy^2;
        uMu(domain_inds,:) = u_m(domain_inds,:).^3 - u_m(domain_inds,:) - gamma*uLaplace(domain_inds,:);
        uMuLaplace(domain_inds,:) = (uMu(domain_inds-Ny,:)-2*uMu(domain_inds,:)+uMu(domain_inds+Ny,:))/dx^2 + (uMu(domain_inds-1,:)-2*uMu(domain_inds,:)+uMu(domain_inds+1,:))/dy^2;

        IBC = "NoFlux";

        if (IBC == "Flux")
            % DUDR.(-n) = J0
            % Jv = KS([Pc-Pi]-sigma[PIp-PIg])
            K = 3e-10; %m^2 / (kPa s) = cm^2 /(barye s)
            S = 4*pi*0.3^2; % surface area of cavity cm^2;
            pi_diff = 26.6664e3; %=20mmHg in interstitium of different tumours, stohrer et al, 2000
            J0 = -K*S*pi_diff;
            %J0 = -0.05;
            [yPos_ind,xPos_ind] = ind2sub([Nx,Ny],interior_bounds);
            yPos_ind = Ny - yPos_ind + 1;
            xPos = x(xPos_ind)'./(sqrt(x(xPos_ind).^2 + y(yPos_ind).^2))';
            yPos = y(yPos_ind)'./(sqrt(x(xPos_ind).^2 + y(yPos_ind).^2))';
            u_m_ip1j = (u_m(interior_bounds+1,:).*yPos - u_m(interior_bounds-1,:).*yPos + u_m(interior_bounds-Ny,:).*xPos - J0)./(xPos);
            u_m_in1j = (u_m(interior_bounds+Ny,:).*xPos + u_m(interior_bounds-1,:).*yPos - u_m(interior_bounds+1,:).*yPos + J0)./(xPos);
            u_m_ijp1 = (u_m(interior_bounds+1,:).*yPos + u_m(interior_bounds-Ny,:).*xPos - u_m(interior_bounds+Ny,:).*xPos - J0)./(yPos);
            u_m_ijn1 = (u_m(interior_bounds+Ny,:).*xPos - u_m(interior_bounds-Ny,:).*xPos + u_m(interior_bounds-1,:).*yPos + J0)./(yPos);
            uLaplace(interior_bounds,:) = (u_m_in1j-2*u_m(interior_bounds,:)+u_m_ip1j)/dx^2 + (u_m_ijn1-2*u_m(interior_bounds,:)+u_m_ijp1)/dy^2;
            uMu(interior_bounds,:) = u_m(interior_bounds,:).^3 - u_m(interior_bounds,:) - gamma*uLaplace(interior_bounds,:);
            uMuLaplace(interior_bounds,:) = (uMu(interior_bounds-Ny,:)-2*uMu(interior_bounds,:)+uMu(interior_bounds+Ny,:))/dx^2 + (uMu(interior_bounds-1,:)-2*uMu(interior_bounds,:)+uMu(interior_bounds+1,:))/dy^2;
        elseif (IBC == "NoFlux")
            [yPos_ind,xPos_ind] = ind2sub([Nx,Ny],interior_bounds);
            yPos_ind = Ny - yPos_ind + 1;
            xPos = x(xPos_ind)'./(sqrt(x(xPos_ind).^2 + y(yPos_ind).^2))';
            yPos = y(yPos_ind)'./(sqrt(x(xPos_ind).^2 + y(yPos_ind).^2))';
            u_m_ip1j = (u_m(interior_bounds+1,:).*yPos - u_m(interior_bounds-1,:).*yPos + u_m(interior_bounds-Ny,:).*xPos)./(xPos);
            u_m_in1j = (u_m(interior_bounds+Ny,:).*xPos + u_m(interior_bounds-1,:).*yPos - u_m(interior_bounds+1,:).*yPos)./(xPos);
            u_m_ijp1 = (u_m(interior_bounds+1,:).*yPos + u_m(interior_bounds-Ny,:).*xPos - u_m(interior_bounds+Ny,:).*xPos)./(yPos);
            u_m_ijn1 = (u_m(interior_bounds+Ny,:).*xPos - u_m(interior_bounds-Ny,:).*xPos + u_m(interior_bounds-1,:).*yPos)./(yPos);
            uLaplace(interior_bounds,:) = (u_m_in1j-2*u_m(interior_bounds,:)+u_m_ip1j)/dx^2 + (u_m_ijn1-2*u_m(interior_bounds,:)+u_m_ijp1)/dy^2;
            uMu(interior_bounds,:) = u_m(interior_bounds,:).^3 - u_m(interior_bounds,:) - gamma*uLaplace(interior_bounds,:);
            uMuLaplace(interior_bounds,:) = (uMu(interior_bounds-Ny,:)-2*uMu(interior_bounds,:)+uMu(interior_bounds+Ny,:))/dx^2 + (uMu(interior_bounds-1,:)-2*uMu(interior_bounds,:)+uMu(interior_bounds+1,:))/dy^2;
        end


        duT = D*uMuLaplace;

        out = duT;

    end 