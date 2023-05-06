function out = cahn_hilliard_1Dspherical_rc(t,u,r,rc_spline,D,gamma) 

        ruDuDr = zeros(numel(r), size(u,2));
        ruLaplace = zeros(numel(r),size(u,2));
        ruMu = zeros(numel(r),size(u,2));
        ruMuLaplace = zeros(numel(r),size(u,2));
        
        dr = r(2)-r(1);
        Nr = length(r);
        ru_m = u;

        %%%%%% BOUNDARY CONDITIONS %%%%%%%

        %BC: r=0, dc/dr = 0, (r_n=1)
        ruLaplace(1,:) = 2*(ru_m(2,:)-ru_m(1,:)) / dr^2; % d2C/dr2 = 2Ci+1/dr2
        ruMu(1,:) = ru_m(1,:).^3 - ru_m(1,:) - gamma*ruLaplace(1,:); % dC/dt = c^3 - c - gamma*gradc
        ruMuLaplace(1,:) = 2*(ruMu(2,:)-ruMu(1,:)) / dr^2; % d2C/dr2 = 2Ci+1/dr2

        %BC: r=Rtot, dc/dr = 0, (r_n=R)
        ruLaplace(Nr,:) = 2*(ru_m(Nr-1,:)-ru_m(Nr,:)) / dr^2; % d2C/dr2 = 2Ci+1/dr2
        ruMu(Nr,:) = ru_m(Nr,:).^3 - ru_m(Nr,:) - gamma*ruLaplace(Nr,:); % dC/dt = c^3 - c - gamma*gradc
        ruMuLaplace(Nr,:) = 2*(ruMu(Nr-1,:)-ruMu(Nr,:)) / dr^2; % d2C/dr2 = 2Ci+1/dr2
       
        %%%%%% EQUATIONS %%%%%%%

        rc = ppval(rc_spline,t);
        interior_bound = find(r>=rc,1)-1;
        domain_inds = 2:interior_bound-2;

        ruDuDr(domain_inds,:) = (ru_m(domain_inds+1,:) - ru_m(domain_inds-1,:)) / (2*dr); % dC/dr = (Ci+1 - Ci-1) / 2dr
        ruLaplace(domain_inds,:) = 2*(ruDuDr(domain_inds,:)./r(domain_inds)') + (ru_m(domain_inds-1,:) - 2*ru_m(domain_inds,:) + ru_m(domain_inds+1,:)) / dr^2;  % nabla2C = 2/r dC/dr + d2C/dr2
      
        ruMu(domain_inds,:) = ru_m(domain_inds,:).^3 - ru_m(domain_inds,:) - gamma*ruLaplace(domain_inds,:);
        ruMuDuDr(domain_inds,:) = (ruMu(domain_inds+1,:) - ruMu(domain_inds-1,:)) / (2*dr);
        ruMuLaplace(domain_inds,:) = 2*ruMuDuDr(domain_inds,:)./r(domain_inds)' + (ruMu(domain_inds-1,:) - 2*ruMu(domain_inds,:) + ru_m(domain_inds+1,:)) / dr^2;
        
        IBC = "ConstantConcentration";

        if (IBC == "NoFlux")
            %ruDuDr(interior_bound,:) = (ru_m(interior_bound+1,:) - ru_m(interior_bound-1,:)) / (2*dr); % dC/dr = (Ci+1 - Ci-1) / 2dr
            ruLaplace(interior_bound,:) = (2*ru_m(interior_bound-1,:) - 2*ru_m(interior_bound,:)) / dr^2;  % nabla2C = 2/r dC/dr + d2C/dr2
          
            ruMu(interior_bound,:) = ru_m(interior_bound,:).^3 - ru_m(interior_bound,:) - gamma*ruLaplace(interior_bound,:);
            ruMuDuDr(interior_bound,:) = (ruMu(interior_bound+1,:) - ruMu(interior_bound-1,:)) / (2*dr);
            %ruMuLaplace(interior_bound,:) = (2*ruMu(interior_bound-1,:) - 2*ruMu(interior_bound,:)) / dr^2;
            ruMuLaplace(interior_bound,:) = 2*ruMuDuDr(interior_bound,:)./r(interior_bound)' + (ruMu(interior_bound-1,:) - 2*ruMu(interior_bound,:) + ru_m(interior_bound+1,:)) / dr^2;

        elseif (IBC == "Flux")
            K = 3e-10; %m^2 / (kPa s) = cm^2 /(barye s)
            S = 4*pi*0.3^2; % surface area of cavity cm^2;
            pi_diff = 26.6664e3; %=20mmHg in interstitium of different tumours, stohrer et al, 2000
            J0 = -K*S*pi_diff;
            %J0 = 0.05;

            ruDuDr(interior_bound,:) = J0; % dC/dr = (Ci+1 - Ci-1) / 2dr
            ruLaplace(interior_bound,:) = 2*(ruDuDr(interior_bound,:)./r(interior_bound)') + (2*ru_m(interior_bound-1,:) - 2*ru_m(interior_bound,:) + 2*J0*dr) / dr^2;  % nabla2C = 2/r dC/dr + d2C/dr2
          
            ruMu(interior_bound,:) = ru_m(interior_bound,:).^3 - ru_m(interior_bound,:) - gamma*ruLaplace(interior_bound,:);
            ruMuDuDr(interior_bound,:) = (ruMu(interior_bound+1,:) - ruMu(interior_bound-1,:)) / (2*dr);
            ruMuLaplace(interior_bound,:) = 2*ruMuDuDr(interior_bound,:)./r(interior_bound)' + (ruMu(interior_bound-1,:) - 2*ruMu(interior_bound,:) + ru_m(interior_bound+1,:)) / dr^2;

        elseif (IBC == "ConstantConcentration")
            C0 = 1;
            ruDuDr(interior_bound,:) = (C0 - ru_m(interior_bound-1,:)) / (2*dr); % dC/dr = (Ci+1 - Ci-1) / 2dr
            ruLaplace(interior_bound,:) = 2*ruDuDr(interior_bound,:)./r(interior_bound)' + (ru_m(interior_bound-1,:) - 2*ru_m(interior_bound) + C0) / dr^2;  % nabla2C = 2/r dC/dr + d2C/dr2
          
            ruMu(interior_bound,:) = ru_m(interior_bound,:).^3 - ru_m(interior_bound,:) - gamma*ruLaplace(interior_bound,:);
            ruMuDuDr(interior_bound,:) = (ruMu(interior_bound+1,:) - ruMu(interior_bound-1,:)) / (2*dr);
            %ruMuLaplace(interior_bound,:) = (2*ruMu(interior_bound-1,:) - 2*ruMu(interior_bound,:)) / dr^2;
            ruMuLaplace(interior_bound,:) = 2*ruMuDuDr(interior_bound,:)./r(interior_bound)' + (ruMu(interior_bound-1,:) - 2*ruMu(interior_bound,:) + ru_m(interior_bound+1,:)) / dr^2;

        end

        rduT = D*ruMuLaplace;

        %%%%%% RESHAPING OUTPUT TO COLUMN VECTOR %%%%%%%

        out = rduT;

    end 