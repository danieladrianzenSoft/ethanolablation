function cryptFunc()

global Df De Ds phiGE phiES c0 w kd kb kl t Cout ratio phie phis kon koff rc rs he hs ncr ner nez nsr nsz  dr dz rgrid zgrid cryptOpen tnoflux


init = zeros(1,length(rgrid)*length(zgrid));


tic;

prevTime = 0;

[~,Cout] = ode45(@dcdt,t,init);


    function out = dcdt(tt,cc)
        
        if (prevTime < floor(tt/60))
            prevTime = floor(tt/60);
            fprintf('Time in minutes: %0.0f\n',tt/60);
            toc;
        end
        
        out =  zeros(length(rgrid)*length(zgrid),1);
        cmat = reshape(cc(1:length(rgrid)*length(zgrid)),length(rgrid),length(zgrid));
        rhere = rgrid'*ones(1,length(zgrid));
        
        Dhere = De*ones(size(cmat));
        Dhere(1:ncr,:) = Df;
        Dhere(ncr+ner+1:end,nez+1:end) = Ds;
        
        Khere = zeros(size(cmat));
        Khere(ncr+ner+1:end,nez+1:end) = kb;
        
        Cu = [cmat(2,:);cmat(1:end-1,:)];
        Cu(ncr+1,:) = (Df*cmat(ncr-1,:)/dr+(De*cmat(ncr+2,:)/dr))/((Df/phiGE)/dr+De/dr);
        Cu(ncr+ner+1,nez+1:end) = (De*cmat(ncr+ner-1,nez+1:end)/dr+(Ds*cmat(ncr+ner+2,nez+1:end)/dr))/((De/phiES)/dr+Ds/dr);
        Cu(ncr+ner+1,1:nez) = (De*cmat(ncr+ner-1,1:nez)/dr+De*cmat(ncr+ner+2,1:nez)/dr)/((De)/dr+De/dr);
        
        Cd = [cmat(2:end,:);cmat(end-1,:)];
        Cd(ncr,:) = (Df*cmat(ncr-1,:)/dr+(De*cmat(ncr+2,:)/dr))/(Df/dr+De*phiGE/dr);
        Cd(ncr+ner,nez+1:end) = (De*cmat(ncr+ner-1,nez+1:end)/dr+(Ds*cmat(ncr+ner+2,nez+1:end)/dr))/(De/dr+Ds*phiES/dr);
        Cd(ncr+ner,1:nez) = (De*cmat(ncr+ner-1,1:nez)/dr+(De*cmat(ncr+ner+2,1:nez)/dr))/(De/dr+De/dr);
         
        Cl = [[c0*ones(ncr,1);c0*phiGE*ones(ner+nsr,1)],cmat(:,1:end-1)];
        Cl(ncr+ner+1:end,nez+1) = (De*cmat(ncr+ner+1:end,nez-1)/dz+(Ds*cmat(ncr+ner+1:end,nez+2)/dz))/((De/phiES)/dz+Ds/dz);
        Cl(ncr+1:ncr+ner,nez+1) = (De*cmat(ncr+1:ncr+ner,nez-1)/dz+(De*cmat(ncr+1:ncr+ner,nez+2)/dz))/((De)/dz+De/dz);
        Cl(1:ncr,nez+1) = (Df*cmat(1:ncr,nez-1)/dz+(Df*cmat(1:ncr,nez+2)/dz))/((Df)/dz+Df/dz);

        Cl(1:ncr,1) = cryptOpen*(c0*ones(ncr,1)-cmat(1:ncr,2))+cmat(1:ncr,2); 

        if(tt > tnoflux)
           Cl(1:end,1) =  cmat(1:end,2); 
        end

        
        Cr = [cmat(:,2:end),cmat(:,end-1)];
        Cr(ncr+ner+1:end,nez) = (De*cmat(ncr+ner+1:end,nez-1)/dz+(Ds*cmat(ncr+ner+1:end,nez+2)/dz))/(De/dz+(Ds*phiES)/dz);
        Cr(ncr+1:ncr+ner,nez) = (De*cmat(ncr+1:ncr+ner,nez-1)/dz+(De*cmat(ncr+1:ncr+ner,nez+2)/dz))/(De/dz+(De)/dz);
        Cr(1:ncr,nez) = (Df*cmat(1:ncr,nez-1)/dz+(Df*cmat(1:ncr,nez+2)/dz))/(Df/dz+Df/dz);
        
        Ru = rhere-dr/2;
        Ru(1,:) = Ru(2,:);
        Rd = rhere+dr/2;
        Rd(end,:) = Rd(end-1,:);

        
        dcTFV = Dhere.*((Ru./rhere.*(Cu-cmat)+Rd./rhere.*(Cd-cmat))/dr/dr...
            +(Cl-2*cmat+Cr)/dz/dz)-Khere.*cmat;
        
        
        out(1:length(rgrid)*length(zgrid)) = reshape(dcTFV,[],1);
        

    end

end