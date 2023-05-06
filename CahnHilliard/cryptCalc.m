
global Df De Ds phiGE phiES c0 w kd kb kl t Cout ratio phie phis kon koff rc rs he hs ncr ner nez nsr nsz  dr dz rgrid zgrid cryptOpen tnoflux

Df = 6e-6;
De = 7e-8*2;
Ds = 4e-7;
phiGE = 0.75;
phiES = 1;

c0 = 1e7;
V = 100;
w = 5;
kd = 1.22/3600;
kb = 0.119/3600;
kl = 1.41/3600;


ratio = 0.1;
phie = 0.95;
phis = 0.1;
kon = 0.693/3600;
koff = 0.00413/3600;

t = linspace(0,10*60,100);


rc = 20/10000;      %crypt radii 20 micron
rs = 20/10000;      %stroma radii 7.5 - 30 micron (crypt spacing 85-130 micron)
he = 15/10000;      %epithelial thickness 15-25 micron
hs = 1000/10000;    %stromal depth 670-1200 micron

dr = 5/10000;     %2.5
dz = 5/10000;       %5

ncr = rc/dr;
ner = he/dr+1;
nez = he/dz+1;
nsr = rs/dr+1;
nsz = hs/dz+1;

rgrid = [(1:ncr)/(ncr)*rc, rc+(0:ner-1)/(ner-1)*he, rc+he+(0:nsr-1)/(nsr-1)*rs];
zgrid = [(0:nez-1)/(nez-1)*he, he+(0:nsz-1)/(nsz-1)*hs];

cryptOpen = 1;
tnoflux = 10*60;

cryptFunc;

avgCf = zeros(size(t));
avgCe = zeros(size(t));
avgCet = zeros(size(t));
avgCs = zeros(size(t));
avgZs = zeros(size(t));

rrmat = rgrid'*ones(1,length(zgrid));
zzmat = ones(length(rgrid),1)*zgrid;

for i=1:length(t)
    crmat = reshape(Cout(i,1:length(rgrid)*length(zgrid)),length(rgrid),length(zgrid)) .* (rgrid'*ones(1,length(zgrid)));
    Mf = 0; Me = 0; Met = 0; Ms = 0; CRZs = 0;
    Rf = 0; Re = 0; Ret = 0; Rs = 0; CRs = 0;
    for j=1:length(rgrid)
       for k = 1:length(zgrid)
           gsize = 1;
           if(j==ncr || j==ncr+1 || j==ncr+ner || j==ncr+ner+1 || j==length(rgrid))         %j==1
               gsize = gsize/2;
           end
           if(k==1 || k==nez || k==nez+1 || k==length(zgrid))
               gsize = gsize/2;
           end
           if(j==ncr+ner && k==nez)
               gsize = 0.75;
           end
           
           if(j <= ncr)
               Mf = Mf + crmat(j,k)*gsize;
               Rf = Rf + rrmat(j,k)*gsize;
               CRZs = CRZs + crmat(j,k)*rrmat(j,k)*zzmat(j,k)*gsize;
               CRs = CRs + crmat(j,k)*rrmat(j,k)*gsize;
           elseif(j > ncr+ner && k > nez)
               Ms = Ms + crmat(j,k)*gsize;
               Rs = Rs + rrmat(j,k)*gsize;
           else
               Me = Me + crmat(j,k)*gsize;
               Re = Re + rrmat(j,k)*gsize;
               if(k <= nez)
                   Met = Met + crmat(j,k)*gsize;
                   Ret = Ret + rrmat(j,k)*gsize;
               end
           end
               
       end
    end
    avgCf(i) = Mf/Rf;
    avgCe(i) = Me/Re;
    avgCet(i) = Met/Ret;
    avgCs(i) = Ms/Rs;
    avgZs(i) = CRZs/CRs;
end


rgridp = [0 rgrid];

[X,Y] = meshgrid(zgrid,rgridp);

plot(0,0);
hold on;

cplot = reshape(Cout(10,:),[],nez+nsz);
cplot(:,nez) = (cplot(:,nez)+cplot(:,nez+1))/2;
cplot(ncr+ner,:) = (cplot(ncr+ner,:)+cplot(ncr+ner+1,:))/2;
cplot(:,nez+1) = []; X(:,nez+1) = []; Y(:,nez+1) = [];
cplot(ncr+ner+1,:)=[]; X(ncr+ner+1,:)=[]; Y(ncr+ner+1,:)=[];

cplot = [cplot(1,:); cplot];

surf(X,Y,cplot,'EdgeColor','none','FaceColor','interp');
surf(X,-fliplr(Y),cplot,'EdgeColor','none','FaceColor','interp');

view(90,90);

lineW = 2;

line([0 zgrid(end)],[rgridp(ncr+1) rgridp(ncr+1)],[1e10 1e10],'Color','w','LineWidth',lineW);
line([0 0],[rgridp(ncr+1) rgridp(end)],[1e10 1e10],'Color','w','LineWidth',lineW);
line([zgrid(nez) zgrid(end)],[rgridp(ncr+ner+1) rgridp(ncr+ner+1)],[1e10 1e10],'Color','w','LineWidth',lineW);
line([zgrid(nez) zgrid(nez)],[rgridp(ncr+ner+1) rgridp(end)],[1e10 1e10],'Color','w','LineWidth',lineW);

line([0 zgrid(end)],-[rgridp(ncr+1) rgridp(ncr+1)],[1e10 1e10],'Color','w','LineWidth',lineW);
line([0 0],-[rgridp(ncr+1) rgridp(end)],[1e10 1e10],'Color','w','LineWidth',lineW);
line([zgrid(nez) zgrid(end)],-[rgridp(ncr+ner+1) rgridp(ncr+ner+1)],[1e10 1e10],'Color','w','LineWidth',lineW);
line([zgrid(nez) zgrid(nez)],-[rgridp(ncr+ner+1) rgridp(end)],[1e10 1e10],'Color','w','LineWidth',lineW);

surf([0 0;-0.0025 -0.0025],[-rgridp(end) rgridp(end);-rgridp(end) rgridp(end)],[c0 c0;c0 c0],'EdgeColor','none');

colorbar;

%axis equal;
%axis([0 0.02 -0.01 0.01]);
axis([-0.0025 0.02 -0.01 0.01]);
axis square;

ylabel('Radial distance from center of crypt (cm)');
xlabel('Depth into tissue (cm)');
%title('Heat map of TFV concentration (ng/mL) around rectal crypt at 6 minutes post dosing');


%volumes of each compartment
vt = pi*(rc+he+rs)^2*(he+hs);
vc = pi*rc^2*(he+hs);
vs = pi*((rc+he+rs)^2-(rc+he)^2)*hs;
ve = vt-vc-vs;
vet = pi*((rc+he+rs)^2-(rc)^2)*he;


