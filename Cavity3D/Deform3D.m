clear
clc

%% Complete Domain Creation (here, cuboid)

[xd, yd, zd] = meshgrid(-4:0.2:4);
Pdomain = [xd(:) yd(:), zd(:)];

%% Initial Cavity Creation

[xc,yc,zc] = sphere(24);
zc(:) = zc(:)+1;

% N = 20;
% thetavec = linspace(0,pi,N);
% phivec = linspace(0,2*pi,2*N);
% [th, ph] = meshgrid(thetavec,phivec);
% R = ones(size(th)); % should be your R(theta,phi) surface in general
%   
% xc = R.*sin(th).*cos(ph);
% yc = R.*sin(th).*sin(ph);
% zc = R.*cos(th)+2;

Pcavity = [xc(:) yc(:) zc(:)];
%Pcavity = unique(Pcavity,'rows');
cavityshp = alphaShape(Pcavity,2); %%creating alphashape of cavity to do mesh subtraction

%% Tumor Domain (removing initial spherical cavity)

xt = xd; yt = yd; zt = zd;
tf = inShape(cavityshp,xd(:),yd(:),zd(:));
xt(tf) = []; yt(tf) = []; zt(tf) = [];
Ptissue = [xt(:) yt(:), zt(:)];

%% Creating overall meshes

tic
%trDomain = boundary(Pdomain(:,1),Pdomain(:,2),Pdomain(:,3));
trCavity = convhull(Pcavity(:,1),Pcavity(:,2),Pcavity(:,3));
trTissue = boundary(Ptissue(:,1),Ptissue(:,2),Ptissue(:,3));
toc
%% Plotting tissue and cavity

plotTumorCavity(trTissue,trCavity,Pdomain,Ptissue,Pcavity,tf)

%% Multiplying by transformation matrix

lambda1 = 0.5;
lambda2 = 0.5;
lambda3 = 1.5;

F = [lambda1, 0, 0; 0, lambda2, 0; 0, 0, lambda3];
Y = Pcavity*F;
cavityshp = alphaShape(Y,2);

tf = inShape(cavityshp,xd(:),yd(:),zd(:));

xt = xd; yt = yd; zt = zd;
xt(tf) = []; yt(tf) = []; zt(tf) = [];
Ptissue = [xt(:) yt(:), zt(:)];
Pcavity = Y;

trTissue = boundary(Ptissue(:,1),Ptissue(:,2),Ptissue(:,3));
trCavity = convhull(Pcavity(:,1),Pcavity(:,2),Pcavity(:,3));

plotTumorCavity(trTissue,trCavity,Pdomain,Ptissue,Pcavity,tf)
