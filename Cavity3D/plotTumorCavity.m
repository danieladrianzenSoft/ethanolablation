function plotTumorCavity(trTissue,trCavity,Pdomain,Ptissue,Pcavity,tf)

h.fig = figure('Color','w') ;
h.ax = axes;
%colorTissue = [0.95,0.8,0.8];
colorTissue = 'none';
meshColor = [0.8,0.8,0.8];
meshAlpha = 0.2;
%meshColor = 'none';
h.patch1 = trisurf(trTissue,Ptissue(:,1),Ptissue(:,2),Ptissue(:,3),'Facecolor',colorTissue,'FaceAlpha',0.2,'EdgeColor',meshColor,'EdgeAlpha',meshAlpha);
h.alpha = 0.5;

hold on
%colorCavity = [0.5,0.9,0.95];
colorCavity = 'none';
meshColor = [0.6,0.75,0.9];
meshAlpha = 0.5;
h.patch1 = trisurf(trCavity,Pcavity(:,1),Pcavity(:,2),Pcavity(:,3),'Facecolor',colorCavity,'FaceAlpha',0.8,'EdgeColor',meshColor,'EdgeAlpha',meshAlpha);
%plot3(Pdomain(tf,1),Pdomain(tf,2),Pdomain(tf,3),'b.')
%plot3(xd(~tf),yd(~tf),zd(~tf),'r.')

xlim(h.ax,[min(Pdomain(:,1)) max(Pdomain(:,1))])
ylim(h.ax,[min(Pdomain(:,2)) max(Pdomain(:,2))])
zlim(h.ax,[min(Pdomain(:,3)) max(Pdomain(:,3))])

axis equal

end