function plotHeatMapRadial(r,t,tvec,C,bound)

%CONVERTING TO CARTESIAN FOR PLOTTING
 
 theta = linspace(0,2*pi,200);
 x = zeros(1+(length(r)-1)*length(theta),1);
 y = zeros(1+(length(r)-1)*length(theta),1);
 z = zeros(1+(length(r)-1)*length(theta),1);
 
 tumorX = r(bound)*cos(theta);
 tumorY = r(bound)*sin(theta);
 
 for j = 1:length(tvec)
     
     for i = 1:length(theta)
         x((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = r(2:end)*cos(theta(i));
         y((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = r(2:end)*sin(theta(i));
         z((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = C(2:end,tvec(j));
     end

     z(1) = C(1,tvec(j));

     data = [x,y,z];

     X = linspace(min(data(:,1)),max(data(:,1)),150);
     Y = linspace(min(data(:,2)),max(data(:,2)),150);

     [XX,YY]=meshgrid(X,Y);
     
     figure()

     F=scatteredInterpolant(data(:,1),data(:,2),data(:,3));
     contourf(XX,YY,F(XX,YY),100,'LineColor','none')
     xlim([-1.2,1.2])
     ylim([-1.2,1.2])
     hold on
     colormap bone
     colorbar
     plot(tumorX,tumorY,'w--','LineWidth',2)
     xlabel('X (cm)','FontSize',18,'FontWeight','Bold')
     ylabel('Y (cm)','FontSize',18,'FontWeight','Bold')
     title(sprintf('t = %.0f mins',floor(t(tvec(j))/60)),'FontSize',18,'FontWeight','Bold')
     XT = get(gca,'XTick');
     set(gca,'FontSize',16)
     YT = get(gca,'YTick');
     set(gca,'FontSize',16)

     
 end