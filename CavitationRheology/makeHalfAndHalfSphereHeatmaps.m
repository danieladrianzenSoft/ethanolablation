function makeHalfAndHalfSphereHeatmaps(t2_top,tvec_top,tvec_bottom,r,C_top,C_bottom,a_vec)

 %CONVERTING TO CARTESIAN FOR PLOTTING
 
 theta_top = linspace(pi/2,3*pi/2,100);
 x_top = zeros(1+(length(r)-1)*length(theta_top),1);
 y_top = zeros(1+(length(r)-1)*length(theta_top),1);
 z_top = zeros(1+(length(r)-1)*length(theta_top),1);
 
 theta_bottom = linspace(3*pi/2,5*pi/2,100);
 x_bottom = zeros(1+(length(r)-1)*length(theta_bottom),1);
 y_bottom = zeros(1+(length(r)-1)*length(theta_bottom),1);
 z_bottom = zeros(1+(length(r)-1)*length(theta_bottom),1);
 
  
 for j = 1:length(tvec_top)
     
     n_a = find(a_vec(:,1) >= tvec_top(j),1);
     if (isempty(n_a))
         n_a = a_vec(end,1);
     end
     bound = find(r>=a_vec(n_a,2),1);
     
     tumorX_top = r(bound)*cos(theta_top);
     tumorY_top = r(bound)*sin(theta_top);
     tumorX_bottom = r(bound)*cos(theta_bottom);
     tumorY_bottom = r(bound)*sin(theta_bottom);
     
     for i = 1:length(theta_top)
         x_top((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = r(2:end)*cos(theta_top(i));
         y_top((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = r(2:end)*sin(theta_top(i));
         z_top((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = C_top(2:end,tvec_top(j));
     end
     for i = 1:length(theta_bottom)
         x_bottom((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = r(2:end)*cos(theta_bottom(i));
         y_bottom((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = r(2:end)*sin(theta_bottom(i));
         z_bottom((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = C_bottom(2:end,tvec_bottom(j));
     end
     

     z_top(1) = C_top(1,tvec_top(j));
     z_bottom(1) = C_top(1,tvec_bottom(j));


     data_top = [x_top,y_top,z_top];
     data_bottom = [x_bottom,y_bottom,z_bottom];


     X_top = linspace(min(data_top(:,1)),max(data_top(:,1)),300);
     Y_top = linspace(min(data_top(:,2)),max(data_top(:,2)),300);
     X_bottom = linspace(min(data_bottom(:,1)),max(data_bottom(:,1)),300);
     Y_bottom = linspace(min(data_bottom(:,2)),max(data_bottom(:,2)),300);

     [XX_top,YY_top]=meshgrid(X_top,Y_top);
     [XX_bottom,YY_bottom]=meshgrid(X_bottom,Y_bottom);

     
     figure()

     F_top=scatteredInterpolant(data_top(:,1),data_top(:,2),data_top(:,3));
     contourf(XX_top,YY_top,F_top(XX_top,YY_top),250,'LineColor','none')
     xlim([-1.2,1.2])
     ylim([-1.2,1.2])
     hold on
     colormap parula
     colorbar
     caxis([0, 1]);

     
     F_bottom=scatteredInterpolant(data_bottom(:,1),data_bottom(:,2),data_bottom(:,3));
     contourf(XX_bottom,YY_bottom,F_bottom(XX_bottom,YY_bottom),250,'LineColor','none')
     
     plot(tumorX_bottom,tumorY_bottom,'w--','LineWidth',2)
     plot(tumorX_top,tumorY_top,'w--','LineWidth',2)
     plot(zeros(1,length(Y_top)+length(Y_bottom)),[Y_top Y_bottom],'w','LineWidth',4)

     set(gca,'FontSize',30);
     xlabel('X (cm)','FontSize',30)
     ylabel('Y (cm)','FontSize',30)
     title(sprintf('t = %.0f mins',floor(t2_top(tvec_top(j))/60)),'FontSize',30)

     XT = get(gca,'XTick');
     set(gca,'FontSize',30)
     YT = get(gca,'YTick');
     set(gca,'FontSize',30)
     set(gcf, 'Color', [1 1 1], 'Position', [333 193 737 564])
     
 end
end

% function makeHalfAndHalfSphereHeatmaps(t2_top,tvec_top,tvec_bottom,r,C_top,C_bottom,a_vec)
%  %CONVERTING TO CARTESIAN FOR PLOTTING
%  
%  theta_top = linspace(0,2*pi,200);
%  theta_bottom = linspace(pi,2*pi,100);
% 
%  x_top = zeros(1+(length(r)-1)*length(theta_top),1);
%  y_top = zeros(1+(length(r)-1)*length(theta_top),1);
%  z_top = zeros(1+(length(r)-1)*floor(length(theta_top)/2),1);
%  z_bottom = zeros(1+(length(r)-1)*floor(length(theta_bottom)/2),1);
% 
%     for j = 1:length(tvec_top)
%      
%          n_a = find(a_vec(:,1) >= tvec_top(j),1);
%          if (isempty(n_a))
%              n_a = a_vec(end,1);
%          end
%          bound = find(r>=a_vec(n_a,2),1);
% 
%          tumorX = r(bound)*cos(theta_top);
%          tumorY = r(bound)*sin(theta_top);
% 
%          for i = 1:length(theta_top)
%              x_top((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = r(2:end)*cos(theta_top(i));
%              y_top((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = r(2:end)*sin(theta_top(i));
%              z_top((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = C_top(2:end,tvec_top(j));
% 
% %              if (theta(i) <= pi)
% %                 z_top((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = C_top(2:end,tvec_top(j));
% %              elseif (theta(i) > pi && theta(i) <= 2*pi)
% %                 z_bottom((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = C_bottom(2:end,tvec_bottom(j));
% %              end
% 
%          end
%          
% 
%          z_top(1) = C_top(1,tvec_top(j));
%          z_bottom(1) = C_bottom(1,tvec_bottom(j));
%          
%          data_top = [x_top,y_top,z_top];
% 
%          X_top = linspace(min(data_top(:,1)),max(data_top(:,1)),300);
%          Y_top = linspace(min(data_top(:,2)),max(data_top(:,2)),300);
% 
%          %data_top = [x,y,z_top];
%          %data_bottom = [x,y,z_bottom];
% 
%          % X_top = linspace(min(x_top),max(x_top),300);
%          % Y_top = linspace(min(y_top),max(y_top),300);
% %          X_bottom = linspace(min(x),max(x),300);
% %          Y_bottom = linspace(min(y),0,300);
% 
%          [XX_top,YY_top]=meshgrid(X_top,Y_top);
% %          [XX_bottom,YY_bottom]=meshgrid(X_bottom,Y_bottom);
% 
% 
%          figure()
% 
%          F_top=scatteredInterpolant(data_top(:,1),data_top(:,2),data_top(:,3));
%          contourf(XX_top,YY_top,F_top(XX_top,YY_top),250,'LineColor','none')
%          xlim([-1.2,1.2])
%          ylim([-1.2,1.2])
% %          hold on
% %          F_bottom = scatteredInterpolant(X_bottom,Y_bottom,Z_bottom);
% %          contourf(XX_bottom,YY_bottom,F_bottom(XX_bottom,YY_bottom),250,'LineColor','none')
% 
%          colormap parula
%          colorbar
%          caxis([0, 1]);
%          plot(tumorX,tumorY,'w--','LineWidth',2)
%          xlabel('X (cm)','FontSize',18,'FontWeight','Bold')
%          ylabel('Y (cm)','FontSize',18,'FontWeight','Bold')
%          title(sprintf('t = %.0f mins',floor(t2_top(tvec_top(j))/60)),'FontSize',18,'FontWeight','Bold')
%          XT = get(gca,'XTick');
%          set(gca,'FontSize',16)
%          YT = get(gca,'YTick');
%          set(gca,'FontSize',16)
%      
%     end
% end
% 
