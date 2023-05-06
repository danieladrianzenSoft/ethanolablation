function plotCavityStrainRadius(tt, lambda, cavity_radius, parameter_name, parameter)

arrayLegend = cell(1,length(parameter));
for i = 1:length(parameter)
    % units = '';
    if (parameter_name == 'k')
        arrayLegend{i} = sprintf('k = %.0d cm^{4} dyn^{-1} s^{-1}',parameter(i));
    elseif (parameter_name == 'o')
        arrayLegend{i} = sprintf('%c = %.3f %s',char(937),parameter(i),'%');
    else
        arrayLegend{i} = sprintf('E = %.0e dyn cm^{-2}',parameter(i));
    end
    
    
end

% figure()
% cm = colormap('colorcube');
% cm = cm([15,40,55],:);
% for i = 1:length(parameter)
% plot(tt(:,i)/60,lambda(:,i),'LineWidth',2,'Color',cm(i,:))
% hold on
% end

% legend(arrayLegend)
% xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
% ylabel('\lambda','FontSize',18,'FontWeight','Bold')
% set(gca,'FontSize',14)

figure()

%cm = colormap('jet');
%cm = cm([15,100,155,205],:);
cm = [[0, 0, 0]; [1, 0, 0]; [0, 0, 1]; [1, 0, 1];];

for i = 1:length(parameter)
plot(tt(:,i)/60,cavity_radius(:,i),'LineWidth',2,'Color',cm(i,:))
hold on
end

legend(arrayLegend)
xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
ylabel('Cavity radius (cm)','FontSize',18,'FontWeight','Bold')
set(gca,'FontSize',14)
