function plotPressures(tt, lambda, Pc, Pup, Pcrit, parameter_name, parameter)

figure()
%cm = colormap('winter');
%cm = cm([15,100,160],:);
cm = [[0, 0, 0]; [1, 0, 0]; [0, 0, 1]; [1, 0, 1];];

plot(nan,nan, '-', 'LineWidth', 2, 'Color', '[0.4, 0.4, 0.4]')
hold on
plot(nan,nan, '--', 'LineWidth', 2, 'Color', '[0.4, 0.4, 0.4]')
plot(nan,nan, ':', 'LineWidth', 2, 'Color', '[0.4, 0.4, 0.4]')
plot(nan,nan, '-w', 'LineWidth', 2)

arrayLegend = cell(1,length(parameter));

for i = 1:length(parameter)
plot(lambda(:,i),Pc(:,i)*10^(-4),'Color', cm(i,:), 'LineWidth',2);
plot(lambda(:,i),Pcrit(:,i)*10^(-4),'Color', cm(i,:), 'LineWidth',2,'LineStyle','--');
plot(lambda(:,i),Pup(:,i)*10^(-4), 'Color', cm(i,:), 'LineWidth', 2,'LineStyle',':');
    if (parameter_name == 'k')
        arrayLegend{i} = sprintf('k = %.0d cm^{4} dyn^{-1} s^{-1}',parameter(i));
    elseif (parameter_name == 'o')
        arrayLegend{i} = sprintf('%c = %.3f %s',char(937),parameter(i),'%');
    else
        arrayLegend{i} = sprintf('E = %.0e dyn cm^{-2}',parameter(i));
    end
    %text(lambda(1,1), Pc(1,1), arrayLegend)
  
end
legendarray = {'Pc','Pcrit','Pup',' '};
for i = 1:length(parameter)
    legendarray{end+1} = arrayLegend{i};
    legendarray{end+1} = '';
    legendarray{end+1} = '';
end

%legend('Pc','Pcrit','Pup',' ',arrayLegend{1},'','',arrayLegend{2},'','',arrayLegend{3},'','')
legend(legendarray)


 
xlabel('\lambda','FontSize',18,'FontWeight','Bold')
%ylabel('P (dyn/cm^{2})','FontSize',18,'FontWeight','Bold')
ylabel('P (kPa)','FontSize',18,'FontWeight','Bold')
%legend('P_{c}','P_{crit}','P_{up}');
set(gca,'FontSize',14)

figure()
for i = 1 : length(parameter)
    plot(tt(:,i)/60, Pc(:,i)*10^(-4), 'Color', cm(i,:), 'LineWidth', 2)
    hold on
end
xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
%ylabel('P (dyn/cm^{2})','FontSize',18,'FontWeight','Bold')
ylabel('P (kPa)','FontSize',18,'FontWeight','Bold')
set(gca,'FontSize',14)

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
legend(arrayLegend);


%text(arrayLegend)

%legendObject = get(gcf).Children(1)

%legendObject.Color = [1, 0, 0];

end