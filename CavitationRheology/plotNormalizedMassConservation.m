function plotNormalizedMassConservation(r, Q, Vol, t_OH, C_OH, C0_OH, a, varargin)

% r_axis = r(nr_cavity)';

% theta = linspace(0,2*pi,length(r))';

% mass_injected_oh_predicted = 4*pi*trapz((r'.^2).*C_OH)/(length(r)-1);
%plot(t_OH/60,mass_injected_oh_predicted)

%      n_a = find(a_vec(:,1) >= tvec1(j),1);
%      if (isempty(n_a))
%          n_a = a_vec(end,1);
%      end
%      bound = find(r>=a_vec(n_a,2),1);

% mass_inside_oh = zeros(numel(t_OH),1);
% mass_outside_oh = zeros(numel(t_OH),1);
% mass_injected_oh = zeros(numel(t_OH),1);
mass_injected_tot = (C0_OH*Vol).*ones(1,length(t_OH));

bound = find(r >= floor(a(1)), 1);

% theta_2 = linspace(0,2*pi,length(r(1:bound)))';
% theta_3 = linspace(0,2*pi,length(r(bound+1:end)))';

mass_inside_oh = 4*pi*trapz((r(1:bound)'.^2).*C_OH(1:bound,:),1)/(length(r(1:bound))-1);
mass_outside_oh = 4*pi*trapz((r(bound+1:end)'.^2).*C_OH(bound+1:end,:),1)/(length(r(1:bound))-1);
mass_injected_oh = mass_inside_oh+mass_outside_oh;

% for i = 1:numel(t_OH)
%     bound = find(r >= a(i));
%     mass_inside_oh(i) = 4*pi*trapz(r(1:bound),(r(1:bound)'.^2).*C_OH(1:bound,i),1);
%     mass_outside_oh(i) = 4*pi*trapz(r(bound+1:end),(r(bound+1:end)'.^2).*C_OH(bound+1:end,i),1);
%     mass_injected_oh(i) = mass_inside_oh(i)+mass_outside_oh(i);
% end


figure()
plot(t_OH/60,mass_inside_oh,'LineWidth',2)
hold on
plot(t_OH/60,mass_outside_oh,'r','LineWidth',2)
plot(t_OH/60,mass_injected_oh,'k--','LineWidth',2)
plot(t_OH/60,mass_injected_tot,'k','LineWidth',2)

ylabel('Ethanol Mass (mg)','FontSize',18,'FontWeight','Bold')
xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
title('Total Mass of Ethanol Inside and Outside Tumor','FontSize',18,'FontWeight','Bold')
XT = get(gca,'XTick');
set(gca,'FontSize',16)
YT = get(gca,'YTick');
set(gca,'FontSize',16)
legend('Inside tumor', 'Outside Tumor', 'Inside + outside', 'Total injected')
%plot(t2/60,TotMassPredic,'g--','LineWidth',2)
%ylim([0,0.1])



end