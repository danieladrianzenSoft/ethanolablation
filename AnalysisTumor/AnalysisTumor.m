clear
clc

%Ethanol + Ethyl Cellulose
[t2_1,PsurfaceC_1,tm_1,MsurfaceC_1,R2_1,D_S1] = twoCompartment_Ethanol021519S1();
%[t2_1,PsurfaceC_1,tm_1,MsurfaceC_1,R2_1,D_S1] = twoCompartment_Ethanol011619();
%[t2_1,PsurfaceC_1,tm_1,MsurfaceC_1,R2_1,D_S1] = twoCompartment_Ethanol021219S1();

%Ethanol alone
[t2_2,PsurfaceC_2,tm_2,MsurfaceC_2,R2_2,D_S2] = twoCompartment_Ethanol112018();
%[t2_2,PsurfaceC_2,tm_2,MsurfaceC_2,R2_2,D_S2] = twoCompartment_Ethanol022019S1();

%PsurfaceC_2 = PsurfaceC_2/(max(PsurfaceC_1));
%MsurfaceC_2 = MsurfaceC_2/(max(MsurfaceC_1));

figure()
p1 = plot(t2_1/60,PsurfaceC_1,'color',[176,36,24]/255,'LineWidth',2);
hold on
plot(tm_1,MsurfaceC_1,'x','color',[176,36,24]/255)
p2 = plot(t2_2/60,PsurfaceC_2,'LineWidth',2,'color',[5,49,245]/255);
hold on
plot(tm_2,MsurfaceC_2,'x','color',[5,49,245]/255)
%plot(TTemp,ConcTemp,'ro')
%legend([p1,p2],'6% EECF','Ethanol Alone','LineWidth',0)
xlabel('Time (mins)','FontSize',20,'FontWeight','Bold')
ylabel('Ethanol Concentration (normalized)','FontSize',20,'FontWeight','Bold')
%title('Ethanol concentration vs time on tumor surface','FontSize',18,'FontWeight','Bold')
XT = get(gca,'XTick');
set(gca,'FontSize',16,'LineWidth',2)
YT = get(gca,'YTick');
set(gca,'FontSize',16,'LineWidth',2)
box off
legend([p1,p2],'6% EtOH-EC','Ethanol Alone')
legend boxoff



D_S1
R2_1
D_S2
R2_2


% Dvec1 = 1:5;
% Dvec2 = 6:10;
% Dvec3 = 11:15;
% Dvec = [Dvec1',Dvec2',Dvec3'];
% labels = ['1';'2';'3'];

%DvecEtOH = [4.5,5.4,8.1,1.2,2.2,6.3,7.0,4.1]';
DvecEtOH = [4.5,5.4,8.1,6.3,7.0,4.1]';
DvecEtEC = [3.6,2.3,2.6,1.8,3.4,5.0,4.4]';
%Dvec = [DvecEtOH,DvecEtEC];
labels = {'EtOH';'6% EECF'};


% figure()
% CategoricalScatterplot(Dvec,'labels',cellstr(labels));
% ylabel('Diffusion coefficient (cm^{2}/s)','FontWeight','Bold','FontSize',14)
% %boxplot(Dvec,cellstr(labels));


%data = bsxfun(@times, rand(5,3), [50 150 100]);                     % Create Data
%DvecEtOH = [4.5,5.4,8.1,5.2,6.1,7.0,4.1]';
%DvecEtEC = [3.6,2.3,2.6,1.8,3.4,5.0,4.4]';
%data = [DvecEtOH,DvecEtEC];

dmean = [mean(DvecEtOH),mean(DvecEtEC)];                                             % Mean
%dci  = std(data)*tinv(0.975,size(data,1)-1);                        % Confidence Intervals
dci = [std(DvecEtOH)/sqrt(length(DvecEtOH)),std(DvecEtEC)/sqrt(length(DvecEtOH))];         % Standard Error
%dci = std(data);                                                    % Standard Deviation
%dci = prctile(data,0.1);
xt = [0.2,0.4];                                                         % X-Ticks
%xtd = repmat(xt, size(data,1), 1);                                  % X-Ticks For Data
xtd = repmat(xt, max([size(DvecEtOH,1),size(DvecEtEC,1)]), 1);                                  % X-Ticks For Data
sb = [xt'-ones(2,1)*0.05,  xt'+ones(2,1)*0.05]; % Short Bar X
lb = [xt'-ones(2,1)*0.08,  xt'+ones(2,1)*0.08]; % Long Bar X

verty = [linspace(dmean(1)-dci(1),dmean(1)+dci(1),10);linspace(dmean(2)-dci(2),dmean(2)+dci(2),10)];
vertx = [xt(1)*ones(10,1),xt(2)*ones(10,1)];

figure()
plot(xt(1), DvecEtOH, 'ks','MarkerFaceColor','k')
hold on
plot(xt(2), DvecEtEC, 'ks','MarkerFaceColor','k')
for k1 = 1:2
    %plot(lb(k1,:), [1 1]*dmean(k1), sb(k1,:), [1 1]*(dmean(k1)-dci(k1)), sb(k1,:), [1 1]*(dmean(k1)+dci(k1)), '-','color','k','linewidth',2)
    plot(sb(k1,:), [1 1]*(dmean(k1)-dci(k1)), sb(k1,:), [1 1]*(dmean(k1)+dci(k1)), '-','color','k','linewidth',2)
end
for k1 = 1:2
    plot(lb(k1,:), [1 1]*dmean(k1),'-','color','k','linewidth',1.5)
end
for k1 = 1:2
    plot(vertx(:,k1),verty(k1,:),'-','color','k','linewidth',1)
end
hold off
box off
%set(gca,'linewidth',6)
set(gca, 'XTick', xt, 'XTickLabel', {'EtOH','6% EtOH-EC'},'FontSize',20,'FontWeight','Bold','LineWidth',2)
set(gca, 'Color',[216,232,207]/255)
%xlabel('Group')
ylabel('D_{eff} (x10^{-6} cm^{2}/s)','FontSize',20,'FontWeight','Bold')

%% 10 and 90 percentile

% dmedian = median(data);                                                 % Mean
% %dci  = std(data)*tinv(0.975,size(data,1)-1);                        % Confidence Intervals
% %dci = std(data)/sqrt(length(data));                                 % Standard Error
% %dci = std(data);                                                    % Standard Deviation
% dci = prctile(data,[10,90],1)
% xt = [0.2,0.4];                                                         % X-Ticks
% xtd = repmat(xt, size(data,1), 1);                                  % X-Ticks For Data
% sb = [xt'-ones(size(data,2),1)*0.05,  xt'+ones(size(data,2),1)*0.05]; % Short Bar X
% lb = [xt'-ones(size(data,2),1)*0.08,  xt'+ones(size(data,2),1)*0.08]; % Long Bar X
% 
% %verty = [linspace(dmean(1)-dci(1),dmean(1)+dci(1),10);linspace(dmean(2)-dci(2),dmean(2)+dci(2),10)];
% vertx = [xt(1)*ones(10,1),xt(2)*ones(10,1)];
% 
% figure()
% plot(xt, data, 'ks','MarkerFaceColor','k')
% hold on
% for k1 = 1:size(data,2)
%     %plot(lb(k1,:), [1 1]*dmean(k1), sb(k1,:), [1 1]*(dmean(k1)-dci(k1)), sb(k1,:), [1 1]*(dmean(k1)+dci(k1)), '-','color','k','linewidth',2)
%     plot(sb(k1,:), [1 1]*dci(1,k1), sb(k1,:), [1 1]*dci(2,k1), '-','color','k','linewidth',2)
% end
% for k1 = 1:size(data,2)
%     plot(lb(k1,:), [1 1]*dmedian(k1),'-','color','k','linewidth',1.5)
% end
% % for k1 = 1:size(data,2)
% %     plot(vertx(:,k1),verty(k1,:),'-','color','k','linewidth',1)
% % end
% hold off
% set(gca, 'XTick', xt, 'XTickLabel', {'EtOH','6%EECF'},'FontSize',20,'FontWeight','Bold')
% %xlabel('Group')
% ylabel('D_{eff} (x10^{-6} cm^{2}/s)','FontSize',20,'FontWeight','Bold')
%%


%[h,p] = ttest(DvecEtOH,DvecEtEC,'alpha',0.05)
%[h,p] = ttest(data(:,1),data(:,2),'tail','right')
% figure()
% 
% mat = [[DvecEtOH;DvecEtEC],[ones(length(DvecEtOH'),1);2*ones(length(DvecEtEC'),1)]];
% X = mat(:,1);
% G = mat(:,2);
% 
% [p,table,stats,terms] = anovan(X, G, 'model','interaction', 'display','on');
% [c,m,h,nms] = multcompare(stats,'CType','tukey-kramer');


%[p,tbl,stats] = anovan();


