clear
clc

%% Run simulations

%%031920s3avg.csv
% paramSpaceName = 'paramSpace2.csv';
% thickness = [1.2,1.12,1.13,1.16,1.01];
% h_S = mean(thickness)/10;
% outputFileName = '031920s3_Perf.csv';
% measurementFile = '031920s3avg.csv';

%%031920s3.csv PEAK
% paramSpaceName = 'paramSpace2.csv';
% thickness = [1.2,1.12,1.13,1.16,1.01];
% h_S = mean(thickness)/10;
% outputFileName = '031920s3_PeakPerf.csv';
% measurementFile = '031920s3.csv';

%%031920s2.csv
% paramSpaceName = 'paramSpace2.csv';
% thickness = [1.1,1.2,1.07,1.17,1.11,1.1];
% h_S = mean(thickness)/10;
% outputFileName = '031920s2_Perf.csv';
% measurementFile = '031920s2avg.csv';

%%031920s2.csv PEAK
% paramSpaceName = 'paramSpace2.csv';
% thickness = [1.1,1.2,1.07,1.17,1.11,1.1];
% h_S = mean(thickness)/10;
% outputFileName = '031920s2_PeakPerf.csv';
% measurementFile = '031920s2.csv';

% %%031920s1.csv
% paramSpaceName = 'paramSpace2.csv';
% thickness = [1.19,1.06,1.2,1.13,1.1,1.26];
% h_S = mean(thickness)/10;
% outputFileName = '031920s1_Perf.csv';
% measurementFile = '031920s1avg.csv';

%%031920s1.csv PEAK
% paramSpaceName = 'paramSpace2.csv';
% thickness = [1.19,1.06,1.2,1.13,1.1,1.26];
% h_S = mean(thickness)/10;
% outputFileName = '031920s1_PeakPerf.csv';
% measurementFile = '031920s1.csv';

%%101419.csv
% paramSpaceName = 'paramSpace2.csv';
% thickness = [1.75,1.7,1.53,1.56,1.74];
% h_S = mean(thickness)/10;
% outputFileName = '101419_Perf.csv';
% measurementFile = '101419avg.csv';

%%101419.csv PEAK
% paramSpaceName = 'paramSpace2.csv';
% thickness = [1.75,1.7,1.53,1.56,1.74];
% h_S = mean(thickness)/10;
% outputFileName = '101419_PeakPerf.csv';
% measurementFile = '101419.csv';

%TEST.csv
% paramSpaceName = 'paramSpace3.csv';
% thickness = [1.19,1.06,1.2,1.13,1.1,1.26];
% h_S = mean(thickness)/10;
% outputFileName = '031920s1_PeakPerfTEST.csv';
% measurementFile = '031920s1.csv';
% 
% twoCompartment_MFFit(paramSpaceName,h_S,measurementFile,outputFileName)

%% Calculate R2 and LSE

perfMatavg = [csvread('031920s3_Perf.csv'),csvread('031920s2_Perf.csv'),csvread('031920s1_Perf.csv'),csvread('101419_Perf.csv')];
perfMatPeak = [csvread('031920s3_PeakPerf.csv'),csvread('031920s2_PeakPerf.csv'),csvread('031920s1_PeakPerf.csv'),csvread('101419_PeakPerf.csv')];

paramSpace = csvread('paramSpace2.csv');

R2avg = [max(perfMatavg(:,1)),max(perfMatavg(:,3)),max(perfMatavg(:,5)),max(perfMatavg(:,7))]
LSEavg = [min(perfMatavg(:,2)),min(perfMatavg(:,4)),min(perfMatavg(:,6)),min(perfMatavg(:,8))]

R2peak = [max(perfMatPeak(:,1)),max(perfMatPeak(:,3)),max(perfMatPeak(:,5)),max(perfMatPeak(:,7))]
LSEpeak = [min(perfMatPeak(:,2)),min(perfMatPeak(:,4)),min(perfMatPeak(:,6)),min(perfMatPeak(:,8))]

R2LSEindavg = [find(perfMatavg(:,1)==max(perfMatavg(:,1)),1),find(perfMatavg(:,2)==min(perfMatavg(:,2)),1);...
            find(perfMatavg(:,3)==max(perfMatavg(:,3)),1),find(perfMatavg(:,4)==min(perfMatavg(:,4)),1);...
            find(perfMatavg(:,5)==max(perfMatavg(:,5)),1),find(perfMatavg(:,5)==min(perfMatavg(:,5)),1);...
            find(perfMatavg(:,7)==max(perfMatavg(:,7)),1),find(perfMatavg(:,7)==min(perfMatavg(:,7)),1)]
R2LSEindpeak = [find(perfMatPeak(:,1)==max(perfMatPeak(:,1)),1),find(perfMatPeak(:,2)==min(perfMatPeak(:,2)),1);...
            find(perfMatPeak(:,3)==max(perfMatPeak(:,3)),1),find(perfMatPeak(:,4)==min(perfMatPeak(:,4)),1);...
            find(perfMatPeak(:,5)==max(perfMatPeak(:,5)),1),find(perfMatPeak(:,5)==min(perfMatPeak(:,5)),1);...
            find(perfMatPeak(:,7)==max(perfMatPeak(:,7)),1),find(perfMatPeak(:,7)==min(perfMatPeak(:,7)),1)]
        
paramFittedavg = [paramSpace(R2LSEindavg(1,1),:);paramSpace(R2LSEindavg(2,1),:);paramSpace(R2LSEindavg(3,1),:);paramSpace(R2LSEindavg(4,1),:)];
paramFittedpeak = [paramSpace(R2LSEindpeak(1,1),:);paramSpace(R2LSEindpeak(2,1),:);paramSpace(R2LSEindpeak(3,1),:);paramSpace(R2LSEindpeak(4,1),:)];

thickness031920s3 = [1.2,1.12,1.13,1.16,1.01];
thickness031920s2 = [1.1,1.2,1.07,1.17,1.11,1.1];
thickness031920s1 = [1.19,1.06,1.2,1.13,1.1,1.26];
thickness101419 = [1.75,1.7,1.53,1.56,1.74];

thickness = [mean(thickness031920s3),mean(thickness031920s2),mean(thickness031920s1),mean(thickness101419)]/10;
measurementFileAVG = {'031920s3avg.csv','031920s2avg.csv','031920s1avg.csv','101419avg.csv'};
measurementFilePEAK = {'031920s3.csv','031920s2.csv','031920s1.csv','101419.csv'};

%% Plot fit for window averages
for i = 1:size(paramFittedavg,1)
  
[tm,tp,MsurfaceC,PsurfaceC] = twoCompartment_HS217evaluate(paramFittedavg(i,:),thickness(i),measurementFileAVG{i});

figure(i)
%plot(tp1/60,PsurfaceC1,'LineWidth',2)
%hold on
%plot(tp2/60,PsurfaceC2,'LineWidth',2)
%hold on
plot(tp/60,PsurfaceC,'LineWidth',2)
hold on
plot(tm,MsurfaceC,'kx')
%plot(mdl)
%legend('Model Predictions','Raman Measurements')
xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
ylabel('HSP90 Tissue Surface Concentration (\muM)','FontSize',18,'FontWeight','Bold')
%title('Ethanol concentration vs time on stromal surface','FontSize',18,'FontWeight','Bold')
XT = get(gca,'XTick');
set(gca,'FontSize',16)
YT = get(gca,'YTick');
set(gca,'FontSize',16)

end

%% Plot fit for window peaks
for i = 1:size(paramFittedpeak,1)
  
[tm,tp,MsurfaceC,PsurfaceC] = twoCompartment_HS217evaluate(paramFittedpeak(i,:),thickness(i),measurementFilePEAK{i});

figure(i+size(paramFittedavg,1))
%plot(tp1/60,PsurfaceC1,'LineWidth',2)
%hold on
%plot(tp2/60,PsurfaceC2,'LineWidth',2)
%hold on
plot(tp/60,PsurfaceC,'LineWidth',2)
hold on
plot(tm,MsurfaceC,'kx')
%plot(mdl)
%legend('Model Predictions','Raman Measurements')
xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
ylabel('HSP90 Tissue Surface Concentration (\muM)','FontSize',18,'FontWeight','Bold')
%title('Ethanol concentration vs time on stromal surface','FontSize',18,'FontWeight','Bold')
XT = get(gca,'XTick');
set(gca,'FontSize',16)
YT = get(gca,'YTick');
set(gca,'FontSize',16)

end

% perfMatTest = csvread('031920s1_PeakPerfTEST.csv');
% paramSpace = csvread('paramSpace3.csv');
% 
% R2test = max(perfMatTest(:,1));
% LSEtest = min(perfMatTest(:,2));
% R2LSEindtest = [find(perfMatTest(:,1)==max(perfMatTest(:,1)),1),find(perfMatTest(:,2)==min(perfMatTest(:,2)),1)];
% paramFittedtest = paramSpace(R2LSEindtest(1,1),:);
% thicknessTEST = mean([1.19,1.06,1.2,1.13,1.1,1.26])/10;
% measurementFileTEST = '031920s3avg.csv';
% [tm,tp,MsurfaceC,PsurfaceC] = twoCompartment_HS217evaluate(paramFittedtest,thicknessTEST,measurementFileTEST);
% figure()
% plot(tp/60,PsurfaceC,'LineWidth',2)
% hold on
% plot(tm,MsurfaceC,'kx')
% xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
% ylabel('HSP90 Tissue Surface Concentration (\muM)','FontSize',18,'FontWeight','Bold')
% %title('Ethanol concentration vs time on stromal surface','FontSize',18,'FontWeight','Bold')
% XT = get(gca,'XTick');
% set(gca,'FontSize',16)
% YT = get(gca,'YTick');
% set(gca,'FontSize',16)

