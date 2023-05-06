
clear;
clc;
% Load Depth Scan spectra 

PDMS = load('PDMSNoSub.mat');
LiverBlobInjectedEthanol = load('LiverBlobEthanolInjectedNoSub.mat');
LiverBlob = load('liverblobNoSub.mat');
Ethanol = load('ethanolNoSub.mat');


wavenumPDMS = PDMS.data.wavenumcrop;
spectraPDMS = PDMS.data.speccrop

% overtime = load('TissueOH.mat');
% wavenumOT = overtime.data.wavenumcrop;
% spectraOT = overtime.data.speccrop;

for k = 1:size(PDMS,2)
    %smoothspecPDMS(:, k) = smooth(spectraPDMS(:, k), sdex, 'sgolay', 2);
    smoothspecPDMS(:, k) = spectraPDMS(:,k);
    normspecPDMS(:, k) = spectraPDMS(:, k)/mean(spectraPDMS(:,k));
end

figure()
plot(wavenumPDMS, normspecPDMS)
ylabel('PDMS','FontSize',16,'FontWeight','Bold')


wavenumEthanol = Ethanol.data.wavenumcrop;
spectraEthanol = Ethanol.data.speccrop

% overtime = load('TissueOH.mat');
% wavenumOT = overtime.data.wavenumcrop;
% spectraOT = overtime.data.speccrop;

for k = 1:size(Ethanol,2)
    %smoothspecPDMS(:, k) = smooth(spectraPDMS(:, k), sdex, 'sgolay', 2);
    smoothspecEthanol(:, k) = spectraEthanol(:,k);
    normspecEthanol(:, k) = spectraEthanol(:, k)/mean(spectraEthanol(:,k));
end

normspecEthanol = fshift(normspecEthanol, 106);

figure()
plot(wavenumEthanol, normspecEthanol)
ylabel('Ethanol','FontSize',16,'FontWeight','Bold')


wavenumLiverBlob = LiverBlob.data.wavenumcrop;
spectraLiverBlob = LiverBlob.data.speccrop

% overtime = load('TissueOH.mat');
% wavenumOT = overtime.data.wavenumcrop;
% spectraOT = overtime.data.speccrop;

for k = 1:size(LiverBlob,2)
    %smoothspecPDMS(:, k) = smooth(spectraPDMS(:, k), sdex, 'sgolay', 2);
    smoothspecLiverBlob(:, k) = spectraLiverBlob(:,k);
    normspecLiverBlob(:, k) = spectraLiverBlob(:, k)/mean(spectraLiverBlob(:,k));
end

%     normspecEthanol(:, k) = fshift(normspecEthanol(:,k),103);


figure()
plot(wavenumLiverBlob, normspecLiverBlob)
ylabel('Liver','FontSize',16,'FontWeight','Bold')


wavenumLiverBlobInjectedEthanol = LiverBlobInjectedEthanol.data.wavenumcrop;
spectraLiverBlobInjectedEthanol = LiverBlobInjectedEthanol.data.speccrop

% overtime = load('TissueOH.mat');
% wavenumOT = overtime.data.wavenumcrop;
% spectraOT = overtime.data.speccrop;

for k = 1:size(LiverBlobInjectedEthanol,2)
    %smoothspecPDMS(:, k) = smooth(spectraPDMS(:, k), sdex, 'sgolay', 2);
    smoothspecLiverBlobInjectedEthanol(:, k) = spectraLiverBlobInjectedEthanol(:,k);
    normspecLiverBlobInjectedEthanol(:, k) = spectraLiverBlobInjectedEthanol(:, k)/mean(spectraLiverBlobInjectedEthanol(:,k));
end

figure()
plot(wavenumLiverBlobInjectedEthanol, normspecLiverBlobInjectedEthanol)
ylabel('Liver And Injected Ethanol','FontSize',16,'FontWeight','Bold')


%% Establish Basis array
% Basis Array for positive fitting 
% C(:, 1:272) = normRB(:,:);
% C(:, 273:282) = normgel(:,:);
% C(:, 283:329) = normL(:,:);
% C(:, 330) = normTFV(:,:);
% 
[RR,CC] = size(normspecLiverBlob);

Livervar = [1:CC];
% Bkgdvar = [CC+1];
%Lvar = [188:250];
OHvar = [CC+2];

%size(RBvar)
%size(normRB)

%[tempsizeR,tempsizeC] = size(normRB);

C = zeros(length(normspecLiverBlob),CC+2);
%size(normRB(speccrop,:))
%normOH
C(:, Livervar) = normspecLiverBlob(:,:);
% C(:, OHvar) = normBkgd(speccrop,:);
% C(:, Lvar) = normL(speccrop,:);  
C(:, OHvar) = normspecEthanol(:,:);

% figure()
% plot(speccrop,C(:, OHvar))

%% Chemometric Fitting
% Find the coeeficients of the matrix computation [Fit] = [Basis][Coeff]

% Scan 1
for k = 1:size(normspecLiverBlobInjectedEthanol,2);
    coeff1(1:size(C,2), k) = lsqnonneg(C,normspecLiverBlobInjectedEthanol(:,k));
end

%% Fitting and Determining Relative Contributions - Positive Fit
% Find the fit spectra and smooth
% Find the tissue Contribution and the TFV contribution

% Scan 1 - Tissue 
for k = 1:size(normspecLiverBlobInjectedEthanol,2);
%for k = 1:5;
    fit1(:,k) = C(:,:)*coeff1(:,k);
    %fit1(:,k) = smooth(fit1(:,k), 20, 'sgolay',2);
    tissueContr1(:,k) = C(:, Livervar)*coeff1(Livervar,k);
    OHbasisContr1(:,k) = C(:, OHvar)*coeff1(OHvar, k);
    OHContr1(:,k) = trapz(OHbasisContr1(:,k))./(trapz(tissueContr1(:,k)));
    Tissue1(:,k) = trapz(tissueContr1(:,k))./trapz(fit1(:,k));

end

%%

%%

% CM = lines(5);
% figure(); clf
% hold on
% plot(x1(:,[1:end])+shift(1), GContr1(:,[1:end]), '-o', 'color',CM(1,:), 'LineWidth', 3)
% plot(x23(:,[1:end])+shift(2), GContr2(:,[1:end]), '-o', 'color',CM(2,:), 'LineWidth', 3)
% plot(x23(:,[1:end])+shift(3), GContr3(:,[1:end]), '-o', 'color',CM(3,:), 'LineWidth', 3)
% xlabel('x Depth (\mum)','FontSize',12,'FontWeight','Bold'); ylabel('Relative Concentration (a.u.)','FontSize',12,'FontWeight','Bold')
% title('Raman Relative Concentrations of Gel in Porcine Vaginal Tissue, 2hr incubation','FontSize',14,'FontWeight','Bold')
% legend('Scan1', 'Scan2', 'Scan3')

%%

% figure(2); clf
% for i = 1:size(Scan3(:,:), 2);
%     plot(wavenum, Scan3(:,i))
%     pause
%     
%     hold on
% end

%% Area Plots of chemometric fit
% Figure for representative ES Zone spectra
%for i =1:size(Scan1crop,2)
%for i =length(x1)-5:length(x1)

for i = 1:size(normspecLiverBlobInjectedEthanol,2)

figure(); clf  
a=area(wavenumLiverBlob, tissueContr1(:,i),'FaceColor',[0 0 1]);
%a.FaceColor = [0 0 1];
hold on

%e.FaceColor = [0 1 1];

b = area(wavenumEthanol, OHbasisContr1(:,i),'FaceColor',[0 1 0]);
%b.FaceColor = [0 1 0];
% 
% c=area(wavenumcrop, TFVbasisContr1(:,i),'FaceColor',[1 0 0]);
% %c.FaceColor = [1 0 0];


plot(wavenumLiverBlobInjectedEthanol, normspecLiverBlobInjectedEthanol(:,i), 'Color',[0, 130/255, 255/255],'LineWidth', 1.5)
plot(wavenumLiverBlobInjectedEthanol, fit1(:,i),'Color',[255/255, 139/255, 0], 'LineWidth', 1.5)
plot(wavenumLiverBlobInjectedEthanol, normspecLiverBlobInjectedEthanol(:,i)-fit1(:,i),'Color',[255/255, 196/255, 0])
legend('Tissue', 'Ethanol','Raw Data', 'Fit', 'Residuals')
xlabel('Raman Shift (cm^{-1})','FontSize',18,'FontWeight','Bold'); ylabel(' Raman Intensity (a.u.)','FontSize',18,'FontWeight','Bold')
title(sprintf('Representative Spectra for Ethanol-injected Tissue'),'FontSize',18,'FontWeight','Bold')
XTickL = get(gca,'XTick');
set(gca,'FontSize',14);
YTickL = get(gca,'YTick');
set(gca,'FontSize',14);
axis('tight')

Residual = mean(abs(normspecLiverBlobInjectedEthanol(:,i)-fit1(:,i))./normspecLiverBlobInjectedEthanol(:,i))*100;

end

% FinalConc = [x1',OHContr1']
% dlmwrite(sprintf('DataVsTime.csv'),FinalConc)

% %% Convert relative concentrations to "true" concentration
% 
% %y = @(x) 10.20 .* x + 0.01158;    % calibration curve (x = raman, y = true)
% y = @(x) 10.35 .* x; 
% 
% TFV_RelConc = [TFVContr1' TFVContr2' TFVContr3']; 
% TFV_Conc = y(TFV_RelConc);
% minnum = min(TFV_Conc);
% minnum = minnum(1).*ones(size(TFV_Conc));
% TFV_Concshift = TFV_Conc-minnum;
% 
% %% Plot linear connected points to see trend
% CM = lines(5);
% figure(6); clf
% hold on
% ylabel('TFV Concentration (%w/w)')
% dex = [2 2 2 2 2];
% err = 0.0613.*ones(size(TFV_Conc));
% for i=1:5;
%     plot(xcrop1(:,[1:end-dex(i)])+shift(i), TFV_Conc([1:end-dex(i)],i), '-o', 'color', CM(i,:), 'LineWidth', 2.5, ...
%         'MarkerSize', 10)
%     %errorbar(xcrop1(:,[dex(i):end]), TFV_Conc([dex(i):end],i), err([dex(i):end],i)) 
% end
% hold on
% plot(xcrop1, ones(size(xcrop1)).*0.0613, 'k--')
% legend('Scan1 (12 min)', 'Scan2 (34 min)', 'Scan3 (57 min)', 'Scan4 (80 min)', 'Scan5(107 min)', 'LOD')%
% xlabel('x Depth (\mum)'); 
% title('TFV Concentrations in Porcine Rectal Tissue')
% axis([-100 800 0 3])
% 
% %% Plot points for polyfit
% clear P
% 
% CM = lines(5);
% xval = xcrop1(:,[1:end-3]);
% yval = TFV_Conc([1:end-3], :)';
% xfit = linspace(0, xval(:,1), 1000);
% 
% 
% for k=1:5;
% P(k,:) = polyfit(xval, yval(k,:), 5);
% yfit(k,:) = polyval(P(k,:),xfit);
% end
% 
% figure(7); clf
% hold on
% ylabel('TFV Conc (%w/w)')
% for i=1:5;
%     plot(xval+shift(i), yval(i,:), 'o', 'color', CM(i,:), 'LineWidth', 2.5)
%     hold on
% end
% for i=1:5;
%     plot(xfit+shift(i), yfit(i,:), '-', 'color', CM(i,:), 'LineWidth', 2.5);
% end
% legend('Scan1 (12 min)', 'Scan2 (34 min)', 'Scan3 (57 min)', 'Scan4 (80 min)', 'Scan5(107 min)')
% xlabel('x Depth (\mum)'); 
% title('TFV Concentrations in Porcine Rectal Tissue')
% axis([0 800 0 2])
% 



% 
% % Import Basis Spectra
% VB1 = load('TissueBasis.mat');
% %VB4 = load('VB4.mat');
% %VB5 = load('VB5.mat');
% %size(VB4)
% %size(VB5)
% %Ethanol = load('Ethanol.mat');
% Ethanol = load('EthanolPrePost.mat');
% systemBkgd = load('SystemBkgd.mat');
% %RB_L = load('/Users/Aubrey/Documents/Katz Lab/Summer/Rectal Basis/RectalBasisSpec.mat');
% 
% %spectraRB = [VB1.data.speccrop VB4.data.speccrop VB5.data.speccrop];
% spectraRB = VB1.data.speccrop;
% waveRB=VB1.data.wavenumcrop(:,:);
% spectraOH = Ethanol.data.speccrop(:,2);                   
% waveOH = Ethanol.data.wavenumcrop(:,:);
% spectraBkgd = systemBkgd.data.speccrop(:,:);
% waveBkgd = systemBkgd.data.wavenumcrop(:,:);
% 
% figure()
% plot(waveRB,spectraRB)
% 
% % Smooth and Normalize basis
% for k = 1:size(spectraRB,2)
%     %smoothRB(:, k) = smooth(spectraRB(:, k), sdex, 'sgolay', 2);
%     smoothRB(:, k) = spectraRB(:,k);
%     normRB(:, k) = smoothRB(:, k)/mean(smoothRB(:,k));
% end
% 
% for k = 1:size(spectraOH,2)
%     %smoothgel(:,k) = smooth(spectragel(:,k), sdex, 'sgolay', 2);
%     smoothOH(:,k) = spectraOH(:,k);
%     normOH(:,k) = smoothOH(:,k)/mean(smoothOH(:,k));
% end
% 
% for k = 1:size(spectraBkgd,2)
%     %smoothBkgd(:, k) = smooth(spectraBkgd(:, k), sdex, 'sgolay', 2);
%     smoothBkgd(:, k) = spectraBkgd(:,k);
%     normBkgd(:, k) = smoothBkgd(:, k)/mean(smoothBkgd(:,k));
% end
% %normBkgd = normBkgd(:,[1,2]);
% 
% %x1 = [3,5.5,7.5,10,12.5,14.5,17.5,20,22.5,25,27,29.5,32.5,34.5,36.5,38.5,40.5,42.5,44.5,51.5]; 
% 
% x1 = [1.5,4,6,8,10,14,16,18,22,24.5,27,29,31,33,35,36.5,38,40,42,43.5,46,50,53.5,55,57,59.5,62.5,64,66]; 
% %x1 = [2.5,6,8.5,10.5,13,15.5,17.5,19.5,21.5,23.5,25.5,27.5,29.5,31.5,33,35,37,39,41,43.5,45.5,47.5,49.5,51.5,53.5,55.5,57,59,61,63,67,69.5]; 
% 
% %size(x1)
% %size(normspecOT)
% 
% %%% Align all the RB basis and experiment spectra using the Phe peak
% List = {normRB};
% t_f = 0;
% for j = 1:size(List,2)
% dexU = find(wavenumOT>=1012, 1); 
% dexL = find(wavenumOT>=965, 1); 
% var = cell2mat(List(j));
% for i = 1:size(var,2)
%     Phepeakdex(i) = (dexL-1) + find(var([dexL:dexU],i)==max(var([dexL:dexU],i)));
%     diff(i) = wavenumOT(Phepeakdex(i)) - 1004;
%     tstart = tic;
%     while abs(diff(i))>=1 && t_f<=1
%         if diff(i)>0
%             var(:,i) = fshift(var(:,i), -0.1);
%         else diff(i)<0;
%             var(:,i) = fshift(var(:,i), 0.1);
%         end
%         t_f = toc(tstart);
% 
%         Phepeakdex(i) = (dexL-1) + find(var([dexL:dexU],i)==max(var([dexL:dexU],i)));
%         diff(i) = wavenumOT(Phepeakdex(i)) - 1004;
%     end
%     Listnew(j) = {var};
% end
% end
% 
% normRB = cell2mat(Listnew(1));
% 
% 
% shiftX = -4;
% normOH = fshift(normOH,shiftX);
% % shiftX = -12;
% % normRB = fshift(normRB,shiftX);
% 
% checkIndeces = [];
% checkIndecesplot = zeros(1,length(checkIndeces));
% 
% for i = 1:length(checkIndeces)
% checkIndecesplot(i) = find(x1==checkIndeces(i));
% end
% 
% normspecOT(:,find(x1==66)) = [];
% x1(:,find(x1==66)) = [];
% %
% % elim = [51.5];
% % normspecOT
% % x1(:,find(x1==35)) = [];
% % normspecOT(:,find(x1==54)) = [];
% % x1(:,find(x1==54)) = [];
% % normspecOT(:,find(x1==60)) = [];
% % x1(:,find(x1==60)) = [];
% % normspecOT(:,find(x1==55.5)) = [];
% % x1(:,find(x1==55.5)) = [];
% 
% % elim = 35;
% % elimInd = zeros(1,length(elim));
% % 
% % for i = 1:length(elim)
% % elimInd(i) = find(x1==elim(i));
% % end
% % %normspecOT(:,12) = [];
% % 
% % normspecOT(:,elimInd) = [];
% % x1(:,elimInd) = [];
% 
% % shiftX = 1;
% % normspecOT(:,31) = fshift(normspecOT(:,31),shiftX);
% % normspecOT(:,32) = fshift(normspecOT(:,32),shiftX);
% 
% %  shiftX = 1;
% %  normspecOT(:,31) = fshift(normspecOT(:,33),shiftX);
% %  normspecOT(:,32) = fshift(normspecOT(:,34),shiftX);
% % 
% % shiftX = 2;
% % normspecOT(:,25) = fshift(normspecOT(:,25),shiftX);
% % normspecOT(:,26) = fshift(normspecOT(:,26),shiftX);
% % normspecOT(:,27) = fshift(normspecOT(:,27),shiftX);
% % normspecOT(:,28) = fshift(normspecOT(:,28),shiftX);
% % normspecOT(:,29) = fshift(normspecOT(:,29),shiftX);
% 
% 
% 
% %shiftX=1;
% %normspecOT(:,31) = fshift(normspecOT(:,31),shiftX);
% 
% %shiftX=1;
% %normspecOT(:,32) = fshift(normspecOT(:,32),shiftX);
% %  
% %  shiftX=4.6;
% %  normspecOT(:,11) = fshift(normspecOT(:,11),shiftX);
% %  
% %  shiftX=4.6;
% %  normspecOT(:,12) = fshift(normspecOT(:,12),shiftX);
% %  
% %  shiftX=-20.7;
% %  normspecOT(:,13) = fshift(normspecOT(:,13),shiftX);
% % 
% %  shiftX=4;
% %  normspecOT(:,14) = fshift(normspecOT(:,14),shiftX);
% %  
% %  shiftX=-21.84;
% %  normspecOT(:,15) = fshift(normspecOT(:,15),shiftX);
% %  
% %  shiftX= -19.84;
% %  normspecOT(:,16) = fshift(normspecOT(:,16),shiftX);
% % 
% %  shiftX= -19.84;
% %  normspecOT(:,17) = fshift(normspecOT(:,17),shiftX);
% %  
% %  shiftX= -18.84;
% %  normspecOT(:,18) = fshift(normspecOT(:,18),shiftX);
% %  
% %  shiftX= 3;
% %  normspecOT(:,19) = fshift(normspecOT(:,19),shiftX);
% %  
% %  shiftX= 3;
% %  normspecOT(:,20) = fshift(normspecOT(:,20),shiftX);
% % 
% %  shiftX= 3;
% %  normspecOT(:,21) = fshift(normspecOT(:,21),shiftX);
%  
% %  shiftX=-21.84;
% %  normspecOT(:,16) = fshift(normspecOT(:,16),shiftX);
%  
% %tempvar = 5;
% 
% % one= normgel(:,2);
% % %two=normL(:,5);
% % three=normRB(:,1); %waveRB
% % four = normTFV(:,1);
% % five = normBkgd(:,1);
% % six = normG(:,1);
% 
% % figure(2); clf
% % plot(wavenuml2, Scan2(:,2))
% % hold on
% % plot(wavenuml2, five, 'LineWidth', 1)
% % xlabel('Raman Shift (cm^{-1})'); ylabel('Raman Intensity (a.u.)')
% % title('Raman Signature of Gel')
% % axis('tight')
% % legend('Scan', 'Basis')
% 
% %plot(wavenum, Scan3(:,1))
% 
% % Scan1 = Scan1(:, [1:3, 8:15, 18:end]);
% % Scan2 = Scan2(:, [1:3, 8:15, 18:end]);
% % Scan3 = Scan3(:, [1:3, 8:15, 18:end]);
% % Scan4 = Scan4(:, [1:3, 8:15, 18:end]);
% 
% %%
% 
% 
%   speccrop = [1:1174];
%   Scan1crop = normspecOT(speccrop,:);
% % % Scan2crop = Scan2(speccrop,:);
% % % Scan3crop = Scan3(speccrop,:);
%  wavenumcrop = wavenumOT(speccrop,:);
% %  xcrop1 = x1(:,:);%(:, [1:3, 8:15, 18:end]);
% % % xcrop2 = x23(:,:);%(:, [1:3, 8:15, 18:end]);
% % % xcrop3 = x23(:,:);%(:, [1:3, 8:15, 18:end]);
% %% seuqenctial polot
% % figure(1); clf
% % for i=1:size(Scan2,2);
% %     plot(wavenum, Scan2(:,i))
% %     hold on
% %     title(num2str(i))
% %     pause
% % end
% % 
% %% establish a stacked plot
% 
% figure(); clf
% hold on
% CM = colormap(winter(28));
% CM = flipud(CM);
% k=0;
% [rr,cc] = size(normspecOT);
% indices = [1,3,4,8,15,18];
% for i = 1:length(indices)
%     k=k+1;
%     plot(wavenumOT, normspecOT(:,indices(i))-12.*k, 'color', CM(k*3,:), 'LineWidth', 2)
% end
% xlabel('Wavenumber (cm^{-1})','FontSize',18,'FontWeight','Bold');
% ylabel('\leftarrow Increasing Time','FontSize',18,'FontWeight','Bold');  
% title('Sample spectra on stromal surface over time','FontSize',18,'FontWeight','Bold')
% set(gca,'YTickLabel',[]);
% XTickL = get(gca,'XTick');
% set(gca,'FontSize',14);
% %set(gca,'XTickLabel','FontSize',16);
% ylim([-76,-5])
% xlim([600,1800])
% labels = {'0 \mum','25 \mum','50 \mum','75 \mum' '100 \mum','125 \mum','150 \mum', '200 \mum', '250 \mum'};
% labels=flipud(labels);
% %lcolorbar(labels,'TitleString','Depth');
% %axis('tight')
% 