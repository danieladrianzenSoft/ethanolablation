
clear;
clc;
% Load Depth Scan spectra 

EtOh_100 = load('100_EtOH.mat');
EtOh_75 = load('75_EtOH.mat');
EtOh_50 = load('50_EtOH.mat');
EtOh_25 = load('25_EtOH.mat');
EtOh_0 = load('0_EtOH.mat');
Cyclo = load('cyclo.mat');

wavenumEtOh_100 = EtOh_100.data.wavenumcrop;
wavenumEtOh_75 = EtOh_75.data.wavenumcrop;
wavenumEtOh_50 = EtOh_50.data.wavenumcrop;
wavenumEtOh_25 = EtOh_25.data.wavenumcrop;
wavenumEtOh_0 = EtOh_0.data.wavenumcrop;
wavenumCyclo = EtOh_0.data.wavenumcrop;
spectraEtOh_100 = EtOh_100.data.speccrop;
spectraEtOh_75 = EtOh_75.data.speccrop;
spectraEtOh_50 = EtOh_50.data.speccrop;
spectraEtOh_25 = EtOh_25.data.speccrop;
spectraEtOh_0 = EtOh_0.data.speccrop;
spectraCyclo = EtOh_0.data.speccrop;

for k = 1:size(EtOh_100,2)
    %smoothspecPDMS(:, k) = smooth(spectraPDMS(:, k), sdex, 'sgolay', 2);
    spectraEtOh_100(:, k) = spectraEtOh_100(:,k);
    normSpecEtOh_100(:, k) = spectraEtOh_100(:, k)/mean(spectraEtOh_100(:,k));
%     normSpecEtOh_100(:, k) = spectraEtOh_100(:, k);

end
for k = 1:size(EtOh_75,2)
    %smoothspecPDMS(:, k) = smooth(spectraPDMS(:, k), sdex, 'sgolay', 2);
    spectraEtOh_75(:, k) = spectraEtOh_75(:,k);
    normSpecEtOh_75(:, k) = spectraEtOh_75(:, k)/mean(spectraEtOh_75(:,k));
%     normSpecEtOh_75(:, k) = spectraEtOh_75(:, k)/mean(spectraEtOh_75(:,k));
end
for k = 1:size(EtOh_50,2)
    %smoothspecPDMS(:, k) = smooth(spectraPDMS(:, k), sdex, 'sgolay', 2);
    spectraEtOh_50(:, k) = spectraEtOh_50(:,k);
    normSpecEtOh_50(:, k) = spectraEtOh_50(:, k)/mean(spectraEtOh_50(:,k));
end
for k = 1:size(EtOh_25,2)
    %smoothspecPDMS(:, k) = smooth(spectraPDMS(:, k), sdex, 'sgolay', 2);
    spectraEtOh_25(:, k) = spectraEtOh_25(:,k);
    normSpecEtOh_25(:, k) = spectraEtOh_25(:, k)/mean(spectraEtOh_25(:,k));
end
for k = 1:size(EtOh_0,2)
    %smoothspecPDMS(:, k) = smooth(spectraPDMS(:, k), sdex, 'sgolay', 2);
    spectraEtOh_0(:, k) = spectraEtOh_0(:,k);
    normSpecEtOh_0(:, k) = spectraEtOh_0(:, k)/mean(spectraEtOh_0(:,k));
end
for k = 1:size(Cyclo,2)
    %smoothspecPDMS(:, k) = smooth(spectraPDMS(:, k), sdex, 'sgolay', 2);
    spectraCyclo(:, k) = spectraCyclo(:,k);
    normSpecCyclo(:, k) = spectraCyclo(:, k)/mean(spectraCyclo(:,k));
end


figure()
plot(wavenumEtOh_100, normSpecEtOh_100)
ylabel('Normalized Spectrum','FontSize',16,'FontWeight','Bold')
xlabel('Wavenumber','FontSize',16,'FontWeight','Bold')
hold on
plot(wavenumEtOh_75, normSpecEtOh_75)
hold on
plot(wavenumEtOh_50, normSpecEtOh_50)
hold on
plot(wavenumEtOh_25, normSpecEtOh_25)
hold on
plot(wavenumEtOh_0, normSpecEtOh_0)
legend('100','75','50','25','0')

figure()
plot(wavenumCyclo, normSpecCyclo)
ylabel('Cyclo','FontSize',16,'FontWeight','Bold')
xlabel('Wavenumber','FontSize',16,'FontWeight','Bold')

%% Establish Basis array
% Basis Array for positive fitting 
% C(:, 1:272) = normRB(:,:);
% C(:, 273:282) = normgel(:,:);
% C(:, 283:329) = normL(:,:);
% C(:, 330) = normTFV(:,:);
% 
% [RR,CC] = size(normspecLiverBlob);
% 
% Livervar = [1:CC];
% % Bkgdvar = [CC+1];
% %Lvar = [188:250];
% OHvar = [CC+2];

%size(RBvar)
%size(normRB)

%[tempsizeR,tempsizeC] = size(normRB);

% C = zeros(length(normspecLiverBlob),CC+2);
% %size(normRB(speccrop,:))
% %normOH
% C(:, Livervar) = normspecLiverBlob(:,:);
% % C(:, OHvar) = normBkgd(speccrop,:);
% % C(:, Lvar) = normL(speccrop,:);  
% C(:, OHvar) = normspecEthanol(:,:);

% figure()
% plot(speccrop,C(:, OHvar))

%% Chemometric Fitting
% Find the coeeficients of the matrix computation [Fit] = [Basis][Coeff]

% Scan 1
% for k = 1:size(normspecLiverBlobInjectedEthanol,2);
%     coeff1(1:size(C,2), k) = lsqnonneg(C,normspecLiverBlobInjectedEthanol(:,k));
% end
% 
% %% Fitting and Determining Relative Contributions - Positive Fit
% % Find the fit spectra and smooth
% % Find the tissue Contribution and the TFV contribution
% 
% % Scan 1 - Tissue 
% for k = 1:size(normspecLiverBlobInjectedEthanol,2);
% %for k = 1:5;
%     fit1(:,k) = C(:,:)*coeff1(:,k);
%     %fit1(:,k) = smooth(fit1(:,k), 20, 'sgolay',2);
%     tissueContr1(:,k) = C(:, Livervar)*coeff1(Livervar,k);
%     OHbasisContr1(:,k) = C(:, OHvar)*coeff1(OHvar, k);
%     OHContr1(:,k) = trapz(OHbasisContr1(:,k))./(trapz(tissueContr1(:,k)));
%     Tissue1(:,k) = trapz(tissueContr1(:,k))./trapz(fit1(:,k));
% 
% end
% 
% %%
