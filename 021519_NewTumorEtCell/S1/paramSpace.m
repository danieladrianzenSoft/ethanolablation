D_S = 4*10^(-7):1*10^(-8):1*10^(-6);
k_E = (0:0.25:2)/3600;
k_D = (0:0.25:2)/3600;

%D_S = 1:1:3;
%k_E = 4:1:7;
%k_D = 8:1:10;

mat = zeros(length(D_S)*length(k_E)*length(k_D),3);

ind = 1;

for i = 1:length(D_S)
    for j = 1:length(k_E)
        for k = 1:length(k_D)
            mat(ind,:) = [D_S(i),k_E(j),k_D(k)];
            ind = ind + 1;
        end
    end
end
 
csvwrite('paramSpace3.csv',mat)