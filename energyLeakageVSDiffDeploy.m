%% TDMA vs SDMA in terms of devices' deployment 
clear; clc;

%% Simulation settings
E = [.004 .004 .004]';                    % energy requirements (J)
K = numel(E);                                  % number of users
Nt = [8 8];                                     % number of antennas
tauMax = 8;                               % maximum charging time
MCruns = 1e4;                                     % Monte Carlo runs
kappa = 10;                                     % Rician factor

%% EH transfer function (Powercast P2110B)
load ginv.mat
load g.mat

% TDMA model (quadratic model)
P2110B = readmatrix('data/P2110B_overallEff@868.csv');
a = [-0.001952 0.663 -0.01453];
wsat = 10.^(P2110B(end,1)/10);
[~, idx] = max(P2110B(:,2));
w0 = 10.^(P2110B(idx,1)/10);

%% Sensor deployments
alpha = [0 .5 .9];

s1 = (1 - alpha).*[0 -5 0]' + alpha.*[5 2.5 2.5]';
s2 = (1 - alpha).*[5 0 -5]' + alpha.*[5 2.5 2.5]';
s3 = (1 - alpha).*[10 10 5]' + alpha.*[5 2.5 2.5]';

Er = zeros(K,K,numel(alpha));

parfor ss = 1:numel(alpha)
    disp(num2str(ss))
    Er_ = zeros(K,K,MCruns);

    sensorPos = [s1(:,ss)'; s2(:,ss)'; s3(:,ss)'];
    for mc = 1:MCruns
        %% Channel model
        [h,~,~] = channel_model_URA(K,Nt,mc,kappa,sensorPos);

        %% Solution
        [~,wOptVecTDMA,tOptVecTDMA] = gpTDMA(h,prod(Nt),K,E,tauMax,a,w0,wsat,ginv);

        for kk = 1:K
            for jj = 1:K
                Er_(kk,jj,mc) = g(min(abs(wOptVecTDMA(:,jj)'*h(:,kk))^2*1e3,12.1269))*1e-3*tOptVecTDMA(jj);
            end
        end
    end

    Er(:,:,ss) = mean(Er_,3,'omitnan');
end

disp('Finished')
save './data/TDMAXtraEnergy.mat'
exit