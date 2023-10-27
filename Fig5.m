% This script can be used to reproduce Fig. 5 in [R1]
clear; clc;

%% Simulation settings
Eth = [.004 .004 .004]';                    % energy requirements (J)
K = numel(Eth);                             % number of users
Nt = [16 64];                               % number of PB's antennas
tauMax = 8;                                 % maximum charging time (s)
MCruns = 2;                                 % Monte Carlo runs (recommended 1e4)
kappa = 10;                                 % Rician LoS factor (linear)

% Powercast P2110B RF-EH circuit data
P2110B = readmatrix('data/P2110B_overallEff@868.csv'); % data file
a = [-0.001952 0.663 -0.01453];                        % coefficients of the quadratic model
wsat = 10.^(P2110B(end,1)/10);                         % saturation level [mW]

[~, idx] = max(P2110B(:,2));                           
w0 = 10.^(P2110B(idx,1)/10);                           % input RF power max. conversion efficiency [mW]

% linear interpolant model (benchmark)
load g.mat                               % RF-EH transfer function
load ginv.mat                            % inverse of RF-EH transfer function

% Sensor deployment
alpha = linspace(0,1,3);                            % step size

s1 = (1-alpha).*[0 -5 0]' + alpha.*[5 2.5 2.5]';    % 1st device Cartesian Coordinates [m]
s2 = (1-alpha).*[5 0 -5]' + alpha.*[5 2.5 2.5]';    % 2nd device Cartesian Coordinates [m]
s3 = (1-alpha).*[10 10 5]' + alpha.*[5 2.5 2.5]';   % 3rd device Cartesian Coordinates [m]

%% Monte Carlo simulation

% Memory allocation
EOptOrderAware= zeros(numel(alpha),numel(Nt));
EOptOrderAgnostic = zeros(numel(alpha),numel(Nt));
EOptAllAtOnce = zeros(numel(alpha),numel(Nt));

% Monte Carlo loop
for nt = 1:numel(Nt)
    for ss = 1:numel(alpha)
        % Memory allocation
        EOptOrderAwareTemp = zeros(MCruns,1);
        EOptOrderAgnosticTemp = zeros(MCruns,1);
        EOptAllAtOnceTemp = zeros(MCruns,1);
    
        sensorPos = [s1(:,ss)'; s2(:,ss)'; s3(:,ss)'];
        for mc = 1:MCruns
            disp(mc)
            % Channel reaization URA
            h = channelModelURA(K,Nt(nt),mc,kappa,sensorPos);
    
            % Optimal solution per channel realization
            % Tip! To reduce simulation time, notice that the function
            % orderAware(-) departs from solving the optimization problem
            % P2 which is implemented by orderAgnostic(-)
            [EOptOrderAwareTemp_,~,~] = orderAware(h,Nt(nt),K,Eth,tauMax,a,w0,wsat,g,ginv);
            [EOptOrderAgnosticTemp_,~,~] = orderAgnostic(h,Nt(nt),K,Eth,tauMax,a,w0,wsat,ginv);
            [EOptAllAtOnceTemp_,~,~] = SDMAHeuristic(h,Nt(nt),K,Eth,tauMax,ginv);
        
            EOptOrderAwareTemp(mc) = EOptOrderAwareTemp_;
            EOptOrderAwareTemp(mc) = EOptOrderAgnosticTemp_;
            EOptAllAtOnceTemp(mc) = EOptAllAtOnceTemp_;    
        end
    
        % Optimal solution (average)
        EOptOrderAware(ss,nt) = mean(nonzeros(EOptOrderAwareTemp),'omitnan');
        EOptOrderAgnostic(ss,nt) = mean(nonzeros(EOptOrderAgnosticTemp),'omitnan');
        EOptAllAtOnce(ss,nt) = mean(nonzeros(EOptAllAtOnceTemp),'omitnan');
    end
end


%% Plot results
% IMPORTANT! Reproducing Fig. 4 as in [R1] requires more Monte Carlo
% iterations. This is just an example to illustrate the structure of the
% code

fig = figure;
plotSettings(fig)
fontSize = 15;

semilogy(alpha,EOptOrderAgnostic,'-ok','LineWidth',1.5); hold on
semilogy(alpha,EOptOrderAware, '-s','LineWidth',1.5,'Color','#77AC30')
semilogy(alpha,EOptAllAtOnce,'--+','LineWidth',1.5,'Color','#A2142F'); 

semilogy(alpha,EOptOrderAgnostic,'-ok','LineWidth',1.5);
semilogy(alpha,EOptOrderAware,'-s','LineWidth',1.5,'Color','#77AC30')
semilogy(alpha,EOptAllAtOnce,'--+','LineWidth',1.5,'Color','#A2142F'); 

% Fig. annotations
annotation('ellipse',...
    [0.755874651810585 0.196467991169978 0.0171058495821729 0.101545253863135],...
    'LineWidth',1);

annotation('ellipse',...
    [0.743339832869081 0.355408388520971 0.0171058495821729 0.101545253863135],...
    'LineWidth',1);

annotation('textbox',...
    [0.662952646239554 0.143487858719647 0.129526462395544 0.0519801324503309],...
    'String','N = 64',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off');

annotation('textbox',...
    [0.747910863509749 0.45916114790287 0.129526462395544 0.051980132450331],...
    'String',{'N = 16'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off');

hold off
grid on
xlabel('$\alpha$','FontSize',fontSize,'Interpreter','latex')
ylabel('E [J]','FontSize',fontSize,'Interpreter','latex')
legend('order-agnostic','order-aware','all-at-once strategy','FontSize',fontSize,...
    'interpreter','latex')

%% Refs.
% [R1] O. M. Rosabal, et al., "Energy-Efficient Analog Beamforming for ...
% RF-WET with Charging Time Constraint", IEEE Trans. Veh. Technol. (submitted), 2023.