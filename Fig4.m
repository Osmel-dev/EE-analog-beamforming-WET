% This script can be used to reproduce Fig. 4 in [R1]
clear; clc;

%% Simulation settings
Eth = [.004 .004 .004]';                    % energy requirements (J)
K = numel(Eth);                             % number of users
Nt = [16 64];                               % number of PB's antennas
tauMaxTimeDiv = linspace(2.1,15,2);         % maximum charging time for time division (s)
tauMaxAllAtOnce = linspace(.7,6,2);         % maximum charging time for all-at-once (s)
MCruns = 2;                                 % Monte Carlo runs
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

% devices' positions matrix, e.g., s-th row corresponds to (x_s, y_s, z_s)
% Cartesian coordinates (in meters)
sensorPos = [0 -5 0; 5 0 -5; 10 10 5];

%% Monte Carlo simulation

% Memory allocation
EOptOrderAware= zeros(numel(tauMaxTimeDiv),numel(Nt));
EOptOrderAgnostic = zeros(numel(tauMaxTimeDiv),numel(Nt));
EOptAllAtOnce = zeros(numel(tauMaxAllAtOnce),numel(Nt));

% Monte Carlo loop
for nt = 1:numel(Nt)
    for tt = 1:numel(tauMaxTimeDiv)
        % Memory allocation
        EOptOrderAwareTemp = zeros(MCruns,1);
        EOptOrderAgnosticTemp = zeros(MCruns,1);
        EOptAllAtOnceTemp = zeros(MCruns,1);
    
        for mc = 1:MCruns
            disp(mc)
            % Channel realization URA
            h = channelModelURA(K,Nt(nt),mc,kappa,sensorPos);
            
            % Optimal solution per channel realization
            % Tip! To reduce simulation time, notice that the function
            % orderAware(-) departs from solving the optimization problem
            % P2 which is implemented by orderAgnostic(-)
            [EOptOrderAwareTemp_,~,~] = orderAware(h,Nt(nt),K,Eth,tauMaxTimeDiv(tt),a,w0,wsat,g,ginv);
            [EOptOrderAgnosticTemp_,~,~] = orderAgnostic(h,Nt(nt),K,Eth,tauMaxTimeDiv(tt),a,w0,wsat,ginv);
            [EOptAllAtOnceTemp_,~,~] = SDMAHeuristic(h,Nt(nt),K,Eth,tauMaxAllAtOnce(tt),ginv);
        
            EOptOrderAwareTemp(mc) = EOptOrderAwareTemp_;
            EOptOrderAwareTemp(mc) = EOptOrderAgnosticTemp_;
            EOptAllAtOnceTemp(mc) = EOptAllAtOnceTemp_;    
        end
    
        % Optimal solution (average)
        EOptOrderAware(tt,nt) = mean(nonzeros(EOptOrderAwareTemp),'omitnan');
        EOptOrderAgnostic(tt,nt) = mean(nonzeros(EOptOrderAgnosticTemp),'omitnan');
        EOptAllAtOnce(tt,nt) = mean(nonzeros(EOptAllAtOnceTemp),'omitnan');
    end
end

%% Plot results
% IMPORTANT! Reproducing Fig. 4 as in [R1] requires more Monte Carlo
% iterations. This is just an example to illustrate the structure of the
% code

fig = figure;
plotSettings(fig)

% Nt = 16
subplot(2,1,1);
semilogy(tauMaxTimeDiv, EOptOrderAgnostic(:,1),'-ok','LineWidth',1.5); hold on
semilogy(tauMaxTimeDiv, EOptOrderAware(:,1),'-s','LineWidth',1.5,'Color','#77AC30');
semilogy(tauMaxAllAtOnce, EOptAllAtOnce(:,1),'--+','LineWidth',1.5,'Color','#A2142F'); 

% annotations
annotation('line',[0.220055710306407 0.221448467966574],...
    [0.920529801324504 0.596026490066225],'LineWidth',1,'LineStyle','-.'); 

annotation('line',[0.448467966573816 0.449860724233983],...
    [0.916114790286976 0.598233995584988],'LineWidth',1,'LineStyle','-.');

annotation('arrow',[0.52924791086351 0.427576601671309],...
    [0.666666666666667 0.737306843267108],'LineWidth',1);

title('N=16','fontsize',12,'fontweight','normal','Interpreter','tex')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold off

grid on
xlim([0 6])
ylim([26 83])
xlabel('$\tau_\mathrm{max}$ [s]', 'Interpreter','latex', 'FontSize',15)
ylabel('E [J]', 'Interpreter','latex', 'FontSize',15)
legend('order-agnostic', 'order-aware', 'charging all-at-once','fontsize',15,...
    'Position',[0.585808064395945 0.683878934535297 0.310198508442837 0.115894036577237],'interpreter','latex')

ax1 = axes('Parent',gcf,'Position',[0.26183844011142 0.664459161147903 0.165738161559889 0.154525386313465]);
semilogy(tauMaxTimeDiv, EOptOrderAgnostic(:,1),'-ok','LineWidth',1.5); hold on
semilogy(tauMaxTimeDiv, EOptOrderAware(:,1),'-s','LineWidth',1.5,'Color','#77AC30'); hold off
grid on
xlim([2.5 6])

% Nt = 64
subplot(2,1,2)
semilogy(tauMaxTimeDiv, EOptOrderAgnostic(:,2),'-ok','LineWidth',1.5); hold on
semilogy(tauMaxTimeDiv, EOptOrderAware(:,2),'-s','LineWidth',1.5,'Color','#77AC30');
semilogy(tauMaxAllAtOnce, EOptAllAtOnce(:,2),'--+','LineWidth',1.5,'Color','#A2142F'); 

annotation('line',[0.220055710306407 0.221448467966574],...
    [0.448123620309051 0.123620309050772],'LineWidth',1,'LineStyle','-.');

annotation('line',[0.444289693593315 0.445682451253482],...
    [0.448123620309051 0.123620309050773],'LineWidth',1,'LineStyle','-.');

annotation('arrow',[0.525069637883008 0.426183844011142],...
    [0.17439293598234 0.278145695364238],'LineWidth',1);

title('N=64','fontsize',12,'fontweight','normal','Interpreter','tex')
ax = gca;
ax.TitleHorizontalAlignment = 'left';

hold off
grid on
xlim([0 6])
ylim([7 23])
xlabel('$\tau_\mathrm{max}$ [s]', 'Interpreter','latex', 'FontSize',15)
ylabel('E [J]', 'Interpreter','latex', 'FontSize',15)
legend('order-agnostic', 'order-aware', 'charging all-at-once','fontsize',15,...
    'Position',[0.581641316033817 0.206902433596809 0.310198508442837 0.115894036577237],'interpreter','latex')

ax2 = axes('Parent',gcf,'Position',[0.261838440111421 0.187637969094923 0.164345403899721 0.161147902869757]);
semilogy(tauMaxTimeDiv, EOptOrderAgnostic(:,2),'-ko','LineWidth',1.5); hold on
    semilogy(tauMaxTimeDiv, EOptOrderAware(:,2),'-s','LineWidth',1.5,'Color','#77AC30'); hold off
grid on
xlim([2.5 6])

%% Refs.
% [R1] O. M. Rosabal, et al., "Energy-Efficient Analog Beamforming for ...
% RF-WET with Charging Time Constraint", IEEE Trans. Veh. Technol. (submitted), 2023.