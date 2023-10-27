function [EOpt,wOpt,tOpt] = orderAgnostic(chn,Nt,K,Eth,tauMax,a,w0,wsat,ginv)
    % This function implements the solution of the optimization problem P2.
    % (see Section II.A "Order-agnostic time division" [R1])
    % In the following, we define a list of arguments and outputs of the
    % function. In square brackets, numerical domain, dimensions, and/or
    % unit of the corresponding variable depending on the case. 
    % =====================================================================
    % List of input arguments:
    % chn       : channel realization [complex matrix, Nt x K]
    % Nt        : number of PB's antennas [integer]
    % K         : number of IoT devices [integer]
    % Eth       : devices' energy requirements [real vector, K x 1, J]
    % tauMax    : maximum charging time [real, s]
    % a         : coefficients of RF-EH circuit model (quadratic model)
    % [real vector, 1 x 3]
    % w0        : sensitivity RF-EH circuit [real, mW]
    % wsat      : saturation RF-EH circuit [real, mW]
    % ginv      : inverse of the RF-EH transfer function (linear 
    %             interpolant model) [funcion handle]
    % =====================================================================
    % List of ouputs:
    % EOpt      : PB's energy consumption [real, J]
    % wOpt      : PB's beamformers [complex matrix, Nt x K]
    % tOpt      : optimum time allocation [real vector, s]
    % =====================================================================
    
    %% Optimization problem P2 (GP formulation)
    cvx_begin gp quiet
        variable p(K) nonnegative
        variable t(K) nonnegative
        minimize ( p'*t )
        subject to
            for kk = 1:K
                Eth(kk)/a(2)*Nt*t(kk)^(-1)*p(kk)^(-1)*norm(chn(:,kk),1)^(-2) - a(3)*1e3/(a(2)*Nt)*p(kk)*norm(chn(:,kk),1)^2 - ...
                    a(1)/a(2)*Nt*p(kk)^(-1)*norm(chn(:,kk),1)^(-2)*1e-3 <= 1
               
                Nt*norm(chn(:,kk),1)^(-2)*1e-3*w0*p(kk)^(-1) <= 1
                1/Nt*norm(chn(:,kk),1)^2*1e3*wsat^(-1)*p(kk) <= 1
            end

            sum(t) <= tauMax;
    cvx_end

    % Re-adjust tx power using linear interpolant model for RF-EH
    % (benchmarking purposes)
    for kk = 1:K
        p(kk) = ginv(min(Eth(kk)/t(kk)*1e3,5.911))*1e-3*Nt/abs(chn(:,kk)'*exp(1i*angle(chn(:,kk))))^2;
    end

    %% Retrieve final solution
    EOpt = p'*t;
    wOpt = sqrt(p/Nt)'.*exp(1i*angle(chn));
    tOpt = t;
end

%% Refs.
% [R1] O. M. Rosabal, et al., "Energy-Efficient Analog Beamforming for ...
% RF-WET with Charging Time Constraint", IEEE Trans. Veh. Technol. (submitted), 2023.