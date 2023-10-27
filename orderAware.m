function [EOpt,wOpt,tOpt] = orderAware(h,Nt,K,Eth,tauMax,a,w0,wsat,g,ginv)
    % This function implements the Algorithm 1 "Order-aware time division".
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
    % g         : RF-EH transfer function (linear interpolant model) 
    %             [funcion handle]
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
                Eth(kk)/a(2)*Nt*t(kk)^(-1)*p(kk)^(-1)*norm(h(:,kk),1)^(-2) - a(3)*1e3/(a(2)*Nt)*p(kk)*norm(h(:,kk),1)^2 - ...
                    a(1)/a(2)*Nt*p(kk)^(-1)*norm(h(:,kk),1)^(-2)*1e-3 <= 1
               
                Nt*norm(h(:,kk),1)^(-2)*1e-3*w0*p(kk)^(-1) <= 1
                1/Nt*norm(h(:,kk),1)^2*1e3*wsat^(-1)*p(kk) <= 1
            end

            sum(t) <= tauMax;
    cvx_end

    % Re-adjust tx power using linear interpolant model for RF-EH (benchmarking purposes)
    p = zeros(K,1);
    for kk = 1:K
        p(kk) = ginv(min(Eth(kk)/t(kk)*1e3,5.911))*1e-3*Nt/abs(h(:,kk)'*exp(1i*angle(h(:,kk))))^2;
    end

    % Temporal beamforming vectors
    wTemp = sqrt(p/Nt)'.*exp(1i*angle(h));

    %% Order-aware apporach
    K_ = 1:K;
    for kk = 1:K-1
        % 1. Compute the interference in non-intended devices
        interf = zeros(numel(K_),1);
        for jj = 1:numel(K_)
            for ii = 1:numel(K_)
                if K_(ii) ~= K_(jj)
                    interf(jj) = interf(jj) + g(min(abs(wTemp(:,K_(jj))'*h(:,K_(ii)))^2*1e3,12.1269))*1e-3*t(K_(jj));
                end
            end
        end

        % 2. Select the beamforming that yields the highest interference (sum)
        [~,idx] = max(interf);
        
        % Remove the corresponding device whose beamforming yields the highest interference 
        usr = K_(idx);
        K_(idx) = [];

        % 3. Harvested energy at non-intended devices from interference
        for jj = K_
            % Remaining energy demand
            Eth(jj) = Eth(jj) - g(min(abs(wTemp(:,usr)'*h(:,jj))^2*1e3,12.1269))*1e-3*t(usr);
            
            % Adjust tx power and/or charging time according to the results Section III
            if t(jj) >= Eth(jj)/(g(w0)*1e-3)
                p(jj) = Nt*norm(h(:,jj),1)^(-2)*w0*1e-3;
                t(jj) = Eth(jj)/(g(w0)*1e-3);
            else
                p(jj) = ginv(min(Eth(jj)/t(jj)*1e3,5.911))*1e-3*Nt/abs(h(:,jj)'*exp(1i*angle(h(:,jj))))^2;
            end

            % Recompute the temporal beamforming for the remaining users
            wTemp(:,jj) = sqrt(p(jj)/Nt)'.*exp(1i*angle(h(:,jj)));
        end
    end

    %% Retrieve final solution
    EOpt = p'*t;                                
    wOpt = sqrt(p/Nt)'.*exp(1i*angle(h));       
    tOpt = t;                                   
end

%% Refs.
% [R1] O. M. Rosabal, et al., "Energy-Efficient Analog Beamforming for ...
% RF-WET with Charging Time Constraint", IEEE Trans. Veh. Technol. (submitted), 2023.