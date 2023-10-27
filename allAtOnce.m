function [EOpt,wOpt,tOpt] = allAtOnce(chn,Nt,K,Eth,tauMax,ginv)
    % This funtion implements the solution of P3 in Section III.C [R1].
    % In the following, we define a list of arguments and outputs of the
    % function. In square brackets, numerical domain, dimensions, and/or
    % unit of the corresponding variable depending on the case. 
    % =====================================================================

    % channel matrix SDP formulation
    H = zeros(Nt,Nt,K);
    for kk = 1:K
        H(:,:,kk) = chn(:,kk)*chn(:,kk)';  
    end
    
    %% Optimization problem P3 (SDP formulation)
    cvx_begin quiet
        variable W(Nt,Nt) complex semidefinite
        variable inRF
        maximize ( inRF )
        subject to
            for kk = 1:K
                real(trace(W*H(:,:,kk))) >= inRF;
            end
    cvx_end

    % Recoverig the SDP's solution via dominant eigenvector
    [V,D] = eig(W);
    [~, maxEigIdx] = max(diag(D));
    wTemp = exp(1i*angle(V(:,maxEigIdx))); % optimal phase shift

    % Solve the power allocation for a fixed time (tau)
    tauMin = max(Eth)/5.911*1e3;
    tauVector = linspace(tauMin,tauMax,1e3)';
    p = zeros(numel(tauVector),1);

    % Re-adjust tx power using linear interpolant model for RF-EH
    % (benchmarking purposes) and eq. (13) in [R1]
    for tt = 1:numel(tauVector)
        % Memory allocation
        p_ = zeros(K,1);
        for kk = 1:K
            p_(kk) = ginv(min(Eth(kk)/tauVector(tt)*1e3,5.911))*1e-3*Nt/abs(chn(:,kk)'*wTemp)^2;
        end
        p(tt) = max(p_);
    end

    %% Retrieve final solution
    [EOpt, idx]= min(p.*tauVector);    
    tOpt = tauVector(idx);
    wOpt = sqrt(p(idx)/Nt)*wTemp;
end

%% Refs.
% [R1] O. M. Rosabal, et al., "Energy-Efficient Analog Beamforming for ...
% RF-WET with Charging Time Constraint", IEEE Trans. Veh. Technol. (submitted), 2023.