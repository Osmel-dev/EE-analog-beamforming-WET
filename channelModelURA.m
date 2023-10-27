function chn = channelModelURA(K,Nt,rndseed,kappa,sensorPos)
    % This function generates realizations of the wireless channel
    % according to the model in [R1,R2] for a uniform rectangular array
    % (URA). In the following, we define a list of arguments and outputs of the
    % function. In square brackets, numerical domain, dimensions, and/or
    % unit of the corresponding variable depending on the case. 
    % =====================================================================
    % List of input arguments:
    % K         : number of IoT devices [integer]
    % Nt        : number of PB's antennas [integer]
    % rndseed   : seed of the random number generator [integer]
    % kappa     : Rician LoS factor [real scalar, linear]
    % sensorPos : Sensor deployment [real matrix, K x 3, m]
    % =====================================================================
    % List of outputs:
    % chn       : channel realization [complex matrix, Nt x K]
    % =====================================================================

    % The URA is placed in the y-z plane: coordinates of the origin (0,0,5)
    Ntxy = sqrt(Nt);
    rng(rndseed)        
    
    % sensor positions in Polar coordinates
    [azimuth,elevation,r] = cart2sph(sensorPos(:,1), sensorPos(:,2), sensorPos(:,3));
    phi_k = azimuth + 2*pi;
    theta_k = pi/2 - elevation;

    % channel model according to [R1,R2]
    dx_k = [5 5 sqrt(5^2 + 10^2) sqrt(5^2 + 10^2)];      % distance between the s-th device and the PB    
    Psi_k = pi*sin(theta_k).*sin(phi_k);                 % y^ component
    Omega_k = pi*cos(theta_k);                           % z^ component    
    beta = 10^-1.6*max(r,1).^(-2.7);                     % log-distance path loss model

    % Rician channel samples
    chn = zeros(Nt,K);
    for kk = 1:K
        % horizontal & vertical channel components
        a_h = exp(2*1i*(0:(Ntxy-1))'*Psi_k(kk));
        a_v = exp(2*1i*(0:(Ntxy-1))'*Omega_k(kk));

        % LoS & nLoS channel components
        h_los = sqrt(kappa/(1+kappa))*kron(a_h,a_v);
        h_nlos = sqrt(1/(2*(1+kappa)))*(randn(Nt,1)+1i*randn(Nt,1));
        chn(:,kk) = sqrt(beta(kk))*(h_los + h_nlos);
    end 
end

%% Refs.
% [R1] Achievable Sum-Rate Analysis for Massive MIMO Systems with Different Array Configurations
% [R2] Balanis, C. A. (2015). Antenna theory: analysis and design. John wiley & sons.