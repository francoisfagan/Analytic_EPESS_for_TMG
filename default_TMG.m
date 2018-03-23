function [ mu, Sigma, chol_Sigma, C, lB, uB ] = default_TMG(dimension)
% Define default settings for the mean, covariances and constraints 
% for the Truncated Multivariate Gaussian

% For details on how to define the truncated region, 
% refer to the code in the epmgp directory.
    
    % Define mean
    mu = zeros(dimension,1);
    
    % Shift the left boundary of truncate region to make sampling harder
    left_boundary = 0;

    % Define covariance, lower and upper bounds
    if dimension == 2
        % Covariance
        Sigma = [1 0; 0 1];
        
        % Bounds
        lB = [left_boundary; -1];           % Lower bound
        uB = [left_boundary + 1; 1];        % Upper bound
        
    else 
        % Covariance
        Sigma = eye(dimension);
        
        % Trunctation parameters
        box_length = 1;                     % Length of the box interval along each dimension
    
        % Bounds
        lB = zeros(dimension, 1);           % Lower bound
        lB(1) = left_boundary;              % Lower bound
        uB = lB + box_length.*ones(dimension, 1);   % Upper bound     
    end

    % Cholesky decomposition of covariance matrix
    chol_Sigma = chol(Sigma);

    % Define box constarints
    C = eye(dimension);

end