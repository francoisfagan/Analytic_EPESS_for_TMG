% Code for sampling a Truncated Multivariate Gaussian using Analytic EPESS
% From the paper
% "Elliptical Slice Sampling with Expectation Propagation." 
% Fagan, Francois, Jalaj Bhandari, and John Cunningham. 
% UAI. 2016.

% Steps:
% 1. Set hyperparameters
% 2. Define truncated region and Gaussian density
% 3. Calculate EP-approximation
% 4. Run ESS given the EP approximation
% 5. Run Analytic ESS given the EP approximation
% 6. Plotting
% 7. Display statistics of interets


%% 1. Set hyperparameters

% Experiment parameters
dimension = 2;              % [2,10,50,100]

% MCMC parameters
number_chains = 4;          % Chains are run independently
number_samples = 100;       % Number samples per chain

% Indicators of which methods to run and whether to plot
RUN_EPESS = true;
RUN_ANALYTIC_EPESS = true;
RUN_PLOTS = false;


%% 2. Define truncated region and Gaussian density
% For details on how to define the truncated region, 
% refer to the code in the epmgp directory.
% Below a default TMG setup is used.

% Use default TMG
[mu, Sigma, chol_Sigma, C, lB, uB ] = default_TMG(dimension);

% Define log-likelihood
logLikelihood = @(x)( logPdfTmg( x, mu, chol_Sigma, C, lB, uB ));


%% 3. Calculate EP-approximation
% Either run epmgp for general polyhedrons or 
% axisepmgp for axis alligned boxes

% General polyhedrons
% [logZ, EP_mean , EP_covariance] = epmgp(mu,Sigma,C',lB,uB);

% Axis alligned boxes
[logZ, EP_mean , EP_covariance] = axisepmgp(mu,Sigma,lB,uB);

% Calculate mean and cholesky decomposition
EP_mean = EP_mean';
EP_chol = chol(EP_covariance);


%% 4. Run ESS given the EP approximation

if RUN_EPESS
    disp('Running EPESS')
    
    % Run the slice sampling code
    [ samples_epess, num_fn_evals ] = epess_sampler( number_samples, dimension, number_chains, logLikelihood, EP_mean, EP_chol);

    % Calculate effective sample size
    effective_sample_size = mpsrf(samples_epess);
    
    % Record performance
    % Stores the effective sample size, number of function evaluations 
    % and ratio of effective sample size /  number of function evaluations
    performance_epess = [effective_sample_size, num_fn_evals];
end


%% 5. Run Analytic ESS given the EP approximation
% Note that in this implementation there is only one slice per iteration (J=1)

if RUN_ANALYTIC_EPESS
    disp('Running Efficient EPESS')

    % Specify the Polyhedron
    F=[C; -C];
    g = [-lB;uB];

    % Calculate inverse of the EP covarianc matrix
    EP_cov_inv = inv(EP_covariance);
    
    % Run the analytic slice sampling code
    [ samples_analytic_epess, num_fn_evals ] = analytic_epess_sampler( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, F, g, EP_cov_inv);

    % Calculate effective sample size
    effective_sample_size = mpsrf(samples_analytic_epess);
    
    % Record performance
    % Stores the effective sample size, number of function evaluations 
    % and ratio of effective sample size /  number of function evaluations
    performance_analytic_epess = [effective_sample_size, num_fn_evals];

end


%% 6. Plotting
if RUN_PLOTS

    subplot(1,2,1);
    plot(samples_epess(:,1), samples_epess(:,2), 'x')
    axis([lB(1) , uB(1), lB(2), uB(2)])
    title('EPESS')

    subplot(1,2,2);
    plot(samples_analytic_epess(:,1), samples_analytic_epess(:,2), 'x')
    axis([lB(1) , uB(1), lB(2), uB(2)])
    title('Analytic EPESS')

%     % The code below plots the multivariate Gaussian density
%         hold on
%         ezcontour(@(x,y)(mvnpdf([x;y], EP_mean' , EP_covariance)) , [lB(1)-.5 , uB(1)+.5, lB(2)-1, uB(2)+1] , 200)
%         axis([lB(1)-.5 , uB(1)+.5, lB(2)-1, uB(2)+1])
%         hold off

end
    

%% 7. Display statistics of interets

disp(' ')

if RUN_EPESS
    disp('Results for EPESS')
    disp('----------------------------------------------')
    disp(['Effective sample size: ',  num2str(performance_epess(1))])
    disp(['Number of function evaluations: ',  num2str(performance_epess(2))])
    disp(['Effective sample size / Number of function evaluations: ',  num2str(performance_epess(1)/performance_epess(2))])
    disp(' ')
end

if RUN_ANALYTIC_EPESS
    disp('Results for Analytic EPESS')
    disp('----------------------------------------------')
    disp(['Effective sample size: ',  num2str(performance_analytic_epess(1))])
    disp(['Number of function evaluations: ',  num2str(performance_analytic_epess(2))])
    disp(['Effective sample size / Number of function evaluations: ',  num2str(performance_analytic_epess(1)/performance_analytic_epess(2))])
    disp(' ')
end


