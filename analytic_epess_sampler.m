function [ samples, num_fn_evals ] = analytic_epess_sampler( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, F, g, EP_cov_inv, initial_point) %axis_interval
% Generate samples using Analytic EPESS

    % Define log-likelihood shifted by the EP mean and covariance
    pseudoLogLikelihoodShifted = @(x)( logLikelihood(x+EP_mean) - logGaussPdfChol(x', zeros(dimension,1), EP_chol));

    % Initialize variables for MCMC
    samples = zeros(number_samples , dimension, number_chains);
    num_fn_evals = 1;
        
    % Run MCMC
    for chain_index = 1:number_chains

        % Initialize chain        
        if (nargin < 12) || isempty(initial_point)
            initial_point = 0;
        end
        samples(1,:,chain_index) = initial_point;
        cur_log_like = pseudoLogLikelihoodShifted(initial_point);
    
        % Run MCMC
        for sample_index = 2:number_samples

            % Generate samples
            [samples(sample_index,:,chain_index), cur_log_like, cur_num_fn_evals] = analytic_slice( samples(sample_index-1,:,chain_index), EP_chol, pseudoLogLikelihoodShifted, cur_log_like, F, g, EP_mean, dimension, EP_cov_inv);
            
            % Update number of function evaluations
            num_fn_evals = num_fn_evals + cur_num_fn_evals; % Make sure that adjust for dividing by dimension
        end
        
    end
    
    % Add the EP mean back to the samples
    samples = samples + repmat( EP_mean , number_samples, 1 , number_chains );
    
end
