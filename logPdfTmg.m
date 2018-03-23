function [ logLikelihood ] = logPdfTmg( x, mu, chol_Sigma, C, lB, uB )
% Computes the log-likelihood of point x inside the polyhedron
% If the point is outside the polyhedron, it returns log(0)

% Check that all constraints are satisfied
n = length(lB);                 % Number of constraints
indicators = zeros(n,1);        % Indicators that constraints are satisfied
for j=1:n
    k = C(j,:)*x';
    indicators(j) = double(k >= lB(j) & uB(j)>= k);
end

% Calculate log-likelihood
if prod(indicators) == 0     % point is outside the polyhedron
    logLikelihood=log(0);
else                         % point is inside the polyhedron
    logLikelihood = logGaussPdfChol(x', mu, chol_Sigma) ;
end

end

