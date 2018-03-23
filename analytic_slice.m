function [xx, cur_log_like, number_fn_evaluations] = analytic_slice(xx, prior, log_like_fn, cur_log_like, F, g, EP_mean, dimension, EP_cov_inv)
% This code computes and samples from the acceptable slices of the ellipse 
% which lie within the truncated region

% Calculate the dimension
D = numel(xx);

% Get the momentum variable from the prior
if numel(prior) == D
    % User provided a prior sample:
    nu = reshape(prior, size(xx));
else
    % User specified Cholesky of prior covariance:
    if ~isequal(size(prior), [D D])
        error('Prior must be given by a D-element sample or DxD chol(Sigma)');
    end
    nu = reshape(prior'*randn(D, 1), size(xx));
end
hh = log(rand) + cur_log_like; % Log-likelihood threshold


%% Calculate where ellipse intersects with the polyhedron
[angle_polyhedron_intersection, number_fn_evaluations] = ellipse_polyhedron_intersection(xx, nu, F, g, EP_mean);

%% Calculate where the likelihood is > threshold
% Given by the roots of a quartic equation

% Defining quartic equation constants
mat = EP_cov_inv - eye(dimension);

a0 = -2*hh + nu*mat*nu';
a1 = -2*xx*EP_mean';
a2 = -2*nu*EP_mean';
a3 = 2*xx*mat*nu';
a4 = xx*mat*xx' - (a0+2*hh);

a = (a4^2) + (a3^2);
b = (2*a1*a4) + (2*a2*a3);
c = a1^2 + 2*a0*a4 - a3^2 +a2^2;
d = (2*a0*a1) - (2*a2*a3);
e = (a0^2) - (a2^2);

% Analytic Solution to the quartic equation:
p = (8*a*c - 3*(b^2))/(8*(a^2));
q = (b^3 - 4*a*b*c + 8*(a^2)*d)/(8*(a^3));

delta0 = (c^2) - 3*b*d + 12*a*e;
delta1 = 2*(c^3) - 9*b*c*d + 27*(b^2)*e + 27*a*(d^2) - 72*a*c*e;

temp = (delta1^2 - (4*(delta0^3)) )^0.5;
Q = ((delta1 + temp)/2)^(1/3);
S = 0.5* ( ( (Q + (delta0/Q))/(3*a) - (2/3)*p)^(0.5));

k1 = ((q/S) - 2*p - 4*(S^2))^0.5;
k2 = (-(q/S) - 2*p - 4*(S^2))^0.5;

x1 = -(b/(4*a)) - S + 0.5*k1;
x2 = -(b/(4*a)) - S - 0.5*k1;
x3 = -(b/(4*a)) + S + 0.5*k2;
x4 = -(b/(4*a)) + S - 0.5*k2;

%     Tested against Matlab solver -- the roots match       
%     syms y;
%     S = vpa(solve(y^4*a + y^3*b + y^2*c + y*d + e,y, 'MaxDegree',4))

% Extract the real roots
roots  = [x1, x2, x3, x4];
TOL=1e-10;
np=abs(imag(roots))<TOL;
roots = sort(real(roots(np)));
real_roots  = roots(logical((roots >= -1).*(roots <= 1)));
thresh_region = @(x)( a4*(cos(x)^2) + a3*cos(x)*sin(x) + a2*sin(x) + a1*cos(x) + a0);

% Find angles where likelihood equal to threshold
angle_likelihood_threshold = [0, 2*pi]; % Default if Likelihood > threshold everywhere
if numel(real_roots) > 0  % Likelihood < threshold somewhere

   unsorted_angles = [];
   for k=1:length(real_roots)
       angle = acos(real_roots(k));
       
       % Compute whether intersects at angle or 2*pi - angle
       potential_angles = [angle , 2*pi - angle];        
       [~, index] = min([ abs(thresh_region(potential_angles(1))), abs(thresh_region(potential_angles(2)))] );           
       
       % Add to unsorted angles
       unsorted_angles = [unsorted_angles, potential_angles(index)];
   end

   % Sort unsorted angles and store as angle_likelihood_threshold
   if numel(unsorted_angles) > 0
       unsorted_angles = [0, sort(unsorted_angles), 2*pi];
       angle_likelihood_threshold = []; 
       for i=1:(length(unsorted_angles)-1)
            % Check midway between consequitive angles whether above
            % threshold, if so add to angle_likelihood_threshold
            if thresh_region((unsorted_angles(i) + unsorted_angles(i+1))/2) > 0 
               angle_likelihood_threshold = [angle_likelihood_threshold, unsorted_angles(i), unsorted_angles(i+1)];
            end
       end
   end
end


%% Calculate slice
% Intersection of "where ellipse intersects with the polyhedron"
% and "where the likelihood is > threshold"

slice = range_intersection(angle_polyhedron_intersection, angle_likelihood_threshold);


%% Slice sampling loop
while true 
    phi = sample_from_slice(slice);
    xx_prop = xx*cos(phi) + nu*sin(phi);
    
    % Shouldn't have to check this, but no harm if we do
    cur_log_like = log_like_fn(xx_prop);
    if cur_log_like > hh
        % New point is on slice, ** EXIT LOOP **
        break;
    end

end
xx = xx_prop;


end
