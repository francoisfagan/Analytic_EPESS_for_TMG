function [ angle_slice, fn_eval] = ellipse_polyhedron_intersection( curr_point, nu, F, g, EP_mean)
% Compute where ellipse intersects with the polyhedron

% Project the ellipse onto each of the constraints
F_nu = F*nu';
F_curr_point = F*curr_point';
        
U = sqrt(F_nu.^2 + F_curr_point.^2); % Radius of 1-dimensional ellipse
phi = atan2(-F_nu,F_curr_point); % Initial angle of 1-dimensional ellipse

g = g + F*EP_mean'; % Projected constraints

% Check if the projected ellipse intersects with the constraints
intersect_indicator = abs(g./U)<1;

if any(intersect_indicator) % At least one constraint intersected
    
    % Extract initial angles of all constraints that intersect
    phn=phi(intersect_indicator);
    
    % Calculate the angles at which the ellipse crosses the constraints
    a1=-phn + acos(-g(intersect_indicator)./U(intersect_indicator));
    a2=-phn + 2*pi - acos(-g(intersect_indicator)./U(intersect_indicator));

    % Combine the angles into intervals where inside the polyhedron
    angle_slice = [0, a1(1), a2(1) ,2*pi];
    if length(a1) >1
        for i=2:length(a1)
            new_range = [0, a1(i), a2(i) ,2*pi];
            angle_slice = range_intersection(angle_slice, new_range); 
        end

    end 

else % No constraints intersected
    angle_slice=[0, 2*pi];
end

fn_eval = 1;

end

