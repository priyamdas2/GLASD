function [lb, ub] = angular_bounds(M)
% M: dimension of the correlation matrix
% returns lb and ub vectors matching the size of Theta in Corr2Theta

num_angles = M*(M-1)/2;
lb = zeros(num_angles, 1);
ub = zeros(num_angles, 1);

idx = 1;

for m = 2:M
    num_angles_in_row = m - 1;

    if m == 2
        lb(idx) = -pi/2;
        ub(idx) =  pi/2;
    else
        lb(idx) = 0;             % First angle: [0, pi/2]
        ub(idx) = pi/2;

        if num_angles_in_row > 2
            % Middle angles: [0, pi]
            lb(idx+1:idx+num_angles_in_row-2) = 0;
            ub(idx+1:idx+num_angles_in_row-2) = pi;
        end

        % Last angle: [0, 2*pi)
        lb(idx + num_angles_in_row - 1) = 0;
        ub(idx + num_angles_in_row - 1) = 2*pi;
    end

    idx = idx + num_angles_in_row;
end
end