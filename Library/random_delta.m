function deltaX = random_delta(dx_min, dx_max, dv_min, dv_max)
% Description: this function creates a random deltaX for the initial state
% while having it constrained by min and max values for position and
% velocity.

% Generate the random deltaX vector
deltaX = [dx_min + (dx_max - dx_min) * rand(3,1);...
          dv_min + (dv_max - dv_min) * rand(3,1)];

end