function failure_times = get_failure_times(t0, tf, TU)
% Description: this function create the initial and final time of the
% temporary engine failure between the interval t0 and tf.

% Uniform Distribution for Failure End Time
failure_end_lower = t0;
failure_end_upper = tf;
end_dist = makedist('uniform', 'Lower', failure_end_lower, 'Upper', failure_end_upper);
failure_end_time = random(end_dist, 1);

% Uniform Distribution for Failure Duration
failure_duration_lower = 0;
failure_duration_upper = 5*60/TU;       % 5 minutes
duration_dist = makedist('uniform', 'Lower', failure_duration_lower, 'Upper', failure_duration_upper);
failure_duration = random(duration_dist, 1);

% Make sure the initial failure time is always after t0
if failure_end_time - failure_duration < t0
    failure_start_time = t0;
else
    failure_start_time = failure_end_time - failure_duration;
end

failure_times = [failure_start_time; failure_end_time];


end