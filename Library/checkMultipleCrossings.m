function multipleCrossings = checkMultipleCrossings(distance, min_dist)
    % Initialize the crossing counter
    crossingCount = 0;
    
    if nargin < 2
        min_dist = 9.8;
    end

    % Initialize the state: 1 if above min_dist, 0 if below or equal to min_dist
    currentState = distance(1) > min_dist;

    % Loop through the distance array
    for i = 2:length(distance)
        % Check the current state
        if currentState && distance(i) <= min_dist
            % If it was above min_dist and now it's min_dist or below, increment the crossing count
            crossingCount = crossingCount + 1;
            % Change the state
            currentState = 0;
        elseif ~currentState && distance(i) > min_dist
            % If it was min_dist or below and now it's above min_dist, change the state
            currentState = 1;
        end
    end

    % Check if crossingCount is greater than 1
    multipleCrossings = crossingCount > 1;
end
