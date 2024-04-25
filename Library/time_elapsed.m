function t_split = time_elapsed(t2, t1, TU, format)
% Description: this function computes time elapsed between two timestamps
% provided in seconds and outputs the result in the provided format.

% Set time unit multiplier to 1 if not specified
if nargin < 3
    TU = 1;
end
if nargin < 4
    format = 'hhmmss';
end

% Define constants for time units
Day = 86400;  % seconds in a day
Hour = 3600;  % seconds in an hour
Min = 60;     % seconds in a minute

% Calculate the time elapsed
if format == "hhmmss"
    elapsedSeconds = (t2 - t1) * TU; % Calculate total elapsed seconds considering the time unit
    tHR = floor(elapsedSeconds / Hour); % Compute whole hours
    elapsedSeconds = elapsedSeconds - tHR * Hour; % Subtract hours part from total
    tMIN = floor(elapsedSeconds / Min); % Compute whole minutes
    tSEC = elapsedSeconds - tMIN * Min; % Remaining seconds after subtracting hours and minutes
    
    t_split = [tHR, tMIN, tSEC]; % Output result as [hours, minutes, seconds]
elseif format == "ddhhmmss"
    elapsedSeconds = (t2 - t1) * TU;
    tDAY = floor(elapsedSeconds / Day);
    elapsedSeconds = elapsedSeconds - tDAY * Day;
    tHR = floor(elapsedSeconds / Hour);
    elapsedSeconds = elapsedSeconds - tHR * Hour;
    tMIN = floor(elapsedSeconds / Min);
    tSEC = elapsedSeconds - tMIN * Min;  % calculate the remaining seconds
    
    t_split = [tDAY, tHR, tMIN, tSEC];
end

end
