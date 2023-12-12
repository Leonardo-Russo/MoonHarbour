function Rdot = R1dot(angle, unit)

if nargin < 2
    unit = "rad";
end

if unit == "deg"
    angle = deg2rad(angle);
end

if unit ~= "deg" && unit ~= "rad"
    error('Wrong Unit provided as input.')
end

% Compute Rotation Matrix
Rdot = [0          0               0;
        0          -sin(angle)     cos(angle);
        0          -cos(angle)     -sin(angle)];


end