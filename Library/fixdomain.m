function angle = fixdomain(angle, min, max)
% Description: this function helps fixing the domain for a trigonometric
% angle.

max_iter = 100;

for i = 1 : max_iter
    
    if angle < min
        angle = angle + 2*pi;
    end

    if angle > max
        angle = angle - 2*pi;
    end

    if angle >= min && angle <= max
        break
    end

end

end