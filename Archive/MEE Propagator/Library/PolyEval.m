function [x, xdot, xddot] = PolyEval(t, timespan, matr)

for index = 1 : length(timespan) - 1
    if (t >= timespan(index)) & (t <= timespan(index+1))
        
        x = matr(index,1) + matr(index,2) * (t-timespan(index)) + ...
            matr(index,3) * (t-timespan(index))^2 + matr(index,4) * ...
            (t-timespan(index))^3;
        
        xdot = matr(index,2) + 2 * matr(index,3) * (t-timespan(index)) + ...
            3 * matr(index,4) * (t-timespan(index))^2;
        
        xddot = 2 * matr(index,3) + 6 * matr(index,4) * (t-timespan(index));
        break;
    end
end



