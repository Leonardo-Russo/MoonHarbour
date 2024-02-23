function [crosses, indexes, directions] = detectCross(vect, threshold)
% HELP
%
% This function determines how many crosses of a given threshold are done 
% by the scalar values stored in "vect".
%
% Inputs:
% - vect : row or column vector 
% - threshold : scalar value whose crosses are evaluated
%
% Outputs:
% - crosses : number of threshold crosses detected
% - indexes : indexes of vect corresponding to the value BEFORE each cross
% - directions : this scalar indicates if the cross occurs from above(1) or
%               from below(-1)

    len = length(vect);
    
    crosses = 0;
    for i = 1 : len-1
        
        if vect(i)>threshold && vect(i+1)<threshold
            % cross from above 
            crosses = crosses + 1;

        elseif vect(i)<threshold && vect(i+1)>threshold
            % cross from below 
            crosses = crosses + 1;
        end

    end

    len2 = crosses;
    indexes = zeros(len2,1);
    directions = zeros(len2,1);
    k = 1;
    for i = 1 : len-1
        
        if vect(i)>threshold && vect(i+1)<threshold
            % cross from above 
            indexes(k) = i;
            directions(k) = 1;
            k = k + 1;

        elseif vect(i)<threshold && vect(i+1)>threshold
            % cross from below 
            indexes(k) = i;
            directions(k) = -1;
            k = k + 1;

        end

    end


end