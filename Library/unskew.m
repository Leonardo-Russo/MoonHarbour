function vect = unskew(vect_skewed)
% Description: this function retrieves a vectors from its skewed
% representation.

if size(vect_skewed, 1) ~= 3 || size(vect_skewed, 2) ~= 3
    error('Wrong dimensions in unskew input!')
end

vect = zeros(3, 1);

vect(1) = -vect_skewed(2, 3);
vect(2) = vect_skewed(1, 3);
vect(3) = -vect_skewed(1, 2);


end