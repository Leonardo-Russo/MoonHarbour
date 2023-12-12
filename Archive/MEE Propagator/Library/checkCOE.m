function checkCOE(COE)
% Description: this function checks if each parameter in the COE vector is
% in the correct interval.

incl = COE(3);
Omega = COE(4);
omega = COE(5);
nu = COE(6);

errs = [];
n = 0;

if incl < 0 || incl > pi
    fprintf('incl is: %f\n', rad2deg(incl));
    n = n + 1;
end

if Omega < -pi || Omega > pi
    fprintf('Omega is: %f\n', rad2deg(Omega));
    n = n + 1;
end

if omega < -pi || omega > pi
    fprintf('omega is: %f\n', rad2deg(omega));
    n = n + 1;
end

if nu < -pi || nu > pi
    fprintf('nu is: %f\n', rad2deg(nu));
    n = n + 1;
end

if n > 0
    error('Conversion Error!')
end


end