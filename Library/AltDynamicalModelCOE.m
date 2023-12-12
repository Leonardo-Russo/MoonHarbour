function dCOE = AltDynamicalModelCOE(t, COE, EarthPPsMCI, SunPPsMCI, muE, ...
                                  muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf)

% Description: this function contains the Dynamical Model for the
% integration of the Modified Equinoctial Elements.
% 
% Inputs:
% t = epoch
% MEE = Equinoctial Elements State
% EarthPPsMCI = pp struct with the interpolation of the Earth State in MCI
% SunPPsMCI = pp struct with the interpolation of the Sun State in MCI
% muE = Earth's gravitational parameter in canonical units
% muS = Sun's gravitational parameter in canonical units
% time = reduced time span
% MoonPPsECI = pp struct with the interpolation of the Moon State in ECI
% deltaE = Earth's Ecliptic Obliquity
% psiM = Moon's Precession Angle
% deltaM = Moon's Equator Obliquity
% t0 = initial time
% tf = final time
% 
% Outputs:
% dMEE = derivatives of the MEE State


% Retrieve Global Variables
global muM Rm DU TU pbar

% Initialize Derivative Vector
dCOE = zeros(6, 1);


% Clock for the Integration
Day = 86400;  % seconds in a day
Hour = 3600;  % seconds in an hour
Min = 60;     % seconds in a minute

tDAY = floor((t - t0) * TU / Day);      % calculate the elapsed time components
tHR = floor(((t - t0) * TU - tDAY * Day) / Hour);
tMIN = floor(((t - t0) * TU - tDAY * Day - tHR * Hour) / Min);

timeStr = sprintf('Time Elapsed: %02d days, %02d hrs, %02d mins', tDAY, tHR, tMIN);     % create a string for the time

waitbarMessage = sprintf('Progress: %.2f%%\n%s', (t-t0)/(tf-t0)*100, timeStr);      % create the waitbar message including the time and progress percentage

waitbar((t-t0)/(tf-t0), pbar, waitbarMessage);      % update the waitbar


% Retrieve COE Values
a = COE(1);
e = COE(2);
incl = COE(3);
Omega = COE(4);
omega = COE(5);
nu = COE(6);

% Introduce Local Variables
theta = omega + nu;
p = a*(1-e^2);
X_MCI = COE2rvPCI(COE', muM)';
r_vect = X_MCI(1:3);
v_vect = X_MCI(4:6);

r = norm(r_vect);
h_vect = cross(r_vect, v_vect);
h = norm(h_vect);

MEE = COE2MEE(COE')';

% Compute Third and Fourth Body Accelerations - ThirdFourthBody()
a34B = ThirdFourthBody(MEE, t, EarthPPsMCI, SunPPsMCI, muE, muS);

% Compute Moon Harmonics Perturbation Acceleration
aG_M = MoonHarmPerts(MEE, MoonPPsECI, t, muM, deltaE, psiM, deltaM);

ap = a34B + aG_M;
fr = ap(1);
ft = ap(2);
fh = ap(3);

% Assign Values to State Derivative
dCOE(1) = 2*a^2/h * (e*fr*sin(nu) + ft*(e*cos(nu) + 1));
dCOE(2) = sqrt(p/muM)*fr*sin(nu) + sqrt(p/muM)*ft*(e+e*cos(nu)^2+2*cos(nu))/(1+e*cos(nu));
dCOE(3) = r*fh/h * cos(theta);
dCOE(4) = r*fh/h *sin(theta)/sin(incl);
dCOE(5) = -r*fh*sin(theta)*cos(incl)/(h*sin(incl)) - fr/e*cos(nu)*sqrt(p/muM) + ft*sqrt(p/muM)*sin(nu)*(e*cos(nu)+2)/(e*(1+e*cos(nu)));
dCOE(6) = sqrt(muM/p^3)*(1+e*cos(nu))^2 + fr/e*cos(nu)*sqrt(p/muM) - ft*sqrt(p/muM)*sin(nu)*(e*cos(nu)+2)/e*(1+e*cos(nu));

end


