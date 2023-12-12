function dCOE = DynamicalModelCOE(t, COE, EarthPPsMCI, SunPPsMCI, muE, ...
                                  muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf)

% Description: this function contains the Dynamical Model for the
% integration of the Modified Equinoctial Elements.
% 
% Inputs:
% t = epoch
% COE = Classical Orbital Elements
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
% dCOE = derivatives of the COE State


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

% Convert COE to MEE
MEE = COE2MEE(COE')';

% Compute Third and Fourth Body Accelerations - ThirdFourthBody()
a34B = ThirdFourthBody(MEE, t, EarthPPsMCI, SunPPsMCI, muE, muS);

% Compute Moon Harmonics Perturbation Acceleration
aG_M = MoonHarmPerts(MEE, MoonPPsECI, t, muM, deltaE, psiM, deltaM);

ap = a34B + aG_M;

% Compute State Derivatives
X_MCI = COE2rvPCI(COE', muM)';
dCOE = get_COEdots(COE, X_MCI, ap)';

end


