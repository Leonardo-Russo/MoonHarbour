function dMEE = DynamicalModelMEE(t, MEE, EarthPPsMCI, SunPPsMCI, muE, ...
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

% Retrieve MEE State Variables
x1 = MEE(1);
x2 = MEE(2);
x3 = MEE(3);
x4 = MEE(4);
x5 = MEE(5);
x6 = MEE(6);

% Retrieve Global Variables
global muM Rm DU TU pbar opt

if opt.show_progress
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
end

% Compute Local Variables
eta = get_eta(MEE);
G = get_G(MEE, muM, eta);

% % Compute Third and Fourth Body Accelerations - ThirdBody()
% a3B_E = ThirdBody(MEE, t, EarthPPsMCI, muE, time);
% a3B_S = ThirdBody(MEE, t, SunPPsMCI, muS, time);
% a34B = a3B_E + a3B_S;

% Compute Third and Fourth Body Accelerations - ThirdFourthBody()
a34B = ThirdFourthBody(MEE, t, EarthPPsMCI, SunPPsMCI, muE, muS);

% Compute Moon Harmonics Perturbation Acceleration
aG_M = MoonHarmPerts(MEE, MoonPPsECI, t, muM, deltaE, psiM, deltaM);

ap = a34B + aG_M;

% Assign State Derivatives
dMEE = G*ap;
dMEE(6) = dMEE(6) + sqrt(muM/x1^3)*eta^2;

end


