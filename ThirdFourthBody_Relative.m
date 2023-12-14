function[a34B_ORB] = ThirdFourthBody_Relative(EquinElem,COE_target, t, ppEarthMCI, ppSunMCI, muE, muS)
% HELP
%
% This function evaluates the perturbation due to third and four body,
% projecting it on the LVLH reference frame centered on the virtual
% vehicle, used in relative motion.

global muM

%-------------------------- CLASSICAL ELEMENTS ----------------------------

ClassElem = Equin2Class(EquinElem); 
CartElem = Class2Cart(ClassElem, muM);

R3Omega = angleToR3(COE_target(4));
R1inc = angleToR1(COE_target(3));
R3theta = angleToR3(COE_target(5) + COE_target(6));

%------------------------------ EPHEMERIS ---------------------------------

% Find Planet States from Interpolation of Ephemeris - ppsval()
XE_MCI = ppsval(ppEarthMCI, t);
rM2E_MCI = XE_MCI(1:3);

XS_MCI = ppsval(ppSunMCI, t);
rM2S_MCI = XS_MCI(1:3);

%------------------------------- VECTORS ----------------------------------

rM2SC_MCI = CartElem(1:3);

%------------------------- 3/4 BODY ACCELERATION --------------------------

qE = (norm(rM2SC_MCI)^2 - 2 * dot(rM2SC_MCI, rM2E_MCI)) / norm(rM2E_MCI)^2;
qS = (norm(rM2SC_MCI)^2 - 2 * dot(rM2SC_MCI, rM2S_MCI)) / norm(rM2S_MCI)^2;

a3B = - muE / (norm(rM2E_MCI)^3 * (qE + 1)^(3/2)) * (rM2SC_MCI + rM2E_MCI ...
    * qE * (qE^2 + 3*qE + 3) / ((qE + 1)^(3/2) + 1));

a4B = - muS / (norm(rM2S_MCI)^3 * (qS + 1)^(3/2)) * (rM2SC_MCI + rM2S_MCI ...
    * qS * (qS^2 + 3*qS + 3) / ((qS + 1)^(3/2) + 1));

a34B_MCI = a3B + a4B;
a34B_ORB = (a34B_MCI' * R3Omega' * R1inc' * R3theta')';


end