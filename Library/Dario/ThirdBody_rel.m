function[a34B_ORB] = ThirdBody_rel(EquinElem,COE_target, t, ppEarthMCI, ppSunMCI, muE, muS, ...
    timespan, muM)
% HELP
%
% This function evaluates the perturbation due to third and four body,
% projecting it on the LVLH reference frame centered on the virtual
% vehicle, used in relative motion.

%-------------------------- CLASSICAL ELEMENTS ----------------------------

ClassElem = Equin2Class(EquinElem); 
CartElem = Class2Cart(ClassElem, muM);

R3Omega = angleToR3(COE_target(4));
R1inc = angleToR1(COE_target(3));
R3theta = angleToR3(COE_target(5) + COE_target(6));

%------------------------------ EPHEMERIS ---------------------------------

rM2E_MCI = [PolyEval(t, timespan, flip(ppEarthMCI(1).coefs, 2)); ...
    PolyEval(t, timespan, flip(ppEarthMCI(2).coefs, 2)); ...
    PolyEval(t, timespan, flip(ppEarthMCI(3).coefs, 2))];
 
% rM2E_MCI = [PolyEval(t, timespan, ppEarthMCI(:,1:4)); ...
%      PolyEval(t, timespan, ppEarthMCI(:,5:8)); ...
%      PolyEval(t, timespan, ppEarthMCI(:,9:12))];

rM2S_MCI = [PolyEval(t, timespan, flip(ppSunMCI(1).coefs, 2)); ...
    PolyEval(t, timespan, flip(ppSunMCI(2).coefs, 2)); ...
    PolyEval(t, timespan, flip(ppSunMCI(3).coefs, 2))];

% rM2S_MCI = [PolyEval(t, timespan, ppSunMCI(:,1:4)); ...
%      PolyEval(t, timespan, ppSunMCI(:,5:8)); ...
%      PolyEval(t, timespan, ppSunMCI(:,9:12))];


%------------------------------- VECTORS ----------------------------------

rM2SC_MCI = CartElem(1:3);
rE2M_MCI = - rM2E_MCI;
rS2M_MCI = - rM2S_MCI;
rE2SC_MCI = rE2M_MCI + rM2SC_MCI;
rS2SC_MCI = rS2M_MCI + rM2SC_MCI;

%------------------------- 3/4 BODY ACCELERATION --------------------------

qE = (norm(rM2SC_MCI)^2 - 2 * dot(rM2SC_MCI, rM2E_MCI)) / norm(rM2E_MCI)^2;
qS = (norm(rM2SC_MCI)^2 - 2 * dot(rM2SC_MCI, rM2S_MCI)) / norm(rM2S_MCI)^2;

a3B = - muE / (norm(rM2E_MCI)^3 * (qE + 1)^(3/2)) * (rM2SC_MCI + rM2E_MCI ...
    * qE * (qE^2 + 3*qE + 3) / ((qE + 1)^(3/2) + 1));

a4B = - muS / (norm(rM2S_MCI)^3 * (qS + 1)^(3/2)) * (rM2SC_MCI + rM2S_MCI ...
    * qS * (qS^2 + 3*qS + 3) / ((qS + 1)^(3/2) + 1));

a34B_MCI = a3B + a4B;
a34B_ORB = (a34B_MCI' * R3Omega' * R1inc' * R3theta')';


% a34B_ORB = zeros(3,1);




