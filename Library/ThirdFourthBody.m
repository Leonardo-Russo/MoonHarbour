function a34B_LVLH = ThirdFourthBody(MEE, t, EarthPPsMCI, SunPPsMCI, muE, muS)

% Recall Global Variables
global muM Rm TU DU

% Retrieve COE and MCI States
COE = MEE2COE(MEE'); 
X_MCI = COE2rvPCI(COE, muM)';
rSC_MCI = X_MCI(1:3);

% Define MCI2LVLH Rotation Matrices
R3Omega = R3(COE(4));
R1inc = R1(COE(3));
R3theta = R3(COE(5) + COE(6));

% Find Planet States from Interpolation of Ephemeris - ppsval()
XE_MCI = ppsval(EarthPPsMCI, t);
rE_MCI = XE_MCI(1:3);

XS_MCI = ppsval(SunPPsMCI, t);
rS_MCI = XS_MCI(1:3);

% Evaluate III and IV Body Perturbing Accelerations
qE = (norm(rSC_MCI)^2 - 2 * dot(rSC_MCI, rE_MCI)) / norm(rE_MCI)^2;
qS = (norm(rSC_MCI)^2 - 2 * dot(rSC_MCI, rS_MCI)) / norm(rS_MCI)^2;

a3B = - muE / (norm(rE_MCI)^3 * (qE + 1)^(3/2)) * (rSC_MCI + rE_MCI ...
    * qE * (qE^2 + 3*qE + 3) / ((qE + 1)^(3/2) + 1));

a4B = - muS / (norm(rS_MCI)^3 * (qS + 1)^(3/2)) * (rSC_MCI + rS_MCI ...
    * qS * (qS^2 + 3*qS + 3) / ((qS + 1)^(3/2) + 1));

a34B_MCI = a3B + a4B;
a34B_LVLH = (a34B_MCI' * R3Omega' * R1inc' * R3theta')';


end