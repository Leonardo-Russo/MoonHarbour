function a3B_LVLH = ThirdBody(MEE, t, ThirdPPsMCI, mu3, time)

% Recall Global Variables
global muM Rm TU DU

% Retrieve COE and MCI States
COE = MEE2COE(MEE); 
X_MCI = COE2rvPCI(COE, muM)';
rSC_MCI = X_MCI(1:3);

% Define MCI2LVLH Rotation Matrices
R3Omega = R3(COE(4));
R1inc = R1(COE(3));
R3theta = R3(COE(5) + COE(6));

% % Find Third Body State from Interpolation of Ephemeris - ppsval()
% X3_MCI = ppsval(ThirdPPsMCI, t);
% r3_MCI = X3_MCI(1:3);

% Find Third Body State from Interpolation of Ephemeris - PolyEval()
r3_MCI = [PolyEval(t, time, flip(ThirdPPsMCI(1).coefs, 2)); ...
    PolyEval(t, time, flip(ThirdPPsMCI(2).coefs, 2)); ...
    PolyEval(t, time, flip(ThirdPPsMCI(3).coefs, 2))];

% Evaluate III and IV Body Perturbing Accelerations in MCI
q = (norm(rSC_MCI)^2 - 2 * dot(rSC_MCI, r3_MCI)) / norm(r3_MCI)^2;

a3B_MCI = - mu3 / (norm(r3_MCI)^3 * (q + 1)^(3/2)) * (rSC_MCI + r3_MCI ...
    * q * (q^2 + 3*q + 3) / ((q + 1)^(3/2) + 1));

% Transform Perturbation into LVLH Frame
a3B_LVLH = (a3B_MCI' * R3Omega' * R1inc' * R3theta')';

end