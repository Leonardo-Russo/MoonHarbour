function rM2DSG_SYN = MCI2SYN(X_MCI, time, MoonPPsECI, deltaE, psiM, deltaM)
% Description: this function find the position of the DSG in the Synodic
% Reference frame.
% 
% Inputs:
% X_MCI = State of the DSG in MCI reference frame
% time = time span vector
% MoonPPsECI = pp struct with the interpolation of the Moon State in ECI
% deltaE = Earth's Ecliptic Obliquity
% psiM = Moon's Precession Angle
% deltaM = Moon's Equator Obliquity
% 
% Outputs:
% rM2DSG_SYN = position vector of the DSG in SYN reference frame

% Initialize rM2DSG_SYN Matrix
rM2DSG_SYN = zeros(length(X_MCI(:,1)), 3);

% Compute rM2DSG_SYN at each epoch in time
for i = 1 : size(X_MCI, 1)
    
    % Retrieve Necessary States
    rDSG_MCI = X_MCI(i, 1:3)';
    
    % Find Moon State in ECI - ppsval()
    XM_ECI = ppsval(MoonPPsECI, time(i));

    % % Find Moon State in ECI - PolyEval()
    % rM_ECI = [PolyEval(t(i), time, flip(ppMoonECI(1).coefs, 2)); ...
    % PolyEval(t(i), time, flip(ppMoonECI(2).coefs, 2)); ...
    % PolyEval(t(i), time, flip(ppMoonECI(3).coefs, 2))];
    % vM_ECI = [PolyEval(t(i), time, flip(ppMoonECI(4).coefs, 2)); ...
    % PolyEval(t(i), time, flip(ppMoonECI(5).coefs, 2)); ...
    % PolyEval(t(i), time, flip(ppMoonECI(6).coefs, 2))];
    
    rM_ECI = XM_ECI(1:3);
    vM_ECI = XM_ECI(4:6);
    hM_ECI = cross(rM_ECI, vM_ECI);
    
    % Define ECI2MCI Rotation Matrices
    R1deltaE = R1(-deltaE);
    R3PsiM = R3(psiM);
    R1deltaM = R1(deltaM);
    
    rM2DSG_ECI = (rDSG_MCI' * R1deltaM * R3PsiM * R1deltaE')';
    
    % Determine SYN2ECI Rotation Matrix
    R_SYN2ECI = [(rM_ECI/norm(rM_ECI))'; ...                                % i_hat
                 cross(hM_ECI/norm(hM_ECI), rM_ECI/norm(rM_ECI))'; ...      % j_hat
                 (hM_ECI/norm(hM_ECI))'];                                   % k_hat

    % Convert DSG ECI State into Synodic State
    rM2DSG_SYN(i, :) = rM2DSG_ECI' * R_SYN2ECI';

end

end