function RHO_MCI = rhoLVLH2MCI(RHO_LVLH, Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM)
% Description: this function converts vectors in the MCI reference frame
% into the LVLH rotating reference frame centered on the Target trajectory.

% Recall Global Variables
global muM Rm DU TU

% Retrieve Data from Input
rho_LVLH = RHO_LVLH(1:3);
rhodot_LVLH = RHO_LVLH(4:6);

% Compute MCI2LVLH Rotation Matrix and its Derivative
[R_LVLH2MCI, Rdot_LVLH2MCI] = get_rotLVLH2MCI(Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);

% Rotate both Position and Velocity
rho_MCI = R_LVLH2MCI*rho_LVLH;
rhodot_MCI = R_LVLH2MCI*rhodot_LVLH + Rdot_LVLH2MCI*rho_LVLH;

% Assemble the MCI State
RHO_MCI = [rho_MCI; rhodot_MCI];

end