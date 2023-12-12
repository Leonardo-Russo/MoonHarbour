function[aG] = MoonHarmPerts(MEE, MoonPPsECI, t, timespan, muM, deltaE, psiM, ...
    deltaM)

% Retrieve COE and MCI States
COE = MEE2COE(MEE'); 
X_MCI = COE2rvPCI(COE, muM)';

% Find Moon State in ECI from Interpolation of Ephemeris - PolyEval()
rM_ECI = [PolyEval(t, timespan, flip(MoonPPsECI(1).coefs, 2)); ...
    PolyEval(t, timespan, flip(MoonPPsECI(2).coefs, 2)); ...
    PolyEval(t, timespan, flip(MoonPPsECI(3).coefs, 2))];

% Compute Geographical Coordinates
GeoElem = Cart2Geo(X_MCI);
r = GeoElem(1);
LonAbs = GeoElem(2);
Lat = GeoElem(3);
zeta = GeoElem(6);

% Define ECI2MCI Rotation Matrices
R1deltaE = R1(-deltaE);
R3PsiM = R3(psiM);
R1deltaM = R1(deltaM);

% Compute r_hat
rE2M_MCI = (rM_ECI' * R1deltaE * R3PsiM' * R1deltaM')';
r_hat = rE2M_MCI / norm(rE2M_MCI);

% Compute Latitude and Longitude
Lat0 = asin(r_hat(3));
Lon0 = 2 * atan((r_hat(2) / cos(Lat0)) / (1 + r_hat(1) / cos(Lat0)));
LonGeo = LonAbs - Lon0;

% Compute Perturbing Accelerations in REN
MoonHarm("LPE", "REDUCED")
[ar1, aE1, aN1] = pertAccTessSect(r, Lat, LonGeo);
[ar2, aE2, aN2] = pertAccZon(r, Lat);

ar = ar1 + ar2;
aE = aE1 + aE2;
aN = aN1 + aN2;

% Convert Perturbing Accelerations into LVLH
R1zeta = R1(zeta);
aG = ([ar, aE, aN] * R1zeta')';

end