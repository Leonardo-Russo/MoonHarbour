function[aG] = GravPert(state, ppMoonECI, t, timespan, muM, deltaE, psiM, ...
    deltaM)

%CARTESIAN STATE IN MCI
CartState = Class2Cart(Equin2Class(state(1:6)), muM);

%MOON STATE IN ECI
rE2M_ECI = [PolyEval(t, timespan, flip(ppMoonECI(1).coefs, 2)); ...
    PolyEval(t, timespan, flip(ppMoonECI(2).coefs, 2)); ...
    PolyEval(t, timespan, flip(ppMoonECI(3).coefs, 2))];
% rE2M_ECI = [PolyEval(t, timespan, ppMoonECI(:,1:4)); ...
%     PolyEval(t, timespan, ppMoonECI(:,5:8)); ...
%     PolyEval(t, timespan, ppMoonECI(:,9:12))];

%GEOGRAPHICAL COORDINATES
GeoElem = Cart2Geo(CartState);
r = GeoElem(1);
LonAbs = GeoElem(2);
Lat = GeoElem(3);
zeta = GeoElem(6);

%REFERENCE MERIDIAN
R1deltaE = angleToR1(-deltaE);
R3PsiM = angleToR3(psiM);
R1deltaM = angleToR1(deltaM);
rE2M_MCI = (rE2M_ECI' * R1deltaE * R3PsiM' * R1deltaM')' ...
    / norm((rE2M_ECI' * R1deltaE * R3PsiM' * R1deltaM')');
Lat0 = asin(rE2M_MCI(3));
Lon0 = 2 * atan((rE2M_MCI(1) / cos(Lat0)) / (1 + rE2M_MCI(2) / cos(Lat0)));
LonGeo = LonAbs - Lon0;

%PERTURBING ACCELERATIONS IN (R, E, N)
MoonHarm
[ar1, aE1, aN1] = pertAccTessSect(r, Lat, LonGeo);
[ar2, aE2, aN2] = pertAccZon(r, Lat);

ar = ar1 + ar2;
aE = aE1 + aE2;
aN = aN1 + aN2;

%PERTURBING ACCELERATIONS IN (R, T, H)
R1zeta = angleToR1(zeta);
aG = ([ar, aE, aN] * R1zeta')';

% aG = zeros(3,1);
