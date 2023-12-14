function[ClassElem] = Equin2Class(EquinElem)

%------------------------- CLASSICAL ELEMENTS -----------------------------

p = EquinElem(1);
l = EquinElem(2);
m = EquinElem(3);
n = EquinElem(4);
s = EquinElem(5);
q = EquinElem(6);

%------------------------ EQUINOCTIAL ELEMENTS ----------------------------

e = sqrt(l^2 + m^2);
a = p/(1 - (l^2 + m^2));
i = 2 * atan(sqrt(n^2 + s^2));
Omega = 2 * atan(s/(n + sqrt(n^2 + s^2)));
omega = 2 * atan(m/(sqrt(l^2 + m^2) + l)) - Omega;
f = q - 2 * atan(m/(l + sqrt(l^2 + m^2)));

ClassElem = [a; e; i; Omega; omega; f];
