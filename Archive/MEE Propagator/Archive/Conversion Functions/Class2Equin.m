function[EquinElem] = Class2Equin(ClassElem)

%------------------------- CLASSICAL ELEMENTS -----------------------------

a = ClassElem(1);
e = ClassElem(2);
i = ClassElem(3);
Omega = ClassElem(4);
omega = ClassElem(5);
f = ClassElem(6);

%------------------------ EQUINOCTIAL ELEMENTS ----------------------------

x1 = a * (1 - e^2);
x2 = e * cos(Omega + omega);
x3 = e * sin(Omega + omega);
x4 = sin(i) / (1 + cos(i)) * cos(Omega);
x5 = sin(i) / (1 + cos(i)) * sin(Omega);
x6 = Omega + omega + f;

EquinElem = [x1; x2; x3; x4; x5; x6];
