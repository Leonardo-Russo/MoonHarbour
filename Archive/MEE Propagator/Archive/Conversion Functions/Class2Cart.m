function[Cart] = Class2Cart(ClassElem, muM)

a = ClassElem(1);
e = ClassElem(2);
i = ClassElem(3);
Omega = ClassElem(4);
omega = ClassElem(5);
f = ClassElem(6);

r = a * (1 - e^2) / (1 + e * cos(f));
x = r * (cos(omega + f) * cos(Omega) -  cos(i) * sin(omega + f) * sin(Omega));
y = r * (cos(omega + f) * sin(Omega) +  cos(i) * sin(omega + f) * cos(Omega));
z = r * sin(i) * sin(omega + f);
vx = sqrt(muM / (a * (1 - e^2))) * (- cos(Omega) * (sin(omega + f) + ...
    e * sin(omega)) - sin(Omega) * cos(i) * (cos(omega + f) + e * cos(omega)));
vy = sqrt(muM / (a * (1 - e^2))) * (cos(Omega) * cos(i) * (cos(omega + f) + ...
    e * cos(omega)) - sin(Omega) * (sin(omega + f) + e * sin(omega)));
vz = sqrt(muM / (a * (1 - e^2))) * sin(i) * (cos(omega + f) + e * cos(omega));


Cart = [x; y; z; vx; vy; vz];
