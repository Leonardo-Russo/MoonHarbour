function[ClassElem] = Cart2Class(CartElem, muM)

x = CartElem(1);
y = CartElem(2);
z = CartElem(3);
vx = CartElem(4);
vy = CartElem(5);
vz = CartElem(6);

r = sqrt(x^2 + y^2 + z^2);
rvect = [x, y, z];
rver = rvect / r;

v = sqrt(vx^2 + vy^2 + vz^2);
vvect = [vx, vy, vz];
vver = vvect / v;

hvect = cross(rvect, vvect);
h = norm(hvect);
hver = hvect / h;

a = muM / (2 * muM/r - v^2);
e = sqrt(1 - h^2 / (muM * a));

p = a * (1 - e^2);
thetaver = cross(hver, rver);
vr = 1 / r * dot(vvect, rvect);

i = acos(hver(3));
Omega = atan2(hver(1) / sin(i),- hver(2) / sin(i));
thetastar = atan2(vr / e * sqrt(p / muM), 1 / e * (p / r - 1));
thetat = atan2(rver(3) / sin(i), thetaver(3) / sin(i));
omega = wrapToPi(thetat - thetastar);

ClassElem = [a; e; i; Omega; omega; thetastar];
