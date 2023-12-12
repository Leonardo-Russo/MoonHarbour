function[GeoElem] = Cart2Geo(CartElem)

%CARTESIAN ELEMENTS
x = CartElem(1);
y = CartElem(2);
z = CartElem(3);
vx = CartElem(4);
vy = CartElem(5);
vz = CartElem(6);

r = [x; y; z];
v = [vx; vy; vz];

%GEOGRAPHICAL ELEMENTS
gamma = asin((dot(v, r)/norm(r)) / norm(v));
rver = r /  norm(r);
h = cross(r, v) / norm(cross(r, v));
theta = cross(h, rver);
Lat = asin(rver(3));
LonAbs = 2 * atan((rver(2) / cos(Lat)) / (1 + rver(1) / cos(Lat)));
zeta = 2 * atan((theta(3) / cos(Lat)) / (1 + h(3) / cos(Lat)));

GeoElem = [norm(r); LonAbs; Lat; gamma; norm(v); zeta];
