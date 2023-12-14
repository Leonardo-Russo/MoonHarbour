function[R2] = angleToR2(angle)

R2 = [cos(angle) 0 -sin(angle); 0 1 0; sin(angle) 0 cos(angle)];

end