function[R1] = angleToR1(angle)

R1 = [1 0 0; 0 cos(angle) sin(angle); 0 -sin(angle) cos(angle)];

end