function a = Chebspace(leftbnd,rightbnd,ang0,angf,n)

% leftbnd estremo sx dell'intervallo
% rightbnd estremo dx dell'intervallo
% ang0 angolo corrispondente a leftbnd; standard è pi
% angf angolo corrispondente a rightbnd; standard è 0
% n numero punti (numero intervalli + 1)

% prima metà: ang0 = pi e angf = pi/2
% seconda metà: ang0 = pi/2 e angf = 0

chebychevgrid=@(leftbnd,rightbnd,n)leftbnd+((rightbnd-leftbnd)/(cos(angf)-cos(ang0)))*(cos(linspace(ang0,angf,n))-cos(ang0));
a = chebychevgrid(leftbnd,rightbnd,n);
% figure()
% plot(a, ones(size(a)),'-o')

end