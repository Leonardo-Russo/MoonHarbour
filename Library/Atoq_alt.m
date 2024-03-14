function [q0, qvec] = Atoq_alt(A)

q0 = 0;
qvec = zeros(3,1);

q0M = sqrt(1 + A(1,1) + A(2,2) + A(3,3))/2;
q1M = sqrt(1 + A(1,1) - A(2,2) - A(3,3))/2;
q2M = sqrt(1 - A(1,1) + A(2,2) - A(3,3))/2;
q3M = sqrt(1 - A(1,1) - A(2,2) + A(3,3))/2;

qMAX = max([q0M q1M q2M q3M]);

if q0M == qMAX
    
    q0 = q0M;
    qvec(1) = (A(2,3) - A(3,2))/(4*q0);
    qvec(2) = (A(3,1) - A(1,3))/(4*q0);
    qvec(3) = (A(1,2) - A(2,1))/(4*q0);
    
elseif q1M == qMAX
    
    qvec(1) = q1M;
    q0 = (A(2,3) - A(3,2))/(4*qvec(1));
    qvec(2) = (A(1,2) + A(2,1))/(4*qvec(1));
    qvec(3) = (A(3,1) + A(1,3))/(4*qvec(1));
    
elseif q2M == qMAX
    
    qvec(2) = q2M;
    q0 = (A(3,1) - A(1,3))/(4*qvec(2));
    qvec(1) = (A(1,2) + A(2,1))/(4*qvec(2));
    qvec(3) = (A(2,3) + A(3,2))/(4*qvec(2));
    
elseif q3M == qMAX
    
    qvec(3) = q3M;
    q0 = (A(1,2) - A(2,1))/(4*qvec(3));
    qvec(1) = (A(3,1) + A(1,3))/(4*qvec(3));
    qvec(2) = (A(2,3) + A(3,2))/(4*qvec(3));

end
   
        
        
        
        
        

% qvec(3) = sqrt(1 - A(1,1) - A(2,2) + A(3,3))/2;
% 
% q0 = (A(1,2) - A(2,1))/(4*qvec(3));
% qvec(1) = (A(3,1) + A(1,3))/(4*qvec(3));
% qvec(2) = (A(2,3) + A(3,2))/(4*qvec(3));

end
