function [q0,qvect]=matrixToQuat(matrr)

global signvect;

qsqrd=[1,1,-1,-1;1,-1,1,-1;1,-1,-1,1;1,1,1,1]^(-1)*[matrr(1,1);matrr(2,2);matrr(3,3);1];

[maxq,indmax]=max(qsqrd);

% matrix=[q0^2+q1^2-q2^2-q3^2,2*(q1*q2+q0*q3),2*(q1*q3-q0*q2);...
%     2*(q1*q2-q0*q3),q0^2+q2^2-q1^2-q3^2,2*(q3*q2+q0*q1);...
%     2*(q1*q3+q0*q2),2*(q3*q2-q0*q1),q0^2+q3^2-q2^2-q1^2];

if indmax==1
    q0=signvect(1)*sqrt(maxq);
    q1=(matrr(2,3)-matrr(3,2))/(4*q0);
    q2=(matrr(3,1)-matrr(1,3))/(4*q0);
    q3=(matrr(1,2)-matrr(2,1))/(4*q0);
    signvect(2)=sign(q1);
    signvect(3)=sign(q2);
    signvect(4)=sign(q3);
elseif indmax==2
    q1=signvect(2)*sqrt(maxq);
    q0=(matrr(2,3)-matrr(3,2))/(4*q1);
    q2=(matrr(1,2)+matrr(2,1))/(4*q1);
    q3=(matrr(1,3)+matrr(3,1))/(4*q1);
    signvect(1)=sign(q0);
    signvect(3)=sign(q2);
    signvect(4)=sign(q3);
elseif indmax==3
    q2=signvect(3)*sqrt(maxq);
    q0=(matrr(3,1)-matrr(1,3))/(4*q2);
    q1=(matrr(1,2)-matrr(2,1))/(4*q2);
    q3=(matrr(2,3)+matrr(3,2))/(4*q2);
    signvect(1)=sign(q0);
    signvect(2)=sign(q1);
    signvect(4)=sign(q3);
elseif indmax==4
    q3=signvect(4)*sqrt(maxq);
    q0=(matrr(1,2)-matrr(2,1))/(4*q3);
    q1=(matrr(1,3)+matrr(3,1))/(4*q3);
    q2=(matrr(2,3)+matrr(3,2))/(4*q3);
    signvect(1)=sign(q0);
    signvect(2)=sign(q1);
    signvect(3)=sign(q2);
end

qvect=[q1;q2;q3];

% matrix=[q0^2+q1^2-q2^2-q3^2,2*(q1*q2+q0*q3),2*(q1*q3-q0*q2);...
%     2*(q1*q2-q0*q3),q0^2+q2^2-q1^2-q3^2,2*(q3*q2+q0*q1);...
%     2*(q1*q3+q0*q2),2*(q3*q2-q0*q1),q0^2+q3^2-q2^2-q1^2];

