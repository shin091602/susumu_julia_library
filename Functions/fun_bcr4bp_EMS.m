<<<<<<< HEAD
function  dx = fun_bcr4bp_EMS(t,x,mu,mS,aS,wS,thetaS0)
% Earth - Moon rotating frame
%       t : non-dimensional time
%       x : non-dimensional position and velocity, x = [x y z vx vy vz]'
%      mu : mass ratio of Moon, mu = m2 / (m1 + m2)
%      mS : nondimensional mass of Sun
%      aS : nondimensional Sun orbit radius
%      wS : nondimensional Sun angular velocity
% thetaS0 : the initial value of Sun angle


%the distances
     r1        = ((x(1)+mu)^2 + x(2)^2 + x(3)^2)^0.5;
     r2        = ((x(1)-1+mu)^2 + x(2)^2 + x(3)^2)^0.5;
     Xs        = aS*cos(thetaS0+wS*t);
     Ys        = aS*sin(thetaS0+wS*t);
     Zs        = 0;
     r3        = ((x(1)-Xs)^2 + (x(2)-Ys)^2 + (x(3)-Zs)^2)^0.5;

     dx        = [x(4); 
                  x(5);
                  x(6);
                2*x(5) + x(1) - (1-mu)*(x(1)+mu)/r1^3 - mu*(x(1)-1+mu)/r2^3 - mS*(x(1)-Xs)/r3^3 - mS/aS^3*Xs;         
               -2*x(4) + x(2) - (1-mu)*x(2)/r1^3      - mu*x(2)/r2^3        - mS*(x(2)-Ys)/r3^3 - mS/aS^3*Ys;         
                              - (1-mu)*x(3)/r1^3      - mu*x(3)/r2^3        - mS*(x(3)-Zs)/r3^3 - mS/aS^3*Zs];
=======
function  dx = fun_bcr4bp_EMS(t,x,mu,mS,aS,wS,thetaS0)
% Earth - Moon rotating frame
%       t : non-dimensional time
%       x : non-dimensional position and velocity, x = [x y z vx vy vz]'
%      mu : mass ratio of Moon, mu = m2 / (m1 + m2)
%      mS : nondimensional mass of Sun
%      aS : nondimensional Sun orbit radius
%      wS : nondimensional Sun angular velocity
% thetaS0 : the initial value of Sun angle


%the distances
     r1        = ((x(1)+mu)^2 + x(2)^2 + x(3)^2)^0.5;
     r2        = ((x(1)-1+mu)^2 + x(2)^2 + x(3)^2)^0.5;
     Xs        = aS*cos(thetaS0+wS*t);
     Ys        = aS*sin(thetaS0+wS*t);
     Zs        = 0;
     r3        = ((x(1)-Xs)^2 + (x(2)-Ys)^2 + (x(3)-Zs)^2)^0.5;

     dx        = [x(4); 
                  x(5);
                  x(6);
                2*x(5) + x(1) - (1-mu)*(x(1)+mu)/r1^3 - mu*(x(1)-1+mu)/r2^3 - mS*(x(1)-Xs)/r3^3 - mS/aS^3*Xs;         
               -2*x(4) + x(2) - (1-mu)*x(2)/r1^3      - mu*x(2)/r2^3        - mS*(x(2)-Ys)/r3^3 - mS/aS^3*Ys;         
                              - (1-mu)*x(3)/r1^3      - mu*x(3)/r2^3        - mS*(x(3)-Zs)/r3^3 - mS/aS^3*Zs];
>>>>>>> 748834e2ef04e45bd3c365491093c08a43a6640e
end 