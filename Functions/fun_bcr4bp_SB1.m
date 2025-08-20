<<<<<<< HEAD
function  dx = fun_bcr4bp_SB1(t,x,muSB1,muEM,aS,aEM,wM,thetaM0)
% Sun - B1 rotating frame
%       t : non-dimensional time
%       x : non-dimensional position and velocity, x = [x y z vx vy vz]'
%   muSB1 : mass ratio in Sun - B1 frame, muSB1 = (m1 + m2) / (m1 + m2 + ms)
%    muEM : mass ratio in Earth - Moon frame, muEM = m2 / (m1 + m2)
%      aS : nondimensional Sun orbit radius
%     aEM : nondimensional distance between Earth and Moon
%      wM : nondimensional Moon angular velocity
% thetaM0 : the initial value of Moon angle


%the distances
    thata   = thetaM0 + wM*t;
    rs      = ((x(1)+muSB1)^2 + x(2)^2 + x(3)^2)^0.5;
    rE      = [(-muEM)*cos(thata),(-muEM)*sin(thata),0].*(aEM/aS)+[1-muSB1,0,0];
    rM      = [(1-muEM)*cos(thata),(1-muEM)*sin(thata),0].*(aEM/aS)+[1-muSB1,0,0];
    r1      = norm([x(1)-rE(1),x(2)-rE(2),x(3)-rE(3)]);
    r2      = norm([x(1)-rM(1),x(2)-rM(2),x(3)-rM(3)]);


    dx      = [x(4); 
               x(5);
               x(6);
               2*x(5) + x(1) - (1-muSB1)*(x(1)+muSB1)/rs^3 - muSB1*(1-muEM)*(x(1)-rE(1))/r1^3 - muSB1*muEM*(x(1)-rM(1))/r2^3;         
              -2*x(4) + x(2) - (1-muSB1)*x(2)/rs^3         - muSB1*(1-muEM)*(x(2)-rE(2))/r1^3 - muSB1*muEM*(x(2)-rM(2))/r2^3;         
                             - (1-muSB1)*x(3)/rs^3         - muSB1*(1-muEM)*(x(3)-rE(3))/r1^3 - muSB1*muEM*(x(3)-rM(3))/r2^3]; 
=======
function  dx = fun_bcr4bp_SB1(t,x,muSB1,muEM,aS,aEM,wM,thetaM0)
% Sun - B1 rotating frame
%       t : non-dimensional time
%       x : non-dimensional position and velocity, x = [x y z vx vy vz]'
%   muSB1 : mass ratio in Sun - B1 frame, muSB1 = (m1 + m2) / (m1 + m2 + ms)
%    muEM : mass ratio in Earth - Moon frame, muEM = m2 / (m1 + m2)
%      aS : nondimensional Sun orbit radius
%     aEM : nondimensional distance between Earth and Moon
%      wM : nondimensional Moon angular velocity
% thetaM0 : the initial value of Moon angle


%the distances
    thata   = thetaM0 + wM*t;
    rs      = ((x(1)+muSB1)^2 + x(2)^2 + x(3)^2)^0.5;
    rE      = [(-muEM)*cos(thata),(-muEM)*sin(thata),0].*(aEM/aS)+[1-muSB1,0,0];
    rM      = [(1-muEM)*cos(thata),(1-muEM)*sin(thata),0].*(aEM/aS)+[1-muSB1,0,0];
    r1      = norm([x(1)-rE(1),x(2)-rE(2),x(3)-rE(3)]);
    r2      = norm([x(1)-rM(1),x(2)-rM(2),x(3)-rM(3)]);


    dx      = [x(4); 
               x(5);
               x(6);
               2*x(5) + x(1) - (1-muSB1)*(x(1)+muSB1)/rs^3 - muSB1*(1-muEM)*(x(1)-rE(1))/r1^3 - muSB1*muEM*(x(1)-rM(1))/r2^3;         
              -2*x(4) + x(2) - (1-muSB1)*x(2)/rs^3         - muSB1*(1-muEM)*(x(2)-rE(2))/r1^3 - muSB1*muEM*(x(2)-rM(2))/r2^3;         
                             - (1-muSB1)*x(3)/rs^3         - muSB1*(1-muEM)*(x(3)-rE(3))/r1^3 - muSB1*muEM*(x(3)-rM(3))/r2^3]; 
>>>>>>> 748834e2ef04e45bd3c365491093c08a43a6640e
end