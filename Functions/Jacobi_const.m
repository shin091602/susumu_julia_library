<<<<<<< HEAD
function C = Jacobi_const(x,mu)
% x  : non-dimensional position and velocity, x = [x y z vx vy vz]'
% mu : mass ratio of the primaries

	%File computes Jacobi Energy for a given state (x,v)

	%the distances
	r1 = sqrt((mu+x(1))^2+(x(2))^2+(x(3))^2);
	r2 = sqrt((1-mu-x(1))^2+(x(2))^2+(x(3))^2);

	%Compute the Jacobi Energy
	U = (x(1)^2 + x(2)^2)/2 + (1-mu)/r1 + mu/r2;
	v = sqrt(x(4)^2 + x(5)^2+x(6)^2);
	C = 2*U - v^2;
=======
function C = Jacobi_const(x,mu)
% x  : non-dimensional position and velocity, x = [x y z vx vy vz]'
% mu : mass ratio of the primaries

	%File computes Jacobi Energy for a given state (x,v)

	%the distances
	r1 = sqrt((mu+x(1))^2+(x(2))^2+(x(3))^2);
	r2 = sqrt((1-mu-x(1))^2+(x(2))^2+(x(3))^2);

	%Compute the Jacobi Energy
	U = (x(1)^2 + x(2)^2)/2 + (1-mu)/r1 + mu/r2;
	v = sqrt(x(4)^2 + x(5)^2+x(6)^2);
	C = 2*U - v^2;
>>>>>>> 748834e2ef04e45bd3c365491093c08a43a6640e
end