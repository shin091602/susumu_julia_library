<<<<<<< HEAD
function [N, R] = fun_null_cr3bp(x0, t0, mu, options)
% x0      : initial guess of x0
% t0      : initial guess of period
% mu      : mass ratio of the primaries
% options : options for ode

  X0 = [x0(1), x0(2), x0(3), x0(4), x0(5), x0(6), reshape(eye(6), 1, [])];
  tspan = [0 t0];
  [~, Y] = ode113(@(t, x) fun_stm_cr3bp(t, x, mu), tspan, X0, options);

  X = Y(end, 1:6);
  Phi = reshape(Y(end, 7:end), 6, 6);
  f_x = fun_cr3bp([], X, mu);
  DF = [Phi(2, 1), Phi(2, 3), Phi(2, 5), f_x(2);
        Phi(4, 1), Phi(4, 3), Phi(4, 5), f_x(4);
        Phi(6, 1), Phi(6, 3), Phi(6, 5), f_x(6)];
  R = rank(DF);
  N = null(DF);   % Jacobian matlix
end
=======
function [N, R] = fun_null_cr3bp(x0, t0, mu, options)
% x0      : initial guess of x0
% t0      : initial guess of period
% mu      : mass ratio of the primaries
% options : options for ode

  X0 = [x0(1), x0(2), x0(3), x0(4), x0(5), x0(6), reshape(eye(6), 1, [])];
  tspan = [0 t0];
  [~, Y] = ode113(@(t, x) fun_stm_cr3bp(t, x, mu), tspan, X0, options);

  X = Y(end, 1:6);
  Phi = reshape(Y(end, 7:end), 6, 6);
  f_x = fun_cr3bp([], X, mu);
  DF = [Phi(2, 1), Phi(2, 3), Phi(2, 5), f_x(2);
        Phi(4, 1), Phi(4, 3), Phi(4, 5), f_x(4);
        Phi(6, 1), Phi(6, 3), Phi(6, 5), f_x(6)];
  R = rank(DF);
  N = null(DF);   % Jacobian matlix
end
>>>>>>> 748834e2ef04e45bd3c365491093c08a43a6640e
