<<<<<<< HEAD
% differential correction for halo family in CR3BP by PAC
function [x02_ast, t02_ast, C, G] = fun_differential_correction_cr3bp_PAC(x01_ast, t01_ast, x02, t02, scale, delta, mu, options)
% x01_ast : previously-converged x0 ( (i)th orbit )
% t01_ast : previously-converged period ( (i)th orbit )
% x02     : initial guess of x0 ( (i+1)th orbit )
% t02     : initial guess of period ( (i+1)th orbit )
% scale   : scaling factor of delta
% delta   : null space of the previously-converged Jacobian matrix
% mu      : mass ratio of the primaries
% options : options for ode

  X0 = [x02(1), x02(2), x02(3), x02(4), x02(5), x02(6), reshape(eye(6), 1, [])];
  tspan = [0 t02];
  [~, Y] = ode113(@(t, x) fun_stm_cr3bp(t, x, mu), tspan, X0, options);

  X = Y(end, 1:6);
  Phi = reshape(Y(end, 7:end), 6, 6);
  f_x = fun_cr3bp([], X, mu);
  F = [X(end, 2); X(end, 4); X(end, 6)];
  F_PAC = ([x02(1); x02(3); x02(5); t02] - [x01_ast(1); x01_ast(3); x01_ast(5); t01_ast])'*delta - scale;
  G = [F; F_PAC];
  DF = [Phi(2, 1), Phi(2, 3), Phi(2, 5), f_x(2);
        Phi(4, 1), Phi(4, 3), Phi(4, 5), f_x(4);
        Phi(6, 1), Phi(6, 3), Phi(6, 5), f_x(6)];
  DG = [DF; delta'];
  X_ast = [x02(1); x02(3); x02(5); t02] - DG\G;

  x02_ast = [X_ast(1); 0; X_ast(2); 0; X_ast(3); 0];
  t02_ast = X_ast(4);
  C = Jacobi_const(x02_ast, mu);
=======
% differential correction for halo family in CR3BP by PAC
function [x02_ast, t02_ast, C, G] = fun_differential_correction_cr3bp_PAC(x01_ast, t01_ast, x02, t02, scale, delta, mu, options)
% x01_ast : previously-converged x0 ( (i)th orbit )
% t01_ast : previously-converged period ( (i)th orbit )
% x02     : initial guess of x0 ( (i+1)th orbit )
% t02     : initial guess of period ( (i+1)th orbit )
% scale   : scaling factor of delta
% delta   : null space of the previously-converged Jacobian matrix
% mu      : mass ratio of the primaries
% options : options for ode

  X0 = [x02(1), x02(2), x02(3), x02(4), x02(5), x02(6), reshape(eye(6), 1, [])];
  tspan = [0 t02];
  [~, Y] = ode113(@(t, x) fun_stm_cr3bp(t, x, mu), tspan, X0, options);

  X = Y(end, 1:6);
  Phi = reshape(Y(end, 7:end), 6, 6);
  f_x = fun_cr3bp([], X, mu);
  F = [X(end, 2); X(end, 4); X(end, 6)];
  F_PAC = ([x02(1); x02(3); x02(5); t02] - [x01_ast(1); x01_ast(3); x01_ast(5); t01_ast])'*delta - scale;
  G = [F; F_PAC];
  DF = [Phi(2, 1), Phi(2, 3), Phi(2, 5), f_x(2);
        Phi(4, 1), Phi(4, 3), Phi(4, 5), f_x(4);
        Phi(6, 1), Phi(6, 3), Phi(6, 5), f_x(6)];
  DG = [DF; delta'];
  X_ast = [x02(1); x02(3); x02(5); t02] - DG\G;

  x02_ast = [X_ast(1); 0; X_ast(2); 0; X_ast(3); 0];
  t02_ast = X_ast(4);
  C = Jacobi_const(x02_ast, mu);
>>>>>>> 748834e2ef04e45bd3c365491093c08a43a6640e
end