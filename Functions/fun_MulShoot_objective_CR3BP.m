<<<<<<< HEAD
function J = fun_MulShoot_objective_CR3BP(X, r0, rf, n, mu,options)
% X       : design variable vector
% r0      : position of the initial point
% rf      : position of the final point
% n       : the number of the patch points
% mu      : mass ratio of the primaries
% options : options for ode

  if(size(r0,2) > 1), r0 = r0'; end
  if(size(rf,2) > 1), rf = rf'; end

  X0 = [r0; X(1:(6*n-9)); rf; X((6*n-8):end)];
  Xf = zeros(6*(n-1),1);
  for i = 1:(n-1)
    tspan = [0 X0(6*n+i)];
    [~,Y] = ode113(@(t,x) fun_cr3bp(t,x,mu),tspan,X0((6*i-5):6*i),options);

    Xf((6*i-5):6*i) = Y(end,1:6)';
  end

  F = zeros(6*(n-1),1);
  for i = 1:(n-1)
    F((6*i-5):6*i) = Xf((6*i-5):6*i) - X0((6*i+1):6*(i+1));
  end

  J = norm(F);
=======
function J = fun_MulShoot_objective_CR3BP(X, r0, rf, n, mu,options)
% X       : design variable vector
% r0      : position of the initial point
% rf      : position of the final point
% n       : the number of the patch points
% mu      : mass ratio of the primaries
% options : options for ode

  if(size(r0,2) > 1), r0 = r0'; end
  if(size(rf,2) > 1), rf = rf'; end

  X0 = [r0; X(1:(6*n-9)); rf; X((6*n-8):end)];
  Xf = zeros(6*(n-1),1);
  for i = 1:(n-1)
    tspan = [0 X0(6*n+i)];
    [~,Y] = ode113(@(t,x) fun_cr3bp(t,x,mu),tspan,X0((6*i-5):6*i),options);

    Xf((6*i-5):6*i) = Y(end,1:6)';
  end

  F = zeros(6*(n-1),1);
  for i = 1:(n-1)
    F((6*i-5):6*i) = Xf((6*i-5):6*i) - X0((6*i+1):6*(i+1));
  end

  J = norm(F);
>>>>>>> 748834e2ef04e45bd3c365491093c08a43a6640e
end