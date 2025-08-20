% multiple shooting
function X_n = fun_multiple_shooting_cr3bp(X0,n,mu,options)
% X0      : initial guess for the design variable vector X
% n       : the number of the patch points
% mu      : mass ratio of the primaries
% options : options for ode

  Xf = zeros(6*(n-1),1);
  Phi = zeros(6,6,n-1);
  for i = 1:(n-1)
    Y0 = [X0((6*i-5):6*i); reshape(eye(6),[],1)];
    tspan = [0 X0(6*n+i)];
    [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,mu),tspan,Y0,options);

    Xf((6*i-5):6*i) = Y(end,1:6)';
    Phi(:,:,i) = reshape(Y(end,7:end),6,6);
  end

  F = zeros(6*(n-1),1);
  for i = 1:(n-1)
    F((6*i-5):6*i) = Xf((6*i-5):6*i) - X0((6*i+1):6*(i+1));
  end

  DF = zeros(6*(n-1),7*n-1);
  I = eye(6);
  for i = 1:(n-1)
    DF((6*i-5):6*i, (6*i-5):6*i) = Phi(:,:,i);
    DF((6*i-5):6*i, (6*i+1):6*(i+1)) = -I;

    DF((6*i-5):6*i, 6*n+i) = fun_cr3bp(X0(6*n+i),Xf((6*i-5):6*i),mu);
  end
  X_ast = [X0(4:6*(n-1)); X0((6*n-2):end)] - [DF(:,4:6*(n-1)) DF(:,(6*n-2):end)]\F;

  X_n = [X0(1:3); X_ast(1:(6*n-9)); X0((6*n-5):(6*n-3)); X_ast((6*n-8):end)];
end