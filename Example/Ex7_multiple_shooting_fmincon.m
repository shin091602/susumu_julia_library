% multiple shooting

%% Start of script
close all;  %close all figures
clear;      %clear all variables
clc;        %clear the command terminal
format long
%warning off

% line width
set(0,'DefaultLineLineWidth',1.5) % default 0.5pt
set(0,'DefaultAxesLineWidth',1.5)
set(0,'DefaultTextLineWidth',1.5)

% font size
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)

% font name
set(0,'DefaultTextFontName','Times New Roman')
set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultTextInterpreter','Latex')
set(0,'DefaultLegendInterpreter','Latex')

% figure color
set(0,'DefaultFigureWindowStyle','docked');
set(gcf,'Color','none');
set(gca,'Color','none');
set(gcf,'InvertHardCopy', 'off');

close

current_pass = pwd;
addpath(replace(current_pass, 'Examples', 'Functions'));
addpath(replace(current_pass, 'Examples', 'planet3D'));

myTimer = tic;        %start timer


%% initial settings
[mu,~,~,~] = parameter(1); % Sun-Earth
[L1,L2,L3,L4,L5] = librationPoints(mu);

options_ODE = odeset('RelTol',3e-14, 'AbsTol',1e-14);

options_planet3D.Position = [1-mu, 0, 0]';
options_planet3D.RefPlane = 'ecliptic';
options_planet3D.RotAngle = 135;
options_planet3D.Units = 6378.1363e3*1000/1.5;



%% differential correction
iteration_DC_max = 100;
threshold = 1e-12;

% L1 Lyapunov orbit
x0_1 = [0.987933580858119 0 0 0 0.018052674638026 0]';
t0_1 = 3.290337439058733/2;
for iteration = 1:iteration_DC_max
  [x_n_1,t_n_1,C_L1] = fun_differential_correction_cr3bp(x0_1,t0_1,mu,options_ODE);

  tspan = [0 2*t_n_1];
  [t_corrected, x_corrected] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan, x_n_1, options_ODE);

  x_error = norm(x_corrected(end,:) - x_corrected(1,:));
  %disp( strcat('x_error = ',num2str(x_error)) );
  if x_error < threshold
    break;
  end

  if x_error > 1e+3
    disp('calculation diverged');
    return;
  end

  if iteration == iteration_DC_max
    disp('do not finish');
    return;
  end

  x0_1 = x_n_1;
  t0_1 = t_n_1;
end
x_L1 = x_corrected; % C_L1 = 3.000599995882406
t_L1 = t_n_1;

% L2 Lyapunov orbit
x0_2 = [1.012164676114102 0 0 0 -0.017947125645568 0]';
t0_2 = 3.319976286115385/2;
for iteration = 1:iteration_DC_max
  [x_n_2,t_n_2,C_L2] = fun_differential_correction_cr3bp(x0_2,t0_2,mu,options_ODE);

  tspan = [0 2*t_n_2];
  [t_corrected, x_corrected] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan, x_n_2, options_ODE);

  x_error = norm(x_corrected(end,:) - x_corrected(1,:));
  %disp( strcat('x_error = ',num2str(x_error)) );
  if x_error < threshold
    break;
  end

  if x_error > 1e+3
    disp('calculation diverged');
    return;
  end

  if iteration == iteration_DC_max
    disp('do not finish');
    return;
  end

  x0_2 = x_n_2;
  t0_2 = t_n_2;
end
x_L2 = x_corrected; % C_L2 = 3.000600000000000
t_L2 = t_n_2;



%% initial trajectory
x1 = [0.985000000000000 0.001092112291766 0 0.022 0.005 0]';
x2 = [1-mu -4e-3 0 0.03 0 0]';
x3 = [1.015000000000000 0.001092400176687 0 0.020764610600732 -0.004967672676993 0]';
delta_t1 = 0.75;
delta_t2 = 0.75;

tspan1 = [0 delta_t1];
[~,x_traj_1] = ode113(@(t,x) fun_cr3bp(t,x,mu),tspan1,x1,options_ODE);

tspan2 = [0 delta_t2];
[~,x_traj_2] = ode113(@(t,x) fun_cr3bp(t,x,mu),tspan2,x2,options_ODE);

f1 = figure();
hold on
planet3D('Earth Cloudy',options_planet3D);
plot(x_L1(:,1),x_L1(:,2),'--k');
plot(x_L2(:,1),x_L2(:,2),'--k');
plot(x_traj_1(:,1),x_traj_1(:,2),'Color','#0072BD');
plot(x_traj_2(:,1),x_traj_2(:,2),'Color','#0072BD');
plot(x1(1),x1(2),'o','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',10);
plot(x2(1),x2(2),'o','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',10);
plot(x3(1),x3(2),'o','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',10);
plot(L1(1),L1(2),'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',15);
plot(L2(1),L2(2),'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',15);
view(gca,0,90);
axis equal
xlabel('$x$[-]');
ylabel('$y$[-]');
xlim([0.985 1.015]);
grid on
hold off
f1_name = strcat('Ex7_initial_guess_mu=',num2str(mu));
f1_name = strrep(f1_name,'.',',');
save_fig(f1,f1_name,[0 90]);

f2 = figure();
hold on
planet3D('Earth Cloudy',options_planet3D);
plot(x_L1(:,1),x_L1(:,2),'--k');
plot(x_L2(:,1),x_L2(:,2),'--k');
f2_p1 = plot(x_traj_1(:,1),x_traj_1(:,2),'Color','#0072BD');
plot(x_traj_2(:,1),x_traj_2(:,2),'Color','#0072BD');
plot(x1(1),x1(2),'o','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',10);
plot(x2(1),x2(2),'o','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',10);
plot(x3(1),x3(2),'o','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',10);
plot(L1(1),L1(2),'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',15);
plot(L2(1),L2(2),'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',15);
hold off



%% multiple shooting
X0 = [x1(4:6); x2; x3(4:6); delta_t1; delta_t2];
r0 = x1(1:3);
rf = x3(1:3);
ToF = delta_t1 + delta_t2;
n = 3;

% Optimization with fmincon
problem.x0          = X0;
problem.A           = [];
problem.b           = [];
problem.Aeq         = [];
problem.beq         = [];
problem.ub          = [];
problem.lb          = [];
problem.nonlcon     = @(x) fun_MulShoot_const_CR3BP(x, r0, rf, n, ToF, mu, options_ODE);
problem.objective   = @(x) fun_MulShoot_objective_CR3BP(x, r0, rf, n, mu, options_ODE);
problem.solver      = 'fmincon';
problem.options     = optimoptions(@fmincon,'Algorithm','sqp', 'ConstraintTolerance',1e-14, 'Display','Iter', 'MaxFunctionEvaluations',1e+4, 'MaxIterations',1000, 'ObjectiveLimit',1e-12, 'StepTolerance',1e-14);
[x_ms, fval, exitflag, output] = fmincon(problem);
if exitflag <= 0
  disp('Error in fmincon');
  return;
end
X_n = [r0; x_ms(1:(6*n-9)); rf; x_ms((6*n-8):end)];
save('Ex7_solution.mat','X_n','fval');

figure(2);
hold on
for i = 1:(n-1)
  tspan = [0 X_n(6*n+i)];
  [~,X] = ode113(@(t,x) fun_cr3bp(t,x,mu),tspan,X_n((6*i-5):6*i),options_ODE);
  f2_p2 = plot(X(:,1),X(:,2),'Color','#D95319');
  plot(X(1,1),X(1,2),'o','MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize',10);
  plot(X(end,1),X(end,2),'o','MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize',10);
end
view(gca,0,90);
xlabel('$x$[-]');
ylabel('$y$[-]');
grid on
legend([f2_p1, f2_p2], {'Initial guess','Multiple shooting'}, 'Location','northeast');
hold off
f2_name = strcat('Ex7_multiple_shooting_mu=',num2str(mu));
f2_name = strrep(f2_name,'.',',');
save_fig(f2,f2_name,[0 90]);



%% End of script
time = strcat('calculation time: ', num2str(toc(myTimer)));
disp(time);
