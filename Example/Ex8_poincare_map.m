% poincare map at x = 1 - mu

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

myTimer = tic;        %start timer


%% initial settings
[mu,~,~,~] = parameter(6); % Sun-Jupiter
[L1,L2,L3,L4,L5] = librationPoints(mu);

options_ODE = odeset('RelTol',3e-14, 'AbsTol',1e-14);
options_ODE_poincare = odeset('RelTol',3e-14, 'AbsTol',1e-14,'Events',@(t,x) fun_odestop_CR3BP_Body2_vertical(t,x,mu) );



%% differential correction
iteration_DC_max = 100;
threshold = 1e-12;

% L1 Lyapunov orbit
x0_1 = [0.926672517344211 0 0  0 0.045684284381193 0]';
t0_1 = 2.915640890603734/2;
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
x_L1 = x_corrected; % C_L1 = 3.036999999999999
t_L1 = t_n_1;

% L2 Lyapunov orbit
x0_2 = [1.072915577273646 0 0 0 -0.025945940908770 0]';
t0_2 = 3.184772418214781/2;
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
x_L2 = x_corrected; % C_L2 = 3.036969030648278
t_L2 = t_n_2;



%% zero velocity curve
[x_curve,y_curve] = meshgrid(0.9:1e-4:1.1,-0.1:1e-4:0.1);
[a,b] = size(x_curve);
z_curve = zeros(a,b);

x_curve = reshape(x_curve, [a*b,1]);
y_curve = reshape(y_curve, [a*b,1]);
z_curve = reshape(z_curve, [a*b,1]);

r1 = sqrt((x_curve+mu).^2+y_curve.^2+z_curve.^2);
r2 = sqrt((x_curve-(1-mu)).^2+y_curve.^2+z_curve.^2);
U = 1/2.*(x_curve.^2+y_curve.^2) + (1-mu)./r1 + mu./r2;
C = 2.*U;

x_curve = reshape(x_curve, [a,b]);
y_curve = reshape(y_curve, [a,b]);
C = reshape(C, [a,b]);



%% manifold
tspan_u = [0 9];
tspan_s = [10 0];
xpert = 1e-10;
N = 40;

% L1
[~,~,~, XU_right,~] = fun_manifold_cr3bp(mu, x_L1(1,:)', 2*t_L1, N, xpert, options_ODE);

% L2
[XS_left,~,~,~,~] = fun_manifold_cr3bp(mu, x_L2(1,:)', 2*t_L2, N, xpert, options_ODE);



%% poincare map
yu_right_f1 = zeros(N,6);
yu_right_f2 = NaN(N,6);
ys_left_f1 = zeros(N,6);
ys_left_f2 = NaN(N,6);

f1 = figure();
hold on
contourf(x_curve, y_curve, -C, [-C_L1 -C_L1],'k');
colormap gray
for i = 1:N
  [~,yu_right,te,xe,ie] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan_u, XU_right(:,i), options_ODE_poincare);
  f1_p1 = plot(yu_right(:,1),yu_right(:,2),'Color','r');
  yu_right_f1(i,:) = xe(1,:);

  if size(xe,1) < 2
    continue;
  end
  yu_right_f2(i,:) = xe(2,:);
end
for i = 1:N
  [~,ys_left,te,xe,ie] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan_s, XS_left(:,i), options_ODE_poincare);
  ys_left = flipud(ys_left);
  f1_p2 = plot(ys_left(:,1),ys_left(:,2),'Color','#77AC30');
  ys_left_f1(i,:) = xe(1,:);

  if size(xe,1) < 2
    continue;
  end
  ys_left_f2(i,:) = xe(2,:);
end
plot(x_L1(:,1),x_L1(:,2),'k');
plot(x_L2(:,1),x_L2(:,2),'k');
plot(L1(1),L1(2),'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',15);
plot(L2(1),L2(2),'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',15);
xline(1-mu,'--k',"$x = 1 - \mu$",'Linewidth',1.5,'Fontsize',24,'Interpreter','latex');
view(gca,0,90);
axis equal
xlabel('$x$[-]');
ylabel('$y$[-]');
xlim([0.9 1.1]);
ylim([-0.08 0.08]);
grid on
legend([f1_p1, f1_p2], {'Unstable manifold','Stable manifold'});
hold off
f1_name = strcat('Ex8_poincare_map_traj_mu=',num2str(mu),'_xpert=',num2str(xpert));
f1_name = strrep(f1_name,'.',',');
save_fig(f1,f1_name,[0 90]);


f2 = figure();
hold on
f2_p1 = plot(yu_right_f1(:,2),yu_right_f1(:,5),'Color','r');
f2_p2 = plot(ys_left_f1(:,2),ys_left_f1(:,5),'Color','#77AC30');
f2_p3 = plot(yu_right_f2(:,2),yu_right_f2(:,5),'Color','m');
f2_p4 = plot(ys_left_f2(:,2),ys_left_f2(:,5),'Color','g');
axis square
xlabel('$y$[-]','Interpreter','latex');
ylabel('$\dot{y}$[-]','Interpreter','latex');
xlim([-0.03 0.05]);
ylim([-0.3 0.6]);
grid on
legend([f2_p1, f2_p2, f2_p3, f2_p4], {'1st insertion of unstable manifold','1st insertion of stable manifold',...
    '2nd insertion of unstable manifold','2nd insertion of stable manifold'},'Location','northeastoutside');
hold off
f2_name = strcat('Ex8_poincare_map_mu=',num2str(mu),'_xpert=',num2str(xpert));
f2_name = strrep(f2_name,'.',',');
save_fig(f2,f2_name,[0 90]);



%% End of script
time = strcat('calculation time: ', num2str(toc(myTimer)));
disp(time);
