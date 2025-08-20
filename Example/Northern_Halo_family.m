%% Pseudo arclength continuation for northern halo family around L2 in Sun-Earth CR3BP
%% By:Soi Yamaguchi
% Arclength continuation

%% Start of script
close all;  %close all figures
clear;      %clear all variables
clc;        %clear the command terminal
format long
%warning off

% line width
set(0, 'DefaultLineLineWidth', 0.8) % default 0.5pt
set(0, 'DefaultAxesLineWidth', 0.8)
set(0, 'DefaultTextLineWidth', 0.8)

% font size
set(0, 'DefaultTextFontSize', 13)
set(0, 'DefaultAxesFontSize', 13)

% font name
set(0, 'DefaultTextFontName', 'Times New Roman')
set(0, 'DefaultAxesFontName', 'Times New Roman')
set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')

% figure color
set(0, 'DefaultFigureWindowStyle', 'docked');
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gcf, 'InvertHardCopy', 'off');

close

current_pass = pwd;
addpath(strcat(current_pass, '/Functions'));

myTimer = tic;        %start timer

p = dictionary();

%% DICTIONARY OF SUN-EARTH CRTBP VARIABLES
p("Ms") = 1.9885e+30;
p("M1") = 5.9724e+24;
p("mu") = 3.00346093142064e-06;%SE
%p("mu") = 0.012153599037880;%NRHO--EM
p("G") = 6.67300000000000e-11;
p("iteration_limit") = 50;
p("threshold") = 1e-10;
p("xpert") = 1e-7;
p("tf_L2") = 10;
p("scale") = 1e-4;
count = 0;
p("count_max")     = 8000;
p("plot_interval") = 200;
p("Jacobi_min") = 3.00018403545658;
p("Jacobi_max") = 3.0008942050035974;
p("chara_length_CR3BP_SE") = 149600000000.000;
dim_L = p("chara_length_CR3BP_SE")/1000;
p("chara_mass_CR3BP_SE") = 1.98850597240000e+30;
p("chara_time_CR3BP_SE") = 5023117.97170552;
p("N_CR3BP_SE") = 1.99079536979392e-07;

%% INITIALIZATION
options_ODE   = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);
color      = jet;
Jacobi_lim = linspace(p("Jacobi_min"), p("Jacobi_max"), size(color,1));
Interp_c   = griddedInterpolant(Jacobi_lim, color);
Jacobi1       = zeros(p("count_max")/p("plot_interval")+1, 1);
x0_corrected1 = zeros(p("count_max")/p("plot_interval")+1, 6);
t0_corrected1 = zeros(p("count_max")/p("plot_interval")+1, 1);

%% LIBRATION POINT L2
[~,L2,~,~,~] = librationPoints(p("mu"));

%% INITIAL PERIODIC ORBIT
x0 = [1.01098476029841;0;0.00425821353748602;0;-0.0110580081896328;0];
T = 3.08315887305965;

%orbital period
t_n_L2 = T/2;
x_n_L2 = x0;

%% PSEUDO-ARCLENGTH CONTINUATION
hh = figure();
hold on
grid on
box on
view([-39 31]);
tspan1  = [0, 2*t_n_L2];
[~, x1] = ode113(@(t, x) fun_cr3bp(t, x, p("mu")), tspan1, x_n_L2, options_ODE);
Jacobi1(1)          = Jacobi_const(x_n_L2, p("mu"));
x0_corrected1(1, :) = x_n_L2';
t0_corrected1(1)    = t_n_L2;
rgb = Interp_c(Jacobi1(1));
plot3(x1(:, 1), x1(:, 2), x1(:, 3), 'Color', rgb);

[delta, ~] = fun_null_cr3bp(x_n_L2, t_n_L2, p("mu"), options_ODE);%extract null space of initial values
delta      = -delta;
x02        = x_n_L2;
x02(1)     = x_n_L2(1) + p("scale")*delta(1);
x02(3)     = x_n_L2(3) + p("scale")*delta(2);
x02(5)     = x_n_L2(5) + p("scale")*delta(3);
t02        = t_n_L2 + p("scale")*delta(4);

while 1
count = count + 1;
disp(strcat('count = ', num2str(count)));

% differential correction
for iteration = 1:p("iteration_limit")
  [x02_ast, t02_ast, C, G] = fun_differential_correction_cr3bp_PAC(x_n_L2, t_n_L2, x02, t02, p("scale"), delta, p("mu"), options_ODE);

  tspan2 = [0, 2*t02_ast];
  [t_corrected, x_corrected] = ode113(@(t, x) fun_cr3bp(t, x, p("mu")), tspan2, x02_ast, options_ODE);

  if norm(G) < p("threshold")
    break;
  end

  if norm(G) > 1e+3
    disp('calculation diverged');
    return;
  end

  if iteration == p("iteration_limit")
    disp('do not finish');
    return;
  end

  x02 = x02_ast;
  t02 = t02_ast;
end

  % plot orbit
  if (count >= p("plot_interval")) && (mod(count, p("plot_interval")) == 0)
      Jacobi1(count/p("plot_interval")+1)          = C;
      x0_corrected1(count/p("plot_interval")+1, :) = x_corrected(1, :); %RHS:initial cond
      t0_corrected1(count/p("plot_interval")+1)    = t_corrected(end);
      rgb = Interp_c(C);
      plot3(x_corrected(:, 1), x_corrected(:, 2), x_corrected(:, 3), 'Color', rgb,'LineWidth',1.5);
  end

  x_n_L2 = x02_ast;
  t_n_L2 = t02_ast;

  if count == p("count_max")
      break
  end
end

%% PAC to the primaries' direction
count = 0;

% orbital period
t_n_L2 = T/2;
x_n_L2 = x0;

% initialization
Jacobi2       = zeros(p("count_max")/p("plot_interval")+1, 1);
x0_corrected2 = zeros(p("count_max")/p("plot_interval")+1, 6);
t0_corrected2 = zeros(p("count_max")/p("plot_interval")+1, 1);

% initial state
tspan1  = [0, 2*t_n_L2];
[~, x1] = ode113(@(t, x) fun_cr3bp(t, x, p("mu")), tspan1, x_n_L2, options_ODE);
Jacobi2(1)          = Jacobi_const(x_n_L2, p("mu"));
x0_corrected2(1, :) = x_n_L2';
t0_corrected2(1)    = t_n_L2;
rgb = Interp_c(Jacobi2(1));

% null space
[delta, ~] = fun_null_cr3bp(x_n_L2, t_n_L2, p("mu"), options_ODE);%extract null space of initial values
x02        = x_n_L2;
x02(1)     = x_n_L2(1) + p("scale")*delta(1);
x02(3)     = x_n_L2(3) + p("scale")*delta(2);
x02(5)     = x_n_L2(5) + p("scale")*delta(3);
t02        = t_n_L2 + p("scale")*delta(4);

while 1
  count = count + 1;
  disp(strcat('count = ', num2str(count)));

  % differential correction
  for iteration = 1:p("iteration_limit")
    [x02_ast, t02_ast, C, G] = fun_differential_correction_cr3bp_PAC(x_n_L2, t_n_L2, x02, t02, p("scale"), delta, p("mu"), options_ODE);

    tspan2 = [0, 2*t02_ast];
    [t_corrected, x_corrected] = ode113(@(t, x) fun_cr3bp(t, x, p("mu")), tspan2, x02_ast, options_ODE);

    if norm(G) < p("threshold")
      break;
    end

    if norm(G) > 1e+3
      disp('calculation diverged');
      return;
    end

    if iteration == p("iteration_limit")
      disp('do not finish');
      return;
    end

    x02 = x02_ast;
    t02 = t02_ast;
  end

  % plot orbit
  if (count >= p("plot_interval")) && (mod(count, p("plot_interval")) == 0)
    Jacobi2(count/p("plot_interval")+1)          = C;
    x0_corrected2(count/p("plot_interval")+1, :) = x_corrected(1, :); %RHS:initial cond
    t0_corrected2(count/p("plot_interval")+1)    = t_corrected(end);
    rgb = Interp_c(C);
    plot3(x_corrected(:, 1), x_corrected(:, 2), x_corrected(:, 3), 'Color', rgb,'LineWidth',1.5);
    % Lyapunov orbit
    if x_corrected(1, 3)<p("threshold")
      idxlim = count/p("plot_interval")+1;
      break;
    end
  end

  x_n_L2 = x02_ast;
  t_n_L2 = t02_ast;

  if count == p("count_max")
    break
  end

end

h_p1 = plot3(L2(1), L2(2), L2(3), '*', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
axis equal
colormap jet;
c = colorbar;
ylabel(c, 'Jacobi constant [-]', 'FontSize', 15);
caxis([p("Jacobi_min") p("Jacobi_max")]);
c.Ticks = linspace(p("Jacobi_min"), p("Jacobi_max"), 6);
xlabel('$x$[-]');
ylabel('$y$[-]');
zlabel('$z$[-]');
%legend([h_p1, h_p2], {'$L_2$', 'Earth'}, 'Location', 'northwest');
hold off

%% Arrange orbit data
orbit_data1 = zeros(p("count_max")/p("plot_interval")+1,6+2);
orbit_data2 = zeros(idxlim-1,6+2);
for i=1:p("count_max")/p("plot_interval")+1
  orbit_data1(i,:) = [x0_corrected1(i,:),t0_corrected1(i),Jacobi1(i)];
end
for i=2:idxlim
  orbit_data2(i-1,:) = [x0_corrected2(i,:),t0_corrected2(i),Jacobi2(i)];
end
orbit_data = [flipud(orbit_data1);orbit_data2];

%% Stability color plot
hSt = figure();
hold on
grid on
box on
xlabel('$x$[-]');
ylabel('$y$[-]');
zlabel('$z$[-]');
view([-39 31]);
axis equal
% libration point
plot3(L2(1), L2(2), L2(3),'*', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
% Monodromy Matrix
M = zeros(6,6,p("count_max")/p("plot_interval"));
% Center manifold --altenatives of quasi-periodic tori
orbit_cent = cell(p("count_max")/p("plot_interval")+1,1);
% fifth component -- eigenvalues
eigen5 = zeros(length(orbit_data),1);
for i=1:length(orbit_data)
    Y0 = [orbit_data(i,1:6),reshape(eye(6),1,[])];
    [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,p("mu")),[0 2*orbit_data(i,7)],Y0',options_ODE);
    M(:,:,i) = reshape(Y(end,7:42),6,6);
    % eigenvectors and eignvalues
    [V,D] = eig(M(:,:,i));
    % stability index
    for j=1:6
      if (abs(abs(D(j,j))-1)<1e-3) && (imag(D(j,j))~=0)
        % corresponds to center direciton
        fcen = plot3(Y(:,1),Y(:,2),Y(:,3),'b','Linewidth',1.5);
        % save initial state
        orbit_cent{i,1} = Y(1,1:6);
        break;
      else
        funcen = plot3(Y(:,1),Y(:,2),Y(:,3),'r','Linewidth',1.5);
      end
    end
end
hold off

%% Bifurcation
hB = figure();
hold on
grid on
axis square
xlim([-10 10]);
ylim([-10 10]);

% Monodromy Matrix
M = zeros(6,6,p("count_max")/p("plot_interval"));
alpha = zeros(p("count_max")/p("plot_interval")+1,1);
beta = zeros(p("count_max")/p("plot_interval")+1,1);
for i=2:p("count_max")/p("plot_interval")+1
  Y0 = [orbit_data(i,1:6),reshape(eye(6),1,[])];
  [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,p("mu")),[0 orbit_data(i,6+1)],Y0',options_ODE);
  M(:,:,i) = reshape(Y(end,7:42),6,6);
  alpha(i) = 2-trace(M(:,:,i));
  beta(i) = 0.5.*(alpha(i)^2+2-trace(M(:,:,i))^2);
end
% patch1
px = [0 -4 -10 -10 -6];
py = [-2 -10 -10 10 10];
patch(px,py,[0.69, 0.88, 0.62]);
hold  on
% patch2
px = [0 -6 6];
py = [-2 10 10];
patch(px,py,[0.95, 0.76, 0.55]);
hold  on
% patch3
px = [0 -4 4];
py = [-2 -10 -10];
patch(px,py,[0.62, 0.6, 0.62]);
hold  on
% patch3
px = [0 6 10 10 4];
py = [-2 10 10 -10 -10];
patch(px,py,[0.91, 0.68, 0.82]);
hold  on

% diagram
plot(alpha(2:end),beta(2:end),'b','LineWidth',2);

% ref
x = -20:1e-1:10;
y1 = -2.*x-2;
y2 = 2.*x-2;
y3 = x.^2/4+2;
y4 = 10;
y5 = x+1;

plot(x,y1,'-k','LineWidth',1.5);
plot(x,y2,'-k','LineWidth',1.5);
plot(x,y3,'-k','LineWidth',1.5);
plot(x,y5,'--k','LineWidth',1.5);
hold off

%% SAVE RESULTS
%cd fig/family/
%h_name = strcat('halofamily_SE_mu=', num2str(p("mu")), '_x0=', num2str(x0_corrected1(1,1)), '_z0=', num2str(x0_corrected1(1,3)),...
%'_ydot0=', num2str(x0_corrected1(1,5)), '_t0=', num2str(t0_corrected1(1)), '_scale=', num2str(p("scale")),'_count_max=', num2str(p("count_max")), '_plot_interval=', num2str(p("plot_interval")),'C_min=',num2str(p("Jacobi_min")),'C_max=',num2str(p("Jacobi_max")));
%h_name = strrep(h_name,'.',',');
%save_fig(hh,h_name,[-39 31]);

%hSt_name = strcat('halofamily_stability_mu=', num2str(p("mu")), '_x0=', num2str(x0_corrected1(1,1)), '_z0=', num2str(x0_corrected1(1,3)),...
%'_ydot0=', num2str(x0_corrected1(1,5)), '_t0=', num2str(t0_corrected1(1)), '_scale=', num2str(p("scale")),'_count_max=', num2str(p("count_max")), '_plot_interval=', num2str(p("plot_interval")),'C_min=',num2str(p("Jacobi_min")),'C_max=',num2str(p("Jacobi_max")));
%hSt_name = strrep(hSt_name,'.',',');
%save_fig(hSt,hSt_name,[-39 31]);

%h_name = strcat('halofamily_hodograph_SE_mu=', num2str(p("mu")), '_x0=', num2str(x0_corrected1(1,1)), '_z0=', num2str(x0_corrected1(1,3)),'_ydot0=', num2str(x0_corrected1(1,5)), '_t0=', num2str(t0_corrected1(1)), '_count_max=', num2str(p("count_max")), '_plot_interval=', num2str(p("plot_interval")),'C_min=',num2str(p("Jacobi_min")),'C_max=',num2str(p("Jacobi_max")));
%h_name = strrep(h_name,'.',',');
%save_fig(hFh,h_name,[0 90]);

%broucke_stability_diagram
%hBf = strcat('broucke_stability_diagram_halofamily_SE_mu=', num2str(p("mu")), '_x0=', num2str(x0_corrected1(1,1)), '_z0=', num2str(x0_corrected1(1,3)),'_ydot0=', num2str(x0_corrected1(1,5)), '_t0=', num2str(t0_corrected1(1)), '_count_max=', num2str(p("count_max")), '_plot_interval=', num2str(p("plot_interval")),'C_min=',num2str(p("Jacobi_min")),'C_max=',num2str(p("Jacobi_max")));
%hBf = strrep(hBf,'.',',');
%save_fig(hB,hBf,[0 90]);
%cd ..
%cd ..

%% End of script
time = strcat('calculation time: ', num2str(toc(myTimer)));
disp(time);
