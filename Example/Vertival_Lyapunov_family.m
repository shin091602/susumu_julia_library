%% Pseudo arclength continuation for northern vertical Lyapunov family around L2 in Sun-Earth CR3BP
%% By:Soi Yamaguchi

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

%% INITIAL PERIODIC ORBIT
x0 = [1.0029239047032787;0;0;0;-4.8187208188144461e-2;-2.9416770480063559e-1];
T = 6.2808090883636165;

p = dictionary();
%% DICTIONARY OF SUN-EARTH CRTBP VARIABLES
p("Ms") = 1.9885e+30;
p("M1") = 5.9724e+24;
p("mu") = 3.00346093142064e-06;%SE
p("G") = 6.67300000000000e-11;
p("iteration_limit") = 100;
p("threshold") = 1e-10;
p("xpert") = 1e-7;
p("tf_L2") = 10;
p("scale") = 5e-5;
count = 0;
p("count_max")     = 3000;
p("plot_interval") = 75;
p("Jacobi_min") = 2.913209300553006;
p("Jacobi_max") = 2.979596689425427;

p("chara_length_CR3BP_SE") = 149600000000.000;
dim_L = p("chara_length_CR3BP_SE")/1000;
p("chara_mass_CR3BP_SE") = 1.98850597240000e+30;
p("chara_time_CR3BP_SE") = 5023117.97170552;
p("N_CR3BP_SE") = 1.99079536979392e-07;
xpert_dim = p("xpert")*p("chara_length_CR3BP_SE")/p("chara_time_CR3BP_SE");

%% INITIALIZATION
options_ODE   = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);
color      = jet;
Jacobi_lim = linspace(p("Jacobi_min"), p("Jacobi_max"), size(color,1));
Interp_c   = griddedInterpolant(Jacobi_lim, color);
Jacobi       = zeros(p("count_max")/p("plot_interval")+1, 1);
x0_corrected = zeros(p("count_max")/p("plot_interval")+1, 6);
t0_corrected = zeros(p("count_max")/p("plot_interval")+1, 1);

%% LIBRATION POINT L2
[~,L2,~,~,~] = librationPoints(p("mu"));

%orbital period
t_n_L2 = T/2;
x_n_L2 = x0;

%% PSEUDO-ARCLENGTH CONTINUATION
hvL = figure();
hold on
grid on
box on
xlabel('$x$[-]');
ylabel('$y$[-]');
zlabel('$z$[-]');
axis square
view([-39 32]);
xlim([0.949274012415507 1.013359974485568]);
ylim([-0.039852272743994 0.035931818165097]);

% color setting
colormap jet;
c = colorbar;
ylabel(c, 'Jacobi constant [-]', 'FontSize', 15);
caxis([p("Jacobi_min") p("Jacobi_max")]);
c.Ticks = linspace(p("Jacobi_min"), p("Jacobi_max"), 6);

tspan1  = [0, 2*t_n_L2];
[~, x1] = ode113(@(t, x) fun_cr3bp(t, x, p("mu")), tspan1, x_n_L2, options_ODE);
Jacobi(1)          = Jacobi_const(x_n_L2, p("mu"));
x0_corrected(1, :) = x_n_L2';
t0_corrected(1)    = t_n_L2;
rgb = Interp_c(Jacobi(1));
%plot3(x1(:, 1), x1(:, 2), x1(:, 3), 'Color', 'g');

[delta, ~] = fun_null_vertical_Lyapunov(x_n_L2, t_n_L2, p("mu"), options_ODE);%extract null space of initial values
delta      = delta;
x02        = x_n_L2;
x02(1)     = x_n_L2(1) + p("scale")*delta(1);
x02(5)     = x_n_L2(5) + p("scale")*delta(2);
x02(6)     = x_n_L2(6) + p("scale")*delta(3);
t02        = t_n_L2 + p("scale")*delta(4);

while 1
count = count + 1;
disp(strcat('count = ', num2str(count)));

% differential correction
for iteration = 1:p("iteration_limit")
  [x02_ast, t02_ast, C, G] = fun_differential_correction_vertical_Lyapunov_PAC(x_n_L2, t_n_L2, x02, t02, p("scale"), delta, p("mu"), options_ODE);

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
      Jacobi(count/p("plot_interval")+1)          = C;
      x0_corrected(count/p("plot_interval")+1, :) = x_corrected(1, :); %RHS:initial cond
      t0_corrected(count/p("plot_interval")+1)    = t_corrected(end);
      rgb = Interp_c(C);
      plot3(x_corrected(:, 1), x_corrected(:, 2), x_corrected(:, 3), 'Color', rgb,'Linewidth',1.5);
  end

  x_n_L2 = x02_ast;
  t_n_L2 = t02_ast;

  if count == p("count_max")
      break
  end
end
hold off

%% Stability color plot
hSt = figure();
hold on
grid on
box on
xlabel('$x$[-]');
ylabel('$y$[-]');
zlabel('$z$[-]');
axis square
view([-39 31]);
xlim([0.949274012415507 1.013359974485568]);
ylim([-0.039852272743994 0.035931818165097]);
% libration point
plot(L2(1), L2(2), '*', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
% Monodromy Matrix
M = zeros(6,6,p("count_max")/p("plot_interval"));
% Center manifold --altenatives of quasi-periodic tori
orbit_cent = cell(p("count_max")/p("plot_interval")+1,1);
for i=2:p("count_max")/p("plot_interval")+1
    Y0 = [x0_corrected(i,:),reshape(eye(6),1,[])];
    [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,p("mu")),[0 t0_corrected(i)],Y0',options_ODE);
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

%% arrange orbit data
orbit_data = zeros(p("count_max")/p("plot_interval")+1,6+2);
for i=1:p("count_max")/p("plot_interval")+1
  orbit_data(i,:) = [x0_corrected(i,:),t0_corrected(i),Jacobi(i)];
end
orbit_data = flipud(orbit_data);
%cd mat/
%save("inidata_torus_vertical_Lyapunov.mat","orbit_data");
%cd ..

%% SAVE RESULTS
cd fig/family/
%hvLf = strcat('vertical_Lyapunov_family_SE_mu=', num2str(p("mu")), '_x0=', num2str(x0_corrected(1,1)), '_y0dot=', num2str(x0_corrected(1,5)),'_zdot0=', num2str(x0_corrected(1,6)), '_t0=', num2str(t0_corrected(1)), '_scale=', num2str(p("scale")),'_count_max=', num2str(p("count_max")), '_plot_interval=', num2str(p("plot_interval")),'C_min=',num2str(p("Jacobi_min")),'C_max=',num2str(p("Jacobi_max")));
%hvLf = strrep(hvLf,'.',',');
%save_fig(hvL,hvLf,[-39 31]);
cd ..
cd ..

cd fig/family/
%hSt = strcat('vertical_Lyapunov_family_SE_mu=', num2str(p("mu")), '_x0=', num2str(x0_corrected(1,1)), '_y0dot=', num2str(x0_corrected(1,5)),'_zdot0=', num2str(x0_corrected(1,6)), '_t0=', num2str(t0_corrected(1)), '_scale=', num2str(p("scale")),'_count_max=', num2str(p("count_max")), '_plot_interval=', num2str(p("plot_interval")),'C_min=',num2str(p("Jacobi_min")),'C_max=',num2str(p("Jacobi_max")));
%hSt = strrep(hSt,'.',',');
%save_fig(hSt,hSt,[-39 31]);
cd ..
cd ..

%% End of script
time = strcat('calculation time: ', num2str(toc(myTimer)));
disp(time);
