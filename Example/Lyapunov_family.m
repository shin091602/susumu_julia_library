%% Natural parameter continuation for Lyapunov family around L2 in Sun-Earth CR3BP
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

p = dictionary();

%% DICTIONARY OF SUN-EARTH CRTBP VARIABLES
p("Ms") = 1.9885e+30;
p("M1") = 5.9724e+24;
p("mu") = 3.00346093142064e-06;%SE
p("G") = 6.67300000000000e-11;
p("iteration_limit") = 50;
p("threshold") = 1e-10;
p("xpert") = 1e-6;
p("tf_L2") = 10;
p("scale") = 5e-6;
p("Bscale") = 1e-3;
p("tol") = 1e3;
count = 0;
p("count_max")     = 1500;
p("plot_interval") = 30;
p("Jacobi_min") = 3.0000038403545658;
p("Jacobi_max") = 3.0008942050035974;

%% INITIALIZATION
color      = jet;
Jacobi_lim = linspace(p("Jacobi_min"), p("Jacobi_max"), size(color,1));
Interp_c   = griddedInterpolant(Jacobi_lim, color);
Jacobi       = zeros(p("count_max")/p("plot_interval")+1, 1);
x0_corrected = zeros(p("count_max")/p("plot_interval")+1, 6);
t0_corrected = zeros(p("count_max")/p("plot_interval")+1, 1);

%% LIBRATION POINT L2
[~,L2,~,~,~] = librationPoints(p("mu"));

%% OPTIONS ODE
options_ODE = odeset('RelTol',1e-13, 'AbsTol',1e-13);

%% INITIAL ORBIT
%L2
x0 = [1.0172140117282569;0;0;0;-3.4905740101734528e-2;0];
t0 = 5.5015410088953178/2;
%L1
%x0 = [9.9189564755779780e-1;0;0;0;-1.1474478284680680e-2;0];
%t0 = 3.0806424724612258;

%orbital period
t_n_L2 = t0;
x_n_L2 = x0;

%% PLOT LYAPUNOV FAMILY
hLy = figure();
hold on
grid on
axis equal
xlabel('$x$[-]');
ylabel('$y$[-]');
xlim([0.989201438045879 1.020798561954121]);
ylim([-0.026332196721286 0.026329676459118]);

% libration point
plot(L2(1), L2(2), '*', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);

% NPC toward the libration point direction
while 1
    count = count + 1;
    disp(strcat('count = ', num2str(count)));

    %differential correction
    for iteration=1:p("iteration_limit")
        [x_n,t_n,C] = fun_differential_correction_cr3bp(x0,t0,p("mu"),options_ODE);
        tspan = [0 2*t_n];
        [t_corrected, x_corrected] = ode113(@(t,x) fun_cr3bp(t,x,p("mu")), tspan, x_n, options_ODE);

        x_error = norm(x_corrected(end,:) - x_corrected(1,:));
        if x_error < p("threshold")
            break;
        end

        if x_error > p("tol")
            disp('calculation diverged');
            return;
        end

        if iteration == p("iteration_limit")
            disp('do not finish');
            return;
        end

        x0 = x_n;
        t0 = t_n;
    end

    x0 = x_n;
    t0 = t_n;
    T = 2*t0;

    % plot orbit
    if (count >= p("plot_interval")) && (mod(count, p("plot_interval")) == 0)
        if C>p("Jacobi_max")
            disp("The value of Jacobi constant reaches the upper limit.");
            break;
        end
        Jacobi(count/p("plot_interval")+1)          = C;
        x0_corrected(count/p("plot_interval")+1, :) = x0;
        t0_corrected(count/p("plot_interval")+1)    = T;
        rgb = Interp_c(C);
        plot(x_corrected(:, 1), x_corrected(:, 2), 'Color', rgb,'LineWidth',1.5);
    end

    if count == p("count_max")
        break
    end

    % NPC
    x0(1) = x0(1)-p("scale");
    x0(5) = x0(5)-p("scale");
end

%color setting
colormap jet;
c = colorbar;
ylabel(c, 'Jacobi constant [-]', 'FontSize', 15);
clim([p("Jacobi_min") p("Jacobi_max")]);
c.Ticks = linspace(p("Jacobi_min"), p("Jacobi_max"), 6);
hold off

%% Stability color plot
hSt = figure();
hold on
grid on
axis equal
xlabel('$x$[-]');
ylabel('$y$[-]');
xlim([0.989201438045879 1.020798561954121]);
ylim([-0.026332196721286 0.026329676459118]);
% libration point
plot(L2(1), L2(2), '*', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
% Monodromy Matrix
M = zeros(6,6,p("count_max")/p("plot_interval"));
% Center manifold --altenatives of quasi-periodic tori
orbit_cent = cell(p("count_max")/p("plot_interval")+1,1);
% fifth component -- eigenvalues
eigen5 = zeros(length(orbit_cent),1);
for i=2:p("count_max")/p("plot_interval")+1
    Y0 = [x0_corrected(i,:),reshape(eye(6),1,[])];
    [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,p("mu")),[0 t0_corrected(i)],Y0',options_ODE);
    M(:,:,i) = reshape(Y(end,7:42),6,6);
    % eigenvectors and eignvalues
    [V,D] = eig(M(:,:,i));
    % stability
    if (imag(D(5,5))~=0)&&(abs(abs(D(5,5))-1)<1e-3)
        fcen = plot(Y(:,1),Y(:,2),'b','Linewidth',1.5);
        % save initial state
        orbit_cent{i,1} = Y(1,1:6);
    else
        % corresponds to center direciton
        funcen = plot(Y(:,1),Y(:,2),'r','Linewidth',1.5);
    end
end
hold off

%% Bifurcation
hB = figure();
hold on
grid on
axis equal
xlabel('$Re$[-]');
ylabel('$Im$[-]');
xlim([-1.5 2.0]);
ylim([-1.5 1.5]);

% circle
th0 = linspace(0,2*pi,360);
plot(cos(th0),sin(th0),'k','LineWidth',1.5);

for i=2:p("count_max")/p("plot_interval")+1
    Y0 = [x0_corrected(i,:),reshape(eye(6),1,[])];
    [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,p("mu")),[0 t0_corrected(i)],Y0',options_ODE);
    M(:,:,i) = reshape(Y(end,7:42),6,6);
    % eigenvectors and eignvalues
    [V,D] = eig(M(:,:,i));
    for j=1:6
        if abs(real(D(j,j))-1)<1e-9
            fbp = scatter(real(D(j,j)),imag(D(j,j)),'ro','filled');
        else
            fb = scatter(real(D(j,j)),imag(D(j,j)),'bo');
            hold on
        end
    end
end
legend(fbp,'a Bifurcarion eignevalue','Location','north');
hold  off

%% SAVE RESULTS
%{
cd fig/family/
h_name = strcat('Lyapunovfamily_SE_mu=', num2str(p("mu")), '_x0=', num2str(x0_corrected(1,1)), '_ydot0=', num2str(x0_corrected(1,5)), '_t0=', num2str(t0_corrected(1)),'_count_max=', num2str(p("count_max")), '_plot_interval=', num2str(p("plot_interval")),'C_min=',num2str(p("Jacobi_min")),'C_max=',num2str(p("Jacobi_max")));
h_name = strrep(h_name,'.',',');
save_fig(hLy,h_name,[0 90]);

h_name = strcat('Lyapunovfamily_stability_SE_mu=', num2str(p("mu")), '_x0=', num2str(x0_corrected(1,1)), '_ydot0=', num2str(x0_corrected(1,5)), '_t0=', num2str(t0_corrected(1)),  '_count_max=', num2str(p("count_max")), '_plot_interval=', num2str(p("plot_interval")),'C_min=',num2str(p("Jacobi_min")),'C_max=',num2str(p("Jacobi_max")));
h_name = strrep(h_name,'.',',');
save_fig(hSt,h_name,[0 90]);

cd ..
cd ..
%}

%% End of script
time = strcat('calculation time: ', num2str(toc(myTimer)));
disp(time);
