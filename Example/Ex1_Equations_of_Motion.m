% calculate trajectories 

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



%% initial conditions
Ms = 1.9885e+30; %[kg], mass of Sun, from https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
M1 = 5.9724e+24; %[kg], mass of Earth, from https://nssdc.gsfc.nasa.gov/planetary/factsheet/
M2 = 7.346e+22; %[kg], mass of Moon, from https://nssdc.gsfc.nasa.gov/planetary/factsheet/
G = 6.673e-11; %[m^3/(kg*s^2)], gravitational constant

% Earth-Moon CR3BP
chara_length_CR3BP = 3.8440e+8; %[m]
chara_mass_CR3BP = M1 + M2; %[kg]
chara_time_CR3BP = sqrt( chara_length_CR3BP^3 / (G * chara_mass_CR3BP) ); % == period / (2*pi)
mu_EM = M2 / (M1 + M2);
N_CR3BP = sqrt( G*chara_mass_CR3BP / chara_length_CR3BP^3 );

% Earth-Moon ER3BP
a = 3.8440e+8; %[m], semimajor axis, from https://nssdc.gsfc.nasa.gov/planetary/factsheet/
N_ER3BP = sqrt( G*(M1+M2)/(a^3) ); %[1/s], mean motion
e = 0.0549; % eccentricity, from https://nssdc.gsfc.nasa.gov/planetary/factsheet/

% Earth-Moon-Sun BCR4BP
chara_length_SB1 = 1.4960e+11; %[m]
chara_mass_SB1 = M1 + M2 + Ms; %[kg]
chara_time_SB1 = sqrt( chara_length_SB1^3 / (G * chara_mass_SB1) ); % == period / (2*pi)
mS = Ms / (M1 + M2); % nondimensional mass of Sun
aS = chara_length_SB1 / chara_length_CR3BP; % nondimensional Sun orbit radius
wS = sqrt( (1+mS) / (aS^3) ) - 1; % nondimensional Sun angular velocity
thetaS0 = pi; %[rad], when a total eclipse of the moon occurs
mu_SB1 = (M1 + M2) / (M1 + M2 + Ms); % mass ratio in Sun - B1 frame
aEM = chara_length_CR3BP / chara_length_SB1; % nondimensional distance between Earth and Moon
wM = chara_time_SB1/chara_time_CR3BP - 1; % nondimensional Moon angular velocity
thetaM0 = 0; %[rad], when a total eclipse of the sun occurs

% calculations
options_ODE = odeset('RelTol',1e-14, 'AbsTol',1e-14);


%% L2 Lyapunov orbit in Earth-Moon CR3BP
load('estimated_initial_values.mat','x_n_CR3BP','t_n_CR3BP');
t_CR3BP = t_n_CR3BP;
x0_CR3BP = x_n_CR3BP;
tspan1 = [0 t_CR3BP];
[~, x_CR3BP] = ode113(@(t,x) fun_cr3bp(t,x,mu_EM), tspan1, x0_CR3BP, options_ODE);

f1 = figure;
hold on
p_moon = plot(1-mu_EM, 0, 'o', 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120', 'MarkerSize', 20);
plot(x_CR3BP(1,1), x_CR3BP(1,2),'o', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', '#0072BD', 'MarkerSize', 4);
plot(x_CR3BP(:,1), x_CR3BP(:,2), 'Color', '#0072BD');
axis equal
xlabel('$x$[-]');
ylabel('$y$[-]');
legend([p_moon], {'Moon'});
hold off
f1_name = strcat('Ex1_L2_Lyapunov_orbit_CR3BP_x0=', num2str(x0_CR3BP(1)), '_vy0=', num2str(x0_CR3BP(5)), '_t=', num2str(t_CR3BP));
f1_name = strrep(f1_name,'.',',');
save_fig(f1,f1_name,[0 90]);



%% trajectory in Earth-Moon ER3BP
t_ER3BP = t_CR3BP * chara_time_CR3BP * N_ER3BP;

% transfomation from CR3BP reference frame to the inertial frame
X0_CR3BP = [x0_CR3BP(1:3).*chara_length_CR3BP; x0_CR3BP(4:6).*chara_length_CR3BP./chara_time_CR3BP];
C_CR3BP = [cos(0) -sin(0) 0;
           sin(0)  cos(0) 0;
                0       0 1];
dtheta_dt_CR3BP = N_CR3BP;
rotating_matrix_CR3BP = fun_rotating_to_inertial_matrix(C_CR3BP, dtheta_dt_CR3BP);
X0_inertial = rotating_matrix_CR3BP * X0_CR3BP;

% transformation from the inertial frame to ER3BP reference frame
C_ER3BP = [cos(0) -sin(0) 0;
           sin(0)  cos(0) 0;
                0       0 1];
dtheta_dt_ER3BP = sqrt( G*(M1+M2) * (1+e*cos(0))^4 / (a*(1-e^2))^3 );
rotating_matrix_ER3BP = fun_rotating_to_inertial_matrix(C_ER3BP, dtheta_dt_ER3BP);
X0_ER3BP = rotating_matrix_ER3BP \ X0_inertial;
x0_ER3BP = [X0_ER3BP(1:3)./(a*(1-e)); X0_ER3BP(4:6)./(a*(1-e))./N_ER3BP];

tspan2 = [0 t_ER3BP];
[~, x_ER3BP] = ode113(@(t,x) fun_er3bp(t,x,mu_EM,e), tspan2, x0_ER3BP, options_ODE);

f2 = figure;
hold on
p_moon = plot(1-mu_EM, 0, 'o', 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120', 'MarkerSize', 20);
plot(x_ER3BP(1,1), x_ER3BP(1,2),'o', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', '#0072BD', 'MarkerSize', 4);
plot(x_ER3BP(:,1), x_ER3BP(:,2), 'Color', '#0072BD');
axis equal
xlabel('$x$[-]');
ylabel('$y$[-]');
legend([p_moon], {'Moon'});
hold off
f2_name = strcat('Ex1_ER3BP_x0=', num2str(x0_ER3BP(1)), '_vy0=', num2str(x0_ER3BP(5)), '_t=', num2str(t_ER3BP));
f2_name = strrep(f2_name,'.',',');
save_fig(f2,f2_name,[0 90]);



%% trajectory in Earth-Moon_Sun BCR4BP in the Earth-Moon rotating frame
tspan3 = [0 t_CR3BP];
[~, x_bcr4bp_EMS] = ode113(@(t,x) fun_bcr4bp_EMS(t,x,mu_EM,mS,aS,wS,thetaS0), tspan3, x0_CR3BP, options_ODE);

f3 = figure;
hold on
p_moon = plot(1-mu_EM, 0, 'o', 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120', 'MarkerSize', 20);
plot(x_bcr4bp_EMS(1,1), x_bcr4bp_EMS(1,2),'o', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', '#0072BD', 'MarkerSize', 4);
plot(x_bcr4bp_EMS(:,1), x_bcr4bp_EMS(:,2), 'Color', '#0072BD');
axis equal
xlabel('$x$[-]');
ylabel('$y$[-]');
legend([p_moon], {'Moon'});
hold off
f3_name = strcat('Ex1_bcr4bp_EMS_x0=', num2str(x0_CR3BP(1)), '_vy0=', num2str(x0_CR3BP(5)), '_t=', num2str(t_CR3BP));
f3_name = strrep(f3_name,'.',',');
save_fig(f3,f3_name,[0 90]);



%% trajectory in Earth-Moon_Sun BCR4BP in the Sun-B1 rotating frame
theta = linspace(0, 2*pi, 1001)';
x_moon = 1 - mu_SB1 + aEM*(1 - mu_EM)*cos(theta);
y_moon = aEM*(1 - mu_EM)*sin(theta);
x_Earth = 1 - mu_SB1 + aEM*mu_EM*cos(theta);
y_Earth = aEM*mu_EM*sin(theta);

t_bcr4bp_SB1 = t_CR3BP * chara_time_CR3BP / chara_time_SB1;

% transformation from the intertial frame to SB1 reference frame
C_SB1 = [cos(0) -sin(0) 0;
         sin(0)  cos(0) 0;
              0       0 1];
dtheta_dt_SB1 = 1/chara_time_SB1;
rotating_matrix_SB1 = fun_rotating_to_inertial_matrix(C_SB1, dtheta_dt_SB1);
X0_bcr4bp_SB1 = rotating_matrix_SB1 \ X0_inertial;
x0_bcr4bp_SB1 = [X0_bcr4bp_SB1(1:3)./chara_length_SB1; X0_bcr4bp_SB1(4:6)./chara_length_SB1.*chara_time_SB1];

% parallel translation
x0_bcr4bp_SB1(1) = x0_bcr4bp_SB1(1) + (1 - mu_SB1);

tspan4 = [0 t_bcr4bp_SB1];
[~, x_bcr4bp_EMS] = ode113(@(t,x) fun_bcr4bp_SB1(t,x,mu_SB1,mu_EM,1,aEM,wM,thetaM0), tspan4, x0_bcr4bp_SB1, options_ODE);

f4 = figure;
hold on
p_moon = plot(x_moon, y_moon, '-', 'Color', '#EDB120');
p_Earth = plot(x_Earth, y_Earth, '-', 'Color', '#7E2F8E');
plot(x_bcr4bp_EMS(1,1), x_bcr4bp_EMS(1,2),'o', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', '#0072BD', 'MarkerSize', 4);
plot(x_bcr4bp_EMS(:,1), x_bcr4bp_EMS(:,2), 'Color', '#0072BD');
axis equal
xlabel('$x$[-]');
ylabel('$y$[-]');
legend([p_moon, p_Earth], {'Moon','Earth'});
hold off
f4_name = strcat('Ex1_bcr4bp_SB1_x0=', num2str(x0_bcr4bp_SB1(1)), '_vy0=', num2str(x0_bcr4bp_SB1(5)), '_t=', num2str(t_bcr4bp_SB1));
f4_name = strrep(f4_name,'.',',');
save_fig(f4,f4_name,[0 90]);



%% End of script
time = strcat('calculation time: ', num2str(toc(myTimer)));
disp(time);
