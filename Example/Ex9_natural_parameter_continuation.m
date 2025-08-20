<<<<<<< HEAD
% natural parameter continuation for Lyapunov family in Sun-Earth CR3BP

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
addpath(replace(current_pass, 'Examples', 'Functions'));

myTimer = tic;        %start timer


%% initial settings
[mu, a_1, ~, ~]      = parameter(1); % Sun-Earth
[L1, L2, L3, L4, L5] = librationPoints(mu);
dim_L                = a_1; % [-] -> [km]

options_ODE   = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);
iteration_max = 100;
threshold     = 1e-10;
delta         = 2e-6;
count         = 0;
count_max     = 4000;
plot_interval = 200;

Jacobi_min = 3.000040465479221;
Jacobi_max = 3.000881506800496;

color      = jet;
Jacobi_lim = linspace(Jacobi_min, Jacobi_max, size(color,1));
Interp_c   = griddedInterpolant(Jacobi_lim, color);

Jacobi       = zeros(count_max/plot_interval, 1);
x0_corrected = zeros(count_max/plot_interval, 6);
t0_corrected = zeros(count_max/plot_interval, 1);


%% initial values
x0 = [1.01003168694194 0 0 0 8.61411407047802e-06 0]';
t0 = 1.47271346920646;


%% natural parameter continuation
f1 = figure();

while 1
	count = count + 1;
	disp(strcat('count = ', num2str(count)));

	% differential correction
	for iteration = 1:iteration_max
		[x_n, t_n, C] = fun_differential_correction_cr3bp(x0, t0, mu, options_ODE);

		tspan = [0 2*t_n];
		[t_corrected, x_corrected] = ode113(@(t, x) fun_cr3bp(t, x, mu), tspan, x_n, options_ODE);

		x_error = norm(x_corrected(end, :) - x_corrected(1, :));
		
    if x_error < threshold
			break;
    end
		
    if x_error > 1e+3
			disp('calculation diverged');
			return;
    end
		
    if iteration == iteration_max
			disp('do not finish');
			return;
    end

		x0 = x_n;
		t0 = t_n;
	end

	% plot orbit
	if (count >= plot_interval) && (mod(count, plot_interval) == 0)
		Jacobi(count/plot_interval, 1) = C;
		x0_corrected(count/plot_interval, :) = x_corrected(1, :);
		t0_corrected(count/plot_interval, 1) = t_corrected(end);
		rgb = Interp_c(C);
		plot(x_corrected(:, 1)*dim_L, x_corrected(:, 2)*dim_L, 'Color', rgb);
		hold on
	end

  x0(1) = x0(1) - delta;

	if count == count_max
		break
	end
end

f1_p1 = plot(L2(1)*dim_L, L2(2)*dim_L, '*', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
f1_p2 = plot((1-mu)*dim_L, 0, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
axis equal
colormap jet;
c = colorbar;
ylabel(c, 'Jacobi constant [-]', 'FontSize', 15);
caxis([Jacobi_min Jacobi_max]);
c.Ticks = linspace(Jacobi_min, Jacobi_max, 6);
xlabel('$x$[km]');
ylabel('$y$[km]');
grid on
legend([f1_p1, f1_p2], {'$L_2$', 'Earth'}, 'Location', 'northwest');
hold off


%% save results
f1_name = strcat('Ex9_natural_parameter_continuation_mu=', num2str(mu), '_x0=', num2str(x0(1)), '_ydot0=', num2str(x0(5)), '_t0=', num2str(t0),...
'_threshold=', num2str(threshold), '_delta=', num2str(delta), '_count_max=', num2str(count_max), '_plot_interval=', num2str(plot_interval));
f1_name = strrep(f1_name,'.',',');
save_fig(f1,f1_name,[0 90]);


%% End of script
time = strcat('calculation time: ', num2str(toc(myTimer)));
disp(time);
=======
% natural parameter continuation for Lyapunov family in Sun-Earth CR3BP

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
addpath(replace(current_pass, 'Examples', 'Functions'));

myTimer = tic;        %start timer


%% initial settings
[mu, a_1, ~, ~]      = parameter(1); % Sun-Earth
[L1, L2, L3, L4, L5] = librationPoints(mu);
dim_L                = a_1; % [-] -> [km]

options_ODE   = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);
iteration_max = 100;
threshold     = 1e-10;
delta         = 2e-6;
count         = 0;
count_max     = 4000;
plot_interval = 200;

Jacobi_min = 3.000040465479221;
Jacobi_max = 3.000881506800496;

color      = jet;
Jacobi_lim = linspace(Jacobi_min, Jacobi_max, size(color,1));
Interp_c   = griddedInterpolant(Jacobi_lim, color);

Jacobi       = zeros(count_max/plot_interval, 1);
x0_corrected = zeros(count_max/plot_interval, 6);
t0_corrected = zeros(count_max/plot_interval, 1);


%% initial values
x0 = [1.01003168694194 0 0 0 8.61411407047802e-06 0]';
t0 = 1.47271346920646;


%% natural parameter continuation
f1 = figure();

while 1
	count = count + 1;
	disp(strcat('count = ', num2str(count)));

	% differential correction
	for iteration = 1:iteration_max
		[x_n, t_n, C] = fun_differential_correction_cr3bp(x0, t0, mu, options_ODE);

		tspan = [0 2*t_n];
		[t_corrected, x_corrected] = ode113(@(t, x) fun_cr3bp(t, x, mu), tspan, x_n, options_ODE);

		x_error = norm(x_corrected(end, :) - x_corrected(1, :));
		
    if x_error < threshold
			break;
    end
		
    if x_error > 1e+3
			disp('calculation diverged');
			return;
    end
		
    if iteration == iteration_max
			disp('do not finish');
			return;
    end

		x0 = x_n;
		t0 = t_n;
	end

	% plot orbit
	if (count >= plot_interval) && (mod(count, plot_interval) == 0)
		Jacobi(count/plot_interval, 1) = C;
		x0_corrected(count/plot_interval, :) = x_corrected(1, :);
		t0_corrected(count/plot_interval, 1) = t_corrected(end);
		rgb = Interp_c(C);
		plot(x_corrected(:, 1)*dim_L, x_corrected(:, 2)*dim_L, 'Color', rgb);
		hold on
	end

  x0(1) = x0(1) - delta;

	if count == count_max
		break
	end
end

f1_p1 = plot(L2(1)*dim_L, L2(2)*dim_L, '*', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
f1_p2 = plot((1-mu)*dim_L, 0, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
axis equal
colormap jet;
c = colorbar;
ylabel(c, 'Jacobi constant [-]', 'FontSize', 15);
caxis([Jacobi_min Jacobi_max]);
c.Ticks = linspace(Jacobi_min, Jacobi_max, 6);
xlabel('$x$[km]');
ylabel('$y$[km]');
grid on
legend([f1_p1, f1_p2], {'$L_2$', 'Earth'}, 'Location', 'northwest');
hold off


%% save results
f1_name = strcat('Ex9_natural_parameter_continuation_mu=', num2str(mu), '_x0=', num2str(x0(1)), '_ydot0=', num2str(x0(5)), '_t0=', num2str(t0),...
'_threshold=', num2str(threshold), '_delta=', num2str(delta), '_count_max=', num2str(count_max), '_plot_interval=', num2str(plot_interval));
f1_name = strrep(f1_name,'.',',');
save_fig(f1,f1_name,[0 90]);


%% End of script
time = strcat('calculation time: ', num2str(toc(myTimer)));
disp(time);
>>>>>>> 748834e2ef04e45bd3c365491093c08a43a6640e
