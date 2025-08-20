<<<<<<< HEAD
% pseudo arclength continuation for halo family in Sun-Earth CR3BP

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
scale         = 1e-4;
count         = 0;
count_max     = 8400;
plot_interval = 100;

Jacobi_min = 3.000212571920089;
Jacobi_max = 3.000818466710653;

xlim_min = 1.485*10^8;
xlim_max = 1.525*10^8;
ylim_min = -2.5*10^6;
ylim_max = 2.5*10^6;
zlim_min = -23*10^5;
zlim_max = 10*10^5;

color      = jet;
Jacobi_lim = linspace(Jacobi_min, Jacobi_max, size(color,1));
Interp_c   = griddedInterpolant(Jacobi_lim, color);

Jacobi       = zeros(count_max/plot_interval+1, 1);
x0_corrected = zeros(count_max/plot_interval+1, 6);
t0_corrected = zeros(count_max/plot_interval+1, 1);


%% initial values
x01_ast = [1.00837693674644 0 0.000200000000000000 0 0.00976542047036519 0]';
t01_ast = 1.55123023270505;


%% pseudo arclength continuation
f1 = figure();
hold on

tspan1  = [0, 2*t01_ast];
[~, x1] = ode113(@(t, x) fun_cr3bp(t, x, mu), tspan1, x01_ast, options_ODE);

Jacobi(1)          = Jacobi_const(x01_ast, mu);
x0_corrected(1, :) = x01_ast';
t0_corrected(1)    = t01_ast;
rgb = Interp_c(Jacobi(1));
plot3(x1(:, 1)*dim_L, x1(:, 2)*dim_L, x1(:, 3)*dim_L, 'Color', rgb);
plot3(x1(:, 1)*dim_L, x1(:, 2)*dim_L, zlim_min*ones(numel(x1(:, 1)), 1), 'Color', (rgb + [1, 1, 1])./2);
plot3(x1(:, 1)*dim_L, ylim_max*ones(numel(x1(:, 1)), 1), x1(:, 3)*dim_L, 'Color', (rgb + [1, 1, 1])./2);
plot3(xlim_max*ones(numel(x1(:, 1)), 1), x1(: ,2)*dim_L, x1(:, 3)*dim_L, 'Color', (rgb + [1, 1, 1])./2);

while 1
	count = count + 1;
	disp(strcat('count = ', num2str(count)));

  [delta, ~] = fun_null_cr3bp(x01_ast, t01_ast, mu, options_ODE);
  delta      = -delta;
  x02        = x01_ast;
  x02(1)     = x01_ast(1) + scale*delta(1);
  x02(3)     = x01_ast(3) + scale*delta(2);
  x02(5)     = x01_ast(5) + scale*delta(3);
  t02        = t01_ast + scale*delta(4);

	% differential correction
	for iteration = 1:iteration_max
		[x02_ast, t02_ast, C, G] = fun_differential_correction_cr3bp_PAC(x01_ast, t01_ast, x02, t02, scale, delta, mu, options_ODE);

		tspan2 = [0, 2*t02_ast];
		[t_corrected, x_corrected] = ode113(@(t, x) fun_cr3bp(t, x, mu), tspan2, x02_ast, options_ODE);

		if norm(G) < threshold
			break;
		end
		
    if norm(G) > 1e+3
			disp('calculation diverged');
			return;
    end
		
    if iteration == iteration_max
			disp('do not finish');
			return;
    end

		x02 = x02_ast;
		t02 = t02_ast;
	end

  % plot orbit
	if (count >= plot_interval) && (mod(count, plot_interval) == 0)
		Jacobi(count/plot_interval+1)          = C;
		x0_corrected(count/plot_interval+1, :) = x_corrected(1, :);
		t0_corrected(count/plot_interval+1)    = t_corrected(end);
		rgb = Interp_c(C);
		plot3(x_corrected(:, 1)*dim_L, x_corrected(:, 2)*dim_L, x_corrected(:, 3)*dim_L, 'Color', rgb);
		plot3(x_corrected(:, 1)*dim_L, x_corrected(:, 2)*dim_L, zlim_min*ones(numel(x_corrected(:, 1)), 1), 'Color', (rgb + [1, 1, 1])./2);
    plot3(x_corrected(:, 1)*dim_L, ylim_max*ones(numel(x_corrected(:, 1)), 1), x_corrected(:, 3)*dim_L, 'Color', (rgb + [1, 1, 1])./2);
    plot3(xlim_max*ones(numel(x_corrected(:, 1)), 1), x_corrected(: ,2)*dim_L, x_corrected(:, 3)*dim_L, 'Color', (rgb + [1, 1, 1])./2);
	end

	x01_ast = x02_ast;
	t01_ast = t02_ast;

  if count == count_max
		break
	end
end

f1_p1 = plot3(L2(1)*dim_L, L2(2)*dim_L, L2(3)*dim_L, '*', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
f1_p2 = plot3((1-mu)*dim_L, 0, 0, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
view_angle = [-30 10];
view(gca, view_angle(1), view_angle(2));
axis equal
colormap jet;
c = colorbar;
ylabel(c, 'Jacobi constant [-]', 'FontSize', 15);
caxis([Jacobi_min Jacobi_max]);
c.Ticks = linspace(Jacobi_min, Jacobi_max, 6);
xlim([xlim_min, xlim_max]);
ylim([ylim_min, ylim_max]);
zlim([zlim_min, zlim_max]);
xlabel('$x$[km]');
ylabel('$y$[km]');
zlabel('$z$[km]');
grid on
legend([f1_p1, f1_p2], {'$L_2$', 'Earth'}, 'Location', 'northwest');
hold off


%% save results
f1_name = strcat('Ex10_pseudo_arclength_continuation_mu=', num2str(mu), '_x0=', num2str(x0_corrected(1,1)), '_z0=', num2str(x0_corrected(1,3)),...
'_ydot0=', num2str(x0_corrected(1,5)), '_t0=', num2str(t0_corrected(1)), '_threshold=', num2str(threshold), '_scale=', num2str(scale),...
'_count_max=', num2str(count_max), '_plot_interval=', num2str(plot_interval));
f1_name = strrep(f1_name,'.',',');
save_fig(f1,f1_name,[0 90]);


%% End of script
time = strcat('calculation time: ', num2str(toc(myTimer)));
=======
% pseudo arclength continuation for halo family in Sun-Earth CR3BP

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
scale         = 1e-4;
count         = 0;
count_max     = 8400;
plot_interval = 100;

Jacobi_min = 3.000212571920089;
Jacobi_max = 3.000818466710653;

xlim_min = 1.485*10^8;
xlim_max = 1.525*10^8;
ylim_min = -2.5*10^6;
ylim_max = 2.5*10^6;
zlim_min = -23*10^5;
zlim_max = 10*10^5;

color      = jet;
Jacobi_lim = linspace(Jacobi_min, Jacobi_max, size(color,1));
Interp_c   = griddedInterpolant(Jacobi_lim, color);

Jacobi       = zeros(count_max/plot_interval+1, 1);
x0_corrected = zeros(count_max/plot_interval+1, 6);
t0_corrected = zeros(count_max/plot_interval+1, 1);


%% initial values
x01_ast = [1.00837693674644 0 0.000200000000000000 0 0.00976542047036519 0]';
t01_ast = 1.55123023270505;


%% pseudo arclength continuation
f1 = figure();
hold on

tspan1  = [0, 2*t01_ast];
[~, x1] = ode113(@(t, x) fun_cr3bp(t, x, mu), tspan1, x01_ast, options_ODE);

Jacobi(1)          = Jacobi_const(x01_ast, mu);
x0_corrected(1, :) = x01_ast';
t0_corrected(1)    = t01_ast;
rgb = Interp_c(Jacobi(1));
plot3(x1(:, 1)*dim_L, x1(:, 2)*dim_L, x1(:, 3)*dim_L, 'Color', rgb);
plot3(x1(:, 1)*dim_L, x1(:, 2)*dim_L, zlim_min*ones(numel(x1(:, 1)), 1), 'Color', (rgb + [1, 1, 1])./2);
plot3(x1(:, 1)*dim_L, ylim_max*ones(numel(x1(:, 1)), 1), x1(:, 3)*dim_L, 'Color', (rgb + [1, 1, 1])./2);
plot3(xlim_max*ones(numel(x1(:, 1)), 1), x1(: ,2)*dim_L, x1(:, 3)*dim_L, 'Color', (rgb + [1, 1, 1])./2);

while 1
	count = count + 1;
	disp(strcat('count = ', num2str(count)));

  [delta, ~] = fun_null_cr3bp(x01_ast, t01_ast, mu, options_ODE);
  delta      = -delta;
  x02        = x01_ast;
  x02(1)     = x01_ast(1) + scale*delta(1);
  x02(3)     = x01_ast(3) + scale*delta(2);
  x02(5)     = x01_ast(5) + scale*delta(3);
  t02        = t01_ast + scale*delta(4);

	% differential correction
	for iteration = 1:iteration_max
		[x02_ast, t02_ast, C, G] = fun_differential_correction_cr3bp_PAC(x01_ast, t01_ast, x02, t02, scale, delta, mu, options_ODE);

		tspan2 = [0, 2*t02_ast];
		[t_corrected, x_corrected] = ode113(@(t, x) fun_cr3bp(t, x, mu), tspan2, x02_ast, options_ODE);

		if norm(G) < threshold
			break;
		end
		
    if norm(G) > 1e+3
			disp('calculation diverged');
			return;
    end
		
    if iteration == iteration_max
			disp('do not finish');
			return;
    end

		x02 = x02_ast;
		t02 = t02_ast;
	end

  % plot orbit
	if (count >= plot_interval) && (mod(count, plot_interval) == 0)
		Jacobi(count/plot_interval+1)          = C;
		x0_corrected(count/plot_interval+1, :) = x_corrected(1, :);
		t0_corrected(count/plot_interval+1)    = t_corrected(end);
		rgb = Interp_c(C);
		plot3(x_corrected(:, 1)*dim_L, x_corrected(:, 2)*dim_L, x_corrected(:, 3)*dim_L, 'Color', rgb);
		plot3(x_corrected(:, 1)*dim_L, x_corrected(:, 2)*dim_L, zlim_min*ones(numel(x_corrected(:, 1)), 1), 'Color', (rgb + [1, 1, 1])./2);
    plot3(x_corrected(:, 1)*dim_L, ylim_max*ones(numel(x_corrected(:, 1)), 1), x_corrected(:, 3)*dim_L, 'Color', (rgb + [1, 1, 1])./2);
    plot3(xlim_max*ones(numel(x_corrected(:, 1)), 1), x_corrected(: ,2)*dim_L, x_corrected(:, 3)*dim_L, 'Color', (rgb + [1, 1, 1])./2);
	end

	x01_ast = x02_ast;
	t01_ast = t02_ast;

  if count == count_max
		break
	end
end

f1_p1 = plot3(L2(1)*dim_L, L2(2)*dim_L, L2(3)*dim_L, '*', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
f1_p2 = plot3((1-mu)*dim_L, 0, 0, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
view_angle = [-30 10];
view(gca, view_angle(1), view_angle(2));
axis equal
colormap jet;
c = colorbar;
ylabel(c, 'Jacobi constant [-]', 'FontSize', 15);
caxis([Jacobi_min Jacobi_max]);
c.Ticks = linspace(Jacobi_min, Jacobi_max, 6);
xlim([xlim_min, xlim_max]);
ylim([ylim_min, ylim_max]);
zlim([zlim_min, zlim_max]);
xlabel('$x$[km]');
ylabel('$y$[km]');
zlabel('$z$[km]');
grid on
legend([f1_p1, f1_p2], {'$L_2$', 'Earth'}, 'Location', 'northwest');
hold off


%% save results
f1_name = strcat('Ex10_pseudo_arclength_continuation_mu=', num2str(mu), '_x0=', num2str(x0_corrected(1,1)), '_z0=', num2str(x0_corrected(1,3)),...
'_ydot0=', num2str(x0_corrected(1,5)), '_t0=', num2str(t0_corrected(1)), '_threshold=', num2str(threshold), '_scale=', num2str(scale),...
'_count_max=', num2str(count_max), '_plot_interval=', num2str(plot_interval));
f1_name = strrep(f1_name,'.',',');
save_fig(f1,f1_name,[0 90]);


%% End of script
time = strcat('calculation time: ', num2str(toc(myTimer)));
>>>>>>> 748834e2ef04e45bd3c365491093c08a43a6640e
disp(time);