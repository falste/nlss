% This example simulates a forced Van-der-Pol oscillator and defines an output.
% https://en.wikipedia.org/wiki/Van_der_Pol_oscillator

%% Init
clear;
syms t;					% Symbolic time variable
x = sym('x', [2, 1]);	% State vector variable
u = sym('u', [1, 1]);	% Input vector variable

%% Settings
x0 = [1; 0];			% Initial condition
T  = 20;				% Simulation time

%% Create and simulate system
a = 0.5;
f = [x(2)                                ; ...
	 a * (1-x(1)^2) * x(2) - x(1) + u(1)];
g = [x(1)*x(2)    ; ...
	(x(1)-x(2))^2];

sys = NLSS(t, x, f, u, g);
sys.sim(x0, T, pi*cos(t));

% Alternative:
% [t_sim, x_sim, y_sim] = sys.sim(...);
