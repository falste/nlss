% This example simulates a simple, unforced Van-der-Pol oscillator using a
% nonlinear state-space model.
% https://en.wikipedia.org/wiki/Van_der_Pol_oscillator

%% Init
clear;
syms t;					% Symbolic time variable
x = sym('x', [2, 1]);	% State vector variable

%% Settings
x0 = [1; 0];			% Initial condition
T  = 20;				% Simulation time

%% Create and simulate system
a = 0.5;
f = [x(2)                         ; ...
	 a * (1-x(1)^2) * x(2) - x(1)];

sys = NLSS(t, x, f);
sys.sim(x0, T);

% Alternative:
% [t_sim, x_sim, y_sim] = sys.sim(...);