% This example simulates two systems connected in series. The first system is
% the system from example 3 and the second system is a linear state-space model.

clear;

%% Create first system
syms t;					% Symbolic time variable
x = sym('x', [2, 1]);	% State vector variable
u = sym('u', [1, 1]);	% Input vector variable

a = 0.5;
f = [x(2)                                ; ...
	 a * (1-x(1)^2) * x(2) - x(1) + u(1)];
g = [x(1)*x(2)    ; ...
	(x(1)-x(2))^2];

sys1 = NLSS(t, x, f, u, g);

%% Create second system
A = [-1 0; 1 -3];
B = [1  0; 0  1];
C = [1 0];
D = [0 0];

sys2 = ss(A, B, C, D);

% Also possible, but does not change anything here:
% sys2 = NLSS(ss(A, B, C, D));

%% Connect in series and simulate

combinedSys = sys1 * sys2;

x0 = [1; 0; 0; 0];			% Initial condition
T  = 20;					% Simulation time
combinedSys.sim(x0, T, pi*cos(t));












