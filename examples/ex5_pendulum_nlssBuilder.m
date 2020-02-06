clear;

%% Setup
x0 = [pi/2; 0];
T = 5;

%% System Definition
syms m l g t Mu c_f alpha(t) y1(t);
constants = [m == 1; l == 0.3; g == 9.81; c_f == 0.3];


nlssBuilder = NLSSBuilder(t);

% Define an axis of motion using a variable for that motion, the inertia on that
% axis, the generalized input force along that axis and the friction.
nlssBuilder.AddAxis(alpha, m * l^2, Mu, c_f);

% The kinetic energy in the whole system.
nlssBuilder.V = -m*g*l*cos(alpha);

% Display information about the current nlssBuilder state.
nlssBuilder.Show();

% Create a nlss. Internally, this uses the lagrange-formalism.
sys = nlssBuilder.GenerateNLSS(Mu, [y1 == alpha], constants);

%% Simulation
sys.sim(x0, T, 0);