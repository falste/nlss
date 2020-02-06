classdef NLSSBuilder < handle
% Helper class for creating a nonlinear state-space system (nlss) from
% generalized coordinates and forces along those coordinates. NLSSBuilder is
% also able to handle proportional friction along the generalized coordinates.
	
	properties (Access = private)
		coordData = cell(0, 5);
		t;
	end

	properties (Dependent)
		q;			% Generalized coordinates
		Q_original; % Generalized forces
		c_f;		% Friction coefficient
		Q;			% Generalized forces + friction force along that axis
		inertia;	
		N;			% Number of generalized coordinates
		T;			% Kinetic energy in the system
	end
	
	properties (SetAccess = private)
		substitutions;
	end
	
	properties
		V;
	end

	methods
		function obj = NLSSBuilder(t)
			obj.t = t;
		end
		
		function AddAxis(obj, coord, inertia, force, frictionCoefficient)
			resForce = force - frictionCoefficient * diff(coord, obj.t);
			
			obj.coordData(obj.N+1, :) = {coord, force, frictionCoefficient, resForce, inertia};
		end
		
		function Show(obj)
			fprintf('\t<strong>Coordinate Data:</strong>\n\n');
			fricchar = obj.c_f.';
			for i = 1:length(obj.c_f)
				if isa(fricchar(i), 'double')
					fricchar(i) = num2str(fricchar(i));
				else
					fricchar(i) = char(fricchar(i));
				end
			end
			disp(table( ...
				arrayfun(@char, obj.q.', 'uniform', 0), ...
				arrayfun(@char, obj.inertia.', 'uniform', 0), ...
				arrayfun(@char, obj.Q_original.', 'uniform', 0), ...
				arrayfun(@char, fricchar, 'uniform', 0), ...
				arrayfun(@char, obj.Q.', 'uniform', 0), ...
				'VariableNames',{'Coordinate', 'Inertia', 'Force', 'Friction', 'Force_sum'}));
		end
		
		function nlss = GenerateNLSS(obj, u, outputEqns, constants)
			% Get movement equations via lagrange formalism
			moveEquations = lagrange(obj.q, obj.Q, obj.T, obj.V);
			
			% Turn movement equations into state-space model
			if nargin >= 3
				[x, f, g, obj.substitutions] = moveEquation2ss(obj.q, moveEquations, rhs(outputEqns));
			else
				[x, f, g, obj.substitutions] = moveEquation2ss(obj.q, moveEquations);
			end
			x = subs(x, lhs(constants), rhs(constants));
			u = subs(u, lhs(constants), rhs(constants));
			f = subs(f, lhs(constants), rhs(constants));
			g = subs(g, lhs(constants), rhs(constants));
			nlss = NLSS(obj.t, x, f, u, g);
			
			% Set variable labels of the state-space model (for plotting)
			nlss.xlabels = cell(length(x),1);
			for i = 1:length(x)
				nlss.xlabels{i} = texlabel(lhs(obj.substitutions(i)));
			end
			inputElements = formula(u);
			for i = 1:length(inputElements)
				nlss.ulabels{i} = texlabel(inputElements(i));
			end
			if nargin >= 3
				outputElements = formula(lhs(outputEqns));
				for i = 1:length(outputElements)
					nlss.ylabels{i} = texlabel(outputElements(i));
				end
			end
		end
	end
	
	methods % Getter methods
		function T = get.T(obj)
			T = 0;
			for i = 1:obj.N
				T = T + 1/2 * obj.coordData{i,5} * diff(obj.coordData{i,1}, obj.t)^2;
			end
		end
		function N = get.N(obj)
			N = size(obj.coordData);
			N = N(1);
		end
		function q = get.q(obj)
			for i = 1:obj.N
				q(i) = obj.coordData{i,1};
			end
		end
		function inertia = get.inertia(obj)
			for i = 1:obj.N
				inertia(i) = obj.coordData{i,5};
			end
		end
		function Q_original = get.Q_original(obj)
			for i = 1:obj.N
				Q_original(i) = obj.coordData{i,2};
			end
		end
		function c_f = get.c_f(obj)
			for i = 1:obj.N
				c_f(i) = obj.coordData{i,3};
			end
		end
		function Q = get.Q(obj)
			for i = 1:obj.N
				Q(i) = obj.coordData{i,4};
			end
		end
	end
end















