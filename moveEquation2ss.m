function [x, f, g, substitutions] = moveEquation2ss(q, moveEqns, outEqns)
% Creates a system with the structure
% dx/dt = f(x,u)
%     y = g(x,u)

	doShowSubstitutions = false;
	doShowSystem = false;
	
	if nargin < 3
		doOutput = false;
	else
		doOutput = true;
	end

	syms t;
	% Create state space model
	
	x = symfunVector('x', [length(q)*2, 1], t);
	NaNx = symfunVector('NaNx', [length(q)*2, 1], t);
	%x = sym('x', [length(q), 1]);
	substitutions = sym('subs', [length(q)*3, 1]);
	f = sym('rhsEq', [length(q)*2, 1]);
	% Replace q1 with x1 and x2, q2 with x3 and x4...
	
	% Setup correct substitutions and define equations
	for i = 1:length(moveEqns)
		f(i) = x(i+length(q));
		f(i+length(q)) = rhs(moveEqns(i));
		
		substitutions(i) = q(i) == x(i);
		substitutions(i+length(q)) = diff(q(i), t) == x(i+length(q));
		substitutions(i+2*length(q)) = diff(q(i), t, t) == f(i+length(q));
	end
	substitutions = lhs(substitutions) == subs(rhs(substitutions), lhs(substitutions), rhs(substitutions));
	
	f = subs(f, lhs(substitutions), rhs(substitutions));
	if doOutput
		g = subs(outEqns, lhs(substitutions), rhs(substitutions));
	else
		g = [];
	end
	
	if doShowSubstitutions
		fprintf('\n\t<strong>Substitutions:</strong>\n');
		for i = 1:length(substitutions)
			fprintf('\t\t%s\n', char(substitutions(i)));
		end

		%disp(table(...
		%	arrayfun(@char, lhs(substitutions), 'uniform', 0),...
		%	arrayfun(@char, rhs(substitutions), 'uniform', 0),...
		%	'VariableNames',{'old', 'new'}));
	end
	
	if doShowSystem
		fprintf('\n\t<strong>System:</strong>\n');
		for i = 1:length(f)
			fprintf('\t\t%s\n', char(diff(x(i),t) == f(i)));
		end
		
		if doOutput
			fprintf('\n');
			gElements = formula(g);

			for i = 1:length(gElements)
				fprintf('\t\t%s\n', ['y' num2str(i) '(t) == ' char(gElements(i))]);
			end

			%disp(table(...
			%	arrayfun(@char, diff(x,t), 'uniform', 0),...
			%	arrayfun(@char, f, 'uniform', 0),...
			%	'VariableNames',{'xdot', 'f'}));
		end
	end
end























