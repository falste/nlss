function rawEqns = lagrange(q, Q, T, V, varargin)
	if nargin < 5
		doPrint = false;
	else
		doPrint = varargin{5};
	end

	syms t;
	qp = sym('qp', [1, length(q)]);
	qv = sym('qv', [1, length(q)]);
	rawEqns = sym('eq', [length(q), 1]);
	
	% Define lagrange equation
	L = T - V;
	L = subs(L, [diff(q, t) q], [qv qp]);
	
	% Evaluate lagrange equation for all i
	for i = 1:length(q)
		Ldq = diff(L, qv(i));
		Lq  = diff(L, qp(i));
		Ldq = subs(Ldq, [qv qp], [diff(q, t) q]);
		Lq  = subs(Lq,  [qv qp], [diff(q, t) q]);
		rawEqns(i) = diff(Ldq, t) - Lq == Q(i);
		rawEqns(i) = isolate(rawEqns(i), diff(q(i), t, t));
	end
	
	% Print input and resulting movement equations
	if doPrint
		fprintf('\t<strong>Lagrange:</strong>\n\t\tT = %s\n', char(T));
		fprintf('\t\tV = %s\n\n', char(V));

		% Print raw equations
		for i = 1:length(rawEqns)
			fprintf('\t\t%s\n', char(rawEqns(i)));
		end
	end
end














