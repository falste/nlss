function rawEqns = lagrange(q, Q, T, V)
	doPrint = true;

	syms t;
	qp = sym('qp', [1, length(q)]);
	qv = sym('qv', [1, length(q)]);
	
	L = T - V;
	L = subs(L, [diff(q, t) q], [qv qp]);
	
% 	dq = symfunVector('dq', length(q), t);
	rawEqns = sym('eq', [length(q), 1]);
	
	for i = 1:length(q)
% 		dLdq_dot = subs(L, diff(q(i),t), dq(i));
% 		dLdq_dot = functionalDerivative(dLdq_dot, dq(i));
% 		dLdq_dot = subs(dLdq_dot, dq(i), diff(q(i),t));
% 		
% 		rawEqns(i) = diff(dLdq_dot, t, 1) - functionalDerivative(L, q(i)) == Q(i);
% 		rawEqns(i) = isolate(rawEqns(i), diff(q(i), t, t));

		Ldq = diff(L, qv(i));
		Lq  = diff(L, qp(i));
		Ldq = subs(Ldq, [qv qp], [diff(q, t) q]);
		Lq  = subs(Lq,  [qv qp], [diff(q, t) q]);
		rawEqns(i) = diff(Ldq, t) - Lq == Q(i);
		rawEqns(i) = isolate(rawEqns(i), diff(q(i), t, t));
	end
	
	if doPrint
		fprintf('\t<strong>Lagrange:</strong>\n\t\tT = %s\n', char(T));
		fprintf('\t\tV = %s\n\n', char(V));

		% Print raw equations
		for i = 1:length(rawEqns)
			fprintf('\t\t%s\n', char(rawEqns(i)));
		end
	end
end














