classdef NLSS < handle
% A nonlinear state-space system with structure
% dx/dt = f(t,x,u)
%     y = g(t,x,u).
%
% Create a nonlinear state-space model using
%     sys = NLSS(t, x, f, [u], [g], [substi])
% or convert a system to nlss using
%     sys = NLSS(linsys)




	properties (SetAccess = private)
		f;
		g;
		t;
		u;
		x;
		
		substitutions;
	end
	
	properties 
		ulabels;
		xlabels;
		ylabels;
	end
	
	properties (Dependent)
		n;
		p;
		m;
	end
	
	methods
		function n = get.n(obj)
			if ~isempty(obj.f)
				n = length(formula(obj.f));
			else
				n = 0;
			end
		end
		function p = get.p(obj)
			if ~isempty(obj.g)
				p = length(formula(obj.g));
			else
				p = 0;
			end
		end
		function m = get.m(obj)
			if ~isempty(obj.u)
				m = length(formula(obj.u));
			else
				m = 0;
			end
		end
		
		function obj = NLSS(varargin)
			if nargin == 1
				sys = varargin{1};
				obj = NLSS.ToNLSS(sys);
				return;
			else
				if nargin == 2 || nargin > 6
					error('Unsupported number of inputs! See ''help nlss'' for usage.');
				end
				t = varargin{1};
				x = varargin{2};
				f = varargin{3};
				
				if nargin > 3, u      = varargin{4}; else u      = sym('u', [0, 1]); end
				if nargin > 4, g      = varargin{5}; else g      = [] ; end
				if nargin > 5, substi = varargin{6}; else substi = []; end
			end
			
			
			if nargin <= 5
				substi = [];
			end
			
			% Setup internal variables for t, u and x
			obj.t = sym('t');
			obj.u = sym('u', [length(formula(u)), 1]);
			obj.x = sym('x', [length(formula(x)), 1]);
			
			% Substitute inputs x and u with the internally generated variables
			% x and u and store these substitutions. This is done to guarantee a
			% consistent internal labelling style.
			obj.substitutions = formula([t == obj.t; u == obj.u; x == obj.x; substi]);
			obj.f = formula(subs(f, lhs(obj.substitutions), rhs(obj.substitutions)));
			obj.g = formula(subs(g, lhs(obj.substitutions), rhs(obj.substitutions)));
			
			% Set labels
			obj.xlabels = cell(obj.n, 1);
			for i = 1:obj.n
				obj.xlabels{i} = ['x_{' num2str(i) '}'];
			end
			obj.ylabels = cell(obj.p, 1);
			for i = 1:obj.p
				obj.ylabels{i} = ['y_{' num2str(i) '}'];
			end
			obj.ulabels = cell(obj.m, 1);
			for i = 1:obj.m
				obj.ulabels{i} = ['u_{' num2str(i) '}'];
			end
		end
		
		function [t_sim, x_sim, y_sim] = sim(obj, x0, T, u)
			% Simulate the system using the initial conditions x0 for the
			% duration T with the inputs u. The results are plotted if
			% nargout == 0.
			
			if nargin == 3
				u = sym(zeros(obj.m, 1));
			end
			if isa(u, 'double')
				u = sym(u);
			end
			
			if length(formula(u))~=obj.m
				error('Wrong number of inputs provided!');
			end
			
			fnum = subs(obj.f, obj.u, u);
			gnum = subs(obj.g, obj.u, u);

			fprintf('\n\tSimulating...');
			% Calculate x using ode45
			odefun = @(ts,xs) double(subs(fnum, [obj.x; obj.t], [xs; ts]));
			[t_sim, x_sim] = ode45(odefun, [0 T], x0);
			
			% Calculate y from x
			y_sim = NaN(length(t_sim), obj.p);
			if obj.p ~= 0
				for i = 1:length(t_sim)
					y_sim(i,:) = subs(gnum.', [obj.x; obj.t], [x_sim(i,:).'; t_sim(i)]);
				end
			end
			
			% Calculate discrete steps of u from input u
			u_sim = NaN(length(t_sim), obj.m);
			if obj.m ~= 0
				for i = 1:length(t_sim)
					u_sim(i,:) = subs(u.', [obj.x; obj.t], [x_sim(i,:).'; t_sim(i)]);
				end
			end
			
			fprintf(' Done.\n');
			
			if nargout == 0
				nplots = 2;
				if obj.p ~= 0
					nplots = nplots + 1;
				end
				if obj.m ~= 0
					nplots = nplots + 1;
				end
				plotindex = 1;
				
				clf;
				if obj.m ~= 0
					subplot(nplots,1,plotindex);
					plot(t_sim, u_sim);
					legend(obj.ulabels);
					title('Input');
					plotindex = plotindex + 1;
				end
				
				if obj.p ~= 0
					subplot(nplots,1,plotindex);
					if ~isempty(obj.g)
						plot(t_sim, y_sim);
						legend(obj.ylabels);
					end
					title('Output');
					plotindex = plotindex + 1;
				end
				
				subplot(nplots,1,plotindex:plotindex+1);
				plot(t_sim, x_sim, '-');
				legend(obj.xlabels);
				xlabel(char(obj.t));
				title('States');
			end
		end
		
		function nlss = mtimes(obj1,obj2)
			% Create one big system out of the serial connection of two systems.
			% The two systems can be of any type that can be converted to ss or
			% nlss.
			
			syses = {obj1, obj2};
			
			% Convert both to NLSS
			for i = 1:length(syses)
				syses{i} = NLSS(syses{i});
			end
			
			if syses{1}.p ~= syses{2}.m
				error('Output dimension of first system must equal input dimension of second system.');
			end
			
			% Combine syses
			n_new = syses{1}.n + syses{2}.n;
			m_new = syses{1}.m;
			%p_new = syses{2}.p;
			
			syms t;
			x_new = sym('x', [n_new 1]);
			u_new = sym('u', [m_new 1]);
			
			gs1 = subs(syses{1}.g, [syses{1}.x; syses{1}.u], [x_new(1:syses{1}.n); u_new]);
			gs2 = subs(syses{2}.g, [syses{2}.x; syses{2}.u], [x_new(syses{1}.n+1:end); gs1]);
			
			fs1 = subs(syses{1}.f, [syses{1}.x; syses{1}.u], [x_new(1:syses{1}.n); u_new]);
			fs2 = subs(syses{2}.f, [syses{2}.x; syses{2}.u], [x_new(syses{1}.n+1:end); gs1]);
			
			f_new = [fs1; fs2];
			g_new = gs2;
			
			substi = [];
			if isa(obj1, 'NLSS')
				substi = [substi; obj1.substitutions];
			end
			if isa(obj2, 'NLSS')
				substi = [substi; obj2.substitutions];
			end
			nlss = NLSS(t, x_new, f_new, u_new, g_new, substi);
			if isa(obj1, 'NLSS')
				nlss.xlabels(1:obj1.n) = obj1.xlabels;
				nlss.ulabels(1:obj1.m) = obj1.ulabels;
			end
			if isa(obj2, 'NLSS')
				nlss.xlabels(n_new-obj2.n + 1:end) = obj2.xlabels;
				nlss.ylabels(1:obj2.p) = obj2.ylabels;
			end
		end
	end
	
	methods (Static, Access = private)
		function nlss = ToNLSS(sys)
			if ~(isa(sys, 'tf') || isa(sys, 'ss') || isa(sys, 'NLSS'))
				error(['System is a ' class(sys) '. Use tf, ss or NLSS!']);
			end

			% Convert to ss if it's a tf
			if isa(sys, 'tf')
				sys = ss(sys);
			end

			% Convert to NLSS if it's a ss
			if isa(sys, 'ss')
				sys = NLSS.ss2nlss(sys);
			end
			
			nlss = sys;
		end
		
		function nlss = ss2nlss(ss_sys)
			n = size(ss_sys.A, 1);
			m = size(ss_sys.B, 2);
			%p = size(ss_sys.C, 1);
			
			syms t;
			x = sym('x', [n, 1]);
			u = sym('u', [m, 1]);
			f = ss_sys.A*x + ss_sys.B*u;
			g = ss_sys.C*x + ss_sys.D*u;
			
			nlss = NLSS(t, x, f, u, g);
		end
	end
end




























