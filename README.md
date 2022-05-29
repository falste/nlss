# nonlinear-statespace-matlab
A nonlinear state-space (nlss) model for MATLAB similar to MATLABs own linear state-space (ss) model.

You can create a nonlinear state-space model

dx/dt = f(t,x,u)

y = g(t,x,u)

as well as connect them in series and simulate them automatically.

Building a nlss can be done by providing the symbolic functions f and g or by using the NLSSBuilder, which is based on the lagrange formalism.
