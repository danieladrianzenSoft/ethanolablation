clear
clc
%% Settings
Nz     = 100;
time = 100;
total_num = Nz + 4;
v0 = zeros(total_num,1);
tspan = (0:0.01:time);
%%Parameter
u_int  = 0.10;
H      = 0.20;
deltaZ = H/(Nz+4);
cFeed  = 10;
%% Torsten
JacFun = J(u_int, deltaZ, Nz, cFeed);
options = odeset( 'Jacobian' , JacFun, 'RelTol',1e-8, 'AbsTol',1e-8);
[t,V] = ode15s(@(t,x)odeTest(t,x,u_int, deltaZ, Nz, cFeed), tspan, v0, options);

figure()
plot(t,V)

function dfdx = J(u_int, deltaZ, Nz, cFeed)
  syms t
  x = sym('x' , [1 Nz+4]);
  f1 = (-u_int .*  (x(1)^2      - cFeed) ./ deltaZ );
  f2 = (-u_int .*  (x(2:Nz+4)   - x(1:Nz+3)) ./ deltaZ );
  f = [f1, f2];
  J = jacobian(f,x);
  dfdx = matlabFunction(J,'Vars',{t, [x.']});
end
function dfdt = odeTest(t, vi, u_int, deltaZ, Nz, cFeed)       
  dfdt = zeros(length(vi),1);   
  %% Variables
  u_int  = 0.10;
  H      = 0.20;
  Nz     = 100;
  deltaZ = H/(Nz+4);
  cFeed  = 10;
  x   = vi(1:Nz+4);
  %% DGL
  dfdt(1)              = (-u_int .*  (x(1)^2   - cFeed) ./ deltaZ );
  dfdt(2: Nz+4)        = (-u_int .*  (x(2:Nz+4)   - x(1:Nz+3)) ./ deltaZ );
end