function FanModel()

a0 = 0.18/10; %originally in mm, converted to cm. Radius of cavity around needle.
R = 0.5; %cm Radius of tumor.
M = 10; %orginally in mmHg, can be converted to dyn/cm2 by *1333.22387
nu = 0.35; %poisson ratio
H = 3*10^(-5); %cm2/mmHg/s
alpha = 0;
lambda = 57; %mmHg
G = lambda*((1-2*nu)/(2*nu));
Pt = 25; %mmHg

r = 0.18:0.01:3;
rspan = [r(1),r(end)]; 

init = [0];
[r,y] = ode45(@diffeqsol,rspan,init);


end

function solution = diffeqsol(r,y,a0,R,M,nu,H,alpha,lambda,G,Pt)

% for r = a0 %y direction BC zero flux (dcdx=0)
%     solution(r) = 2*D_S.*(C(i+1)-C(i))./(dx.^2)- k_E.*C(i);
% end
if Pinf<Pt
    Lp = 0;
else
    Lp = Lp0*(exp(y(3)-Pt/(beta*(2*G+lambda)))-1);
end

solution1 = y(1);
solution2 = -1/(2*G+lambda)*v/k-2*(1/r*y(2)-(y(1)/r^2));
solution3 = -(1/k)*Lp*(y(3)-Pt)*((a0^2)/(r^2));
solution = [solution1;solution2;solution3];

end