clear all;close all;
% The equation we wish to solve is 
% f''' + (1/2)*f*f'' = 0 
% with f(0) = 0, 
%      f'(0) = 0,
%      f'(inf) = 1. 
% We recast this problem as a system of first-order ODEs: 
%       f  = [f ; f' ; f'' ] = [f(1); f(2); f(3)] so that 
%       f' = [f'; f''; f'''] = [f(2); f(3); -(1/2)*f(1)*f(3)] 
% with f(1)(0) = 0, f(2)(0) = 0, f(2)(inf) = 1.
%
% The result should show that:
% f''(0) = 0.332034
%
% Adapted by Philipp Schlatter 2020
%
% Correct results:
% d0  =    4.908980724033209
% d1  =    1.720787657520373
% d2  =    0.664114672430402
% H12 =    2.59110019542704     

[y,f] = fsc(30);

disp(sprintf('f(0) = %f',f(1,3)))
disp(sprintf('cf   = %f Re_x^(-1/2)',2*f(1,3)))

u=f(:,2);
d99= interp1(u,y,0.99);
d1i = trapz(y,1-u);
d1 = y(end)-f(end,1);
d2i = trapz(y,u.*(1-u));
d2 =2*f(1,3);
H12= d1/d2;
H12i= d1i/d2i;
disp(sprintf('H12 = %f',H12))
disp(sprintf('H12i= %f',H12i))
disp(sprintf('d1/d99 = %f',d1/d99))
disp(sprintf('d2/d99 = %f',d2/d99))
disp(sprintf('d1   = %f xRe_x^(-1/2)',d1))
disp(sprintf('d1i  = %f xRe_x^(-1/2)',d1i))
disp(sprintf('d2   = %f xRe_x^(-1/2)',d2))
disp(sprintf('d2i  = %f xRe_x^(-1/2)',d2i))
disp(sprintf('d99  = %f xRe_x^(-1/2)',d99))


figure; hold on;
plot(f(:,1),y,'k-','Linewidth',2)
plot(f(:,2),y,'r-','Linewidth',2)
plot(f(:,3),y,'b-','Linewidth',2)
box on;grid on
axis([0 2.5 0 10])
xlabel('f,f^\prime,f^{\prime\prime}')
ylabel('y/\Delta')
legend('f','f^\prime','f^{\prime\prime}')
title('Blasius solution')


function [y,f] = fsc(yl)
% Use fsolve to ensure the boundary function g = f'(inf)-1 = 0
opt = optimset('Display','off','TolFun',1E-10);
F = fsolve(@(F) eval_boundary(F,yl),0.3,opt);
% Solve the ODE-IVP with the converged initial condition F
[y,f] = solve_ode(F,yl);
end


function [y,f] = solve_ode(F,yl)
%[y,f] = ode45(@(y,f) [f(2); f(3); -0.5*f(1)*f(3)],[0 yl],[0 0 F]);                
[y,f] = ode45(@(y,f) [f(2); f(3); -0.5*f(1)*f(3)],[0:0.01:yl],[0 0 F]);                
end 


function [g] = eval_boundary(F,yl)
% Get the solution to the ODE with inital condition F
[y,f] = solve_ode(F,yl);
% Evaluate the boundary function f'(inf) - 1
g = f(end,2)-1;
end