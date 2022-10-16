% Solution to the Orr-Sommerfeld equation in channel flow
%
% Scaling as in lecture notes.
% Philipp Schlatter 2021

clear all
close all

Re = 5000;    % Re=Re_cl=U_cl*h/nu=3/2*Reb
alpha=1.0255; % streamwise wavenumber
N=200;        % resolution (number of GLC points)

% 2nd- and 4th-order differentiation matrices:
[y,D]=chebdif(N,1);
D2 = D^2; D4=D2^2;

% Orr-Sommerfeld operators A,B and generalized eigenvalues:
% A*v = lambda*B*v
I = eye(N);
A = -(D4-2*D2*alpha^2+I*alpha^4)/Re + 2i*alpha*I + 1i*alpha*diag(1-y(1:N).^2)*(D2-I*alpha^2);
B = 1i*(D2-I*alpha^2);

% A has a rank deficit of 4 (4th order derivatives), B of 2 (2nd order)

% Boundary conditions
% move the eigenvalues of the BC to Q
Q=-999i;
% two homogeneous Dirichlet conditions
A(1,:) = Q*I(1,:);
B(1,:) = I(1,:);
A(N,:) = Q*I(N,:);
B(N,:) = I(N,:);

% two homogeneous Neumann conditions
A(2,:) = Q*D(1,:);
B(2,:) = D(1,:);
A(N-1,:) = Q*D(N,:);
B(N-1,:) = D(N,:);

% compute eigenvalues and eigenvectors (A*v = lambda*B*v)
[vv omega] = eig(A,B);

% alternatively, formulate regular eigenvalue problem M*v=lambda*v
%[vv omega] = eig(inv(B)*A);

omega=diag(omega);
c=omega/alpha;

% extract leading eigenvalue (TS wave)
[cmax i]=max(imag(c));
cmax = c(i);
v=vv(:,i);
u=-D*v/(1i*alpha);   % use continuity equation i alpha u + Dv = 0
disp(sprintf('Leading eigenvlue c=%f + i*%f',real(cmax),imag(cmax)))

% scale/rotate eigenvectors to set phase
[w i]=max(abs(u));
w=u(i);
u=u/w;
v=v/w;

% plot results
figure; hold on
plot(real(c),imag(c),'k.','markersize',12)
plot(real(cmax),imag(cmax),'ro','markersize',12)
axis([ 0 1 -1 0.1])
grid on
xlabel('c_r');ylabel('c_i') 
box on
title(sprintf('Spectrum for Re=%5.0f and alpha=%2.3f',Re,alpha))

figure; hold on
plot(y,real(u),'r')
plot(y,imag(u),'r--')
plot(y,abs(u),'k-')
xlabel('y')
ylabel('u_{TS}')
box on
title(sprintf('TS-wave with c=%f + i*%f',real(cmax),imag(cmax)))

figure; hold on
plot(y,real(v),'r')
plot(y,imag(v),'r--')
plot(y,abs(v),'k-')
xlabel('y')
ylabel('v_{TS}')
box on
title(sprintf('TS-wave with c=%f + i*%f',real(cmax),imag(cmax)))
 


  
  
