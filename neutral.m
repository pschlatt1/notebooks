% Neutral curve for channel flow
%
% Philipp Schlatter 2021

clear all
close all

N=50;      % resolution

% 2nd- and 4th-order differentiation matrices:
[y,D]=chebdif(N,1);
D2 = D^2; D4=D2^2;
I = eye(N);

[RE,ALPHA] = meshgrid(1000:500:40000, 0.5:0.05:1.2);

for ii=1:size(RE,1)
    for jj=1:size(ALPHA,2)
        Re=RE(ii,jj);alpha=ALPHA(ii,jj);
        A = -(D4-2*D2*alpha^2+I*alpha^4)/Re + 2i*alpha*I + 1i*alpha*diag(1-y(1:N).^2)*(D2-I*alpha^2);
        B = 1i*(D2-I*alpha^2);

        Q=-999i;
        A(1,:) = Q*I(1,:);
        B(1,:) = I(1,:);
        A(N,:) = Q*I(N,:);
        B(N,:) = I(N,:);
        A(2,:) = Q*D(1,:);
        B(2,:) = D(1,:);
        A(N-1,:) = Q*D(N,:);
        B(N-1,:) = D(N,:);

        [vv omega] = eig(A,B);
        omega=diag(omega);

        g(ii,jj)=max(imag(omega));
    end
end

figure;
hold on
pcolor(RE,ALPHA,g);shading interp;
contour(RE,ALPHA,g,[0 0],'k')
plot(5772,1.0255,'ro') % Orszag 1971
xlabel('Re')
ylabel('\alpha')
set(gca,'layer','top');grid on; box on
title(sprintf('Neutral curve channel flow'))

