clear all
close all

% Plots the stability region of RK1-RK4
%     RK4 stability criterion:
%     along imaginary axis: sqrt(8.)=2.828
%     along real axis:      2.785
%
%     RK3 stability criterion:
%     along imaginary axis: sqrt(3.)=1.73
%     along real axis:      2.52
%
%   G(z) = sum_k=0^N z^k/k!
%   |G(z)|=1 defines the neutral curve
%
%   G(z) = 1 + z + z^2/2
%   roots ( z^2/2 + z +1-phi = 0)


figure 
hold on

% simple way using contour plot

[X,Y]=meshgrid(-5:.2:5, -5:.2:5);

for nn=1:4

  G=0*X;
  for j=0:nn
    rksum = 1/factorial(j);
    G=G+rksum*(X+i*Y).^(j);
  end

  contour(X,Y,abs(G),[1 1],'g')
  
end

  
% proper way yielding explicit curve
N=100;
phi = exp(i*linspace(0,2*pi,N+1));
phi = phi(1:end-1);

for nn=1:4

  % compute the factors [1/fact(nn) 1/fact(nn-1) ... 1]
  % of the characteristic polynomial, for all factors that include z.
  rksum = 1;
  for j=2:nn
    rksum = [1/factorial(j) rksum];
  end
  
  % and the solutions for each angle of the unit circle, including the
  % final 1.
  p0 = [];
  for j=1:length(phi)
    a = roots([rksum 1-phi(j)]);
    p0 = [p0 a'];
  end
  
  % now reorder all solutions to a continuous line
  pp=zeros(1,length(p0));
  w2=p0;
  w1=sort(w2);
  pp(1) = w1(1);
  % trick: minimise distance from one to the next point
  for j=2:length(p0)
    w2=sort(w2-pp(j-1))+pp(j-1);
    pp(j)=w2(1);
    w2=w2(2:end);
  end
  % close the line
  pp=[pp pp(1)];
  
  if (nn==1)
    plot(pp,'m-')
  elseif (nn==2)
    plot(pp,'b-')
  elseif (nn==3)
    plot(pp,'r-')
  elseif (nn==4)
    plot(pp,'k-')
  end      
    
end

axis equal
axis([-3 0.5 -3 3])


  

figure
[X,Y]=meshgrid(-3:.01:1, -3:.01:3);
Z=X+1i*Y;
G=Z*0;
% simple way by just putting the characteristic polynomial
% (here for AB3)
for i=1:size(Z,1)
    for j=1:size(Z,2)
        a=roots([1;-1-23/12*Z(i,j);16/12*Z(i,j);-5/12*Z(i,j)]);
        G(i,j) = max(abs(a));
    end
end

contour(X,Y,abs(G),[1 1],'g')
axis equal
axis([-3 0.5 -3 3])



