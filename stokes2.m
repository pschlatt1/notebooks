% Stokes second problem
%
% Philipp Schlatter 2021

clear all;close all;
omega=2*pi/6;    % angular frequency of the plate
V=1;             % plate velocity   
nu=10;           % kinematic viscosity

ymax=30;         % domain height
tmax=2*pi/omega; % period
dy=0.25; 
dt=tmax/50;
y=0:dy:ymax;
t=0:dt:5*tmax;
k=sqrt(omega/2/nu);
um = V*exp(-k*y);
delta=4.6/k;

fprintf('boundary layer thickness: %f\n',delta)

%vidfile = VideoWriter('stokes2.mp4','MPEG-4');
%open(vidfile);

for i=1:length(t)

  u = V*exp(-k*y).*cos(k*y-omega*t(i));
  clf;hold on
  plot(u,y*k,'b');
  plot(um,y*k,'g--')
  plot(-um,y*k,'g--')
  plot([-0.1 0.1],[delta*k delta*k],'g--')
  axis([-V V 0 ymax*k]);
  box on;grid on
  xlabel('u/V')
  ylabel('y*k')
  title('Stokes second problem')
  
  % writeVideo(vidfile,getframe(gcf));
  pause(0.1);
      
end

%close(vidfile)
