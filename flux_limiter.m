clear all
close all

% Test code for Burger equation:
Nx=501;
Nt=1001;
Xend=10;
Tend=5;
x=linspace(0,Xend,Nx);
dx=x(2)-x(1);
t=linspace(0,Tend,Nt);
dt=t(2)-t(1);

% initialize velocity vectors
u_fw=zeros(numel(t),numel(x));
u_lw=zeros(numel(t),numel(x));
u_fl=zeros(numel(t),numel(x));
u_lf=zeros(numel(t),numel(x));
u_mc=zeros(numel(t),numel(x));



% initial condition:
u0=zeros(1,numel(x));
u0(x>=1&x<2)=1;
u0(x>8&x<9)=-1;    
u0(x>=2&x<3)=0;
u0(x>=3&x<4)=2;

u_fw(1,:)=u0;
u_lw(1,:)=u0;
u_fl(1,:)=u0;
u_an(1,:)=u0;
u_lf(1,:)=u0;
u_mc(1,:)=u0;

up = u_mc(1,:);
f_p = zeros(1,length(up));
f_m = zeros(1,length(up));
phi = zeros(1,length(up));
r=phi;

figure
for i=1:numel(t)-1
    
    jj=3:numel(x)-2;
    
    % upwind
    % f_fw_1 = u_fw(i,jj-1).^2/2;
    % f_fw_2 = u_fw(i,jj).^2/2;
    % central
    % f_fw_1 = 0.5*(u_lf(i,jj-1).^2/2+u_lf(i,jj).^2/2);
    % f_fw_2 = 0.5*(u_lf(i,jj).^2/2+u_lf(i,jj+1).^2/2);
    % flux splitting
    f_p(jj) = 0.5*(u_fw(i,jj)+abs(u_fw(i,jj))).*u_fw(i,jj)*0.5;
    f_m(jj) = 0.5*(u_fw(i,jj)-abs(u_fw(i,jj))).*u_fw(i,jj)*0.5;
    f_fw_1  = f_p(jj-1)+f_m(jj);
    f_fw_2  = f_p(jj)+f_m(jj+1);
    
    u_fw(i+1,jj) = u_fw(i,jj) - dt/dx*(f_fw_2-f_fw_1);
    
    % Lax-Friedrichs
    f_lf_1 = 0.5*(u_lf(i,jj-1).^2/2+u_lf(i,jj).^2/2) - dx/2/dt*(u_lf(i,jj)-u_lf(i,jj-1));
    f_lf_2 = 0.5*(u_lf(i,jj).^2/2+u_lf(i,jj+1).^2/2) - dx/2/dt*(u_lf(i,jj+1)-u_lf(i,jj));
    u_lf(i+1,jj) = u_lf(i,jj) - dt/dx*(f_lf_2-f_lf_1);
    
    % 2nd order Lax-Wendroff (Richtmyer)
    up(jj) = 0.5*(u_lw(i,jj+1)+u_lw(i,jj))-dt/2/dx*(u_lw(i,jj+1).^2/2-u_lw(i,jj).^2/2);
    f_lw_1 = up(jj-1).^2/2;
    f_lw_2 = up(jj).^2/2;
    u_lw(i+1,jj) = u_lw(i,jj) - dt/dx*(f_lw_2-f_lw_1);
    
    % 2nd order MacCormack
    up(jj) = u_mc(i,jj) - dt/dx*(u_mc(i,jj).^2/2-u_mc(i,jj-1).^2/2);
    fp = up.^2/2;
    f_mc_1 = 0.5*(fp(jj-1) + u_mc(i,jj).^2/2);
    f_mc_2 = 0.5*(fp(jj) + u_mc(i,jj+1).^2/2);
    u_mc(i+1,jj) = u_mc(i,jj) - dt/dx*(f_mc_2-f_mc_1);
    
    % Flux limiters
    r(jj)=(u_fl(i,jj)-u_fl(i,jj-1))./(u_fl(i,jj+1)-u_fl(i,jj));
    r(isnan(r))=1;
    r=min(r,100000);
    r=max(r,-100000);
    phi = (r+abs(r))./(1+abs(r));    % van Leer
    
    % Lax-Wendroff
    up(jj) = 0.5*(u_fl(i,jj+1)+u_fl(i,jj))-dt/2/dx*(u_fl(i,jj+1).^2/2-u_fl(i,jj).^2/2);
    f_high_1 = up(jj-1).^2/2;
    f_high_2 = up(jj).^2/2;
    % Upwind 
%    f_low_1 = u_fl(i,jj-1).^2/2;
%    f_low_2 = u_fl(i,jj).^2/2;
    % Upwind with flux splitting
    f_p(jj) = 0.5*(u_fl(i,jj)+abs(u_fl(i,jj))).*u_fl(i,jj)*0.5;
    f_m(jj) = 0.5*(u_fl(i,jj)-abs(u_fl(i,jj))).*u_fl(i,jj)*0.5;
    f_low_1  = f_p(jj-1)+f_m(jj);
    f_low_2  = f_p(jj)+f_m(jj+1);
    % Lax Friedrichs
%    f_low_1 = 0.5*(u_fl(i,jj-1).^2/2+u_fl(i,jj).^2/2) - dx/2/dt*(u_fl(i,jj)-u_fl(i,jj-1));
%    f_low_2 = 0.5*(u_fl(i,jj).^2/2+u_fl(i,jj+1).^2/2) - dx/2/dt*(u_fl(i,jj+1)-u_fl(i,jj));
    f_fl_1 = f_low_1-phi(jj-1).*(f_low_1-f_high_1);
    f_fl_2 = f_low_2-phi(jj).*(f_low_2-f_high_2);
    u_fl(i+1,jj) = u_fl(i,jj) - dt/dx*(f_fl_2-f_fl_1);
    
    % sum(u_fl(i+1,:))
    
    
    % Plot
    if floor(i/50)==i/50
        cla;
        plot(x,u_lf(i,:),'b.-')
        hold on
        plot(x,u_lw(i,:),'g.-')
        plot(x,u_fw(i,:),'c.-')
        %plot(x,u_mc(i,:),'k.-')
        plot(x,u_fl(i,:),'m.-')
        %       plot(x(jj),phi,'y')
        
        tt=t(i+1);
        %plot([0,1,1+tt,2+tt/2,2+tt/2,4],[0,0,1,1,0,0],'r')
        
        legend('Lax-Friedrichs','Lax-Wendroff','upwind','FluxLimiters')
        title(['t=' num2str(t(i))])
        drawnow;
        %pause(0.01)
        axis([min(x) max(x) -1.1 2.1])
        xlabel('x')
        ylabel('u')
    end
end

figure
pcolor(x,t,u_fl)
title('Flux limiter')
shading flat
xlabel('x')
ylabel('t')
caxis([-2 2])

figure
pcolor(x,t,u_lw)
title('Lax-Wendroff')
shading flat
xlabel('x')
ylabel('t')
caxis([-2 2])

figure
pcolor(x,t,u_lf)
title('Lax-Friedrichs')
shading flat
xlabel('x')
ylabel('t')
caxis([-2 2])

figure
pcolor(x,t,u_fw)
title('Upwind flux splitting')
shading flat
xlabel('x')
ylabel('t')
caxis([-2 2])

figure
waterfall(x,t(1:10:500),u_fl(1:10:500,:))
xlabel('x')
ylabel('t')
