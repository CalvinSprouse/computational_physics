%Program: predator_prey
%Predator-prey model for lions and gazelles

%Author: Darci Snowden
%Date: October 22, 2020

clear all;
close all;

%Define constants
bg= 1.1;     %birthrate of gazelles, 1/years
bl= 0.00025; %birthrate of lions, 1/years
dg= 0.0005;  %deathrate of gazelles, 1/years
dl= 0.7;     %deathrate of lions, 1/years

%Define initial values
nl0=500;      %number of lions at t=0
ng0=3000;     %number of gazelles at t=0

%Define h, tmin, tmax
tmin=0;
tmax=25;     %max time, years
h=0.1;       %time step, years

%Define ODE's
dnldt=@(nl,ng) bl*nl*ng-dl*nl;   %Rate of change lions
dngdt=@(nl,ng) bg*ng-dg*ng*nl;   %Rate of change of gazzeles

%[t,nl,ng]=odeEuler(dnldt,dngdt,tmax,tmin,h,nl0,ng0);
[t,nl,ng]=rk4(dnldt,dngdt,tmax,tmin,h,nl0,ng0);


plot(t,nl,'r',t,ng,'b','LineWidth',2);
legend('Lions','Gazelles');
set(gca,'FontSize',14.);
xlabel('Time (years)','FontSize',14.)
ylabel('Number of Animals','FontSize',14.)


%Use Euler Method to study how the populations change with time
%y1 is number of lions,nl
%y2 is the number of gazelles,ng
function [t,y1,y2]=odeEuler(ode1,ode2,tmax,tmin,h,y10,y20)

    %Number of iterations
    N=round(tmax-tmin)/h;

    %Initial Values
    y1(1)=y10;
    y2(1)=y20;
    t(1)=tmin;

    %Loop to update values using Euler's method
    for i=1:N
        y1(i+1)=y1(i)+ode1(y1(i),y2(i))*h;
        y2(i+1)=y2(i)+ode2(y1(i),y2(i))*h;
        t(i+1)=t(i)+h;
    end
end

%Use RK2 Method to study how the populations change with time
%y1 is number of lions,nl
%y2 is the number of gazelles,ng
function [t,y1,y2] = rk4(ode1,ode2,tmax,tmin,h,y10,y20)

%Number of iterations
N=round((tmax-tmin)/h);

%Initial Values
t(1)=tmin;
y1(1)=y10;
y2(1)=y20;

%Loop to solve ODE
for i=1:N
    
    %Find midpoint values
    k1=h*ode1(y1(i),y2(i));
    j1=h*ode2(y1(i),y2(i));
    
    k2=h*ode1(y1(i)+.5*k1,y2(i)+.5*j1);
    j2=h*ode2(y1(i)+.5*k1,y2(i)+.5*j1);
    
    k3=h*ode1(y1(i)+.5*k2,y2(i)+.5*j2);
    j3=h*ode2(y1(i)+.5*k2,y2(i)+.5*j2);
    
    k4=h*ode1(y1(i)+k3,y2(i)+j3);
    j4=h*ode2(y1(i)+k3,y2(i)+j3);
    
    %Find average of midpoint values and update solution
    y1(i+1)=y1(i)+1/6*(k1+2*k2+2*k3+k4);
    y2(i+1)=y2(i)+1/6*(j1+2*j2+2*j3+j4);
    t(i+1)=t(i)+h;

end

end