%Program: predator_prey
%Predator-prey model for lions and gazelles

%Author: Darci Snowden
%Date: October 22, 2020

clear all;
close all;

%Define initial values
nl0=500;  %number of lions at t=0
ng0=3000; %number of gazelles at t=0

%Define h, tmin, tmax
tmin=0;
tmax=25;     %max time, years
h=0.1;       %time step, years

[t,u]=rk4mat(@predatorprey,tmax,tmin,h,[nl0 ng0]);

%Extract number of lions, 1st column of u array
nl=u(1,:);
%Extract number of gazallels, 2nd column of u array
ng=u(2,:);

plot(t,nl,'r',t,ng,'b','LineWidth',2);
legend('Lions','Gazelles');
set(gca,'FontSize',14.);
xlabel('Time (years)','FontSize',14.)
ylabel('Number of Animals','FontSize',14.)


function dudt = predatorprey(t,u)

    %u(1)=nl;
    %u(2)=ng;

    %Define constants
    bg= 1.1;%birthrate of gazelles, 1/years
    bl= 0.00025; %birthrate of lions, 1/years
    dg= 0.0005; %deathrate of gazelles, 1/years
    dl= 0.7; %deathrate of lions, 1/years
    
    dudt=zeros(2,1);
    %Rate of change of lions
    dudt(1)=bl*u(1)*u(2)-dl*u(1);
    %Rate of change of gazelles
    dudt(2)=bg*u(2)-dg*u(2)*u(1);

end


%Use RK4 method in matrix form to study how the populations change with time
%u(1) is number of lions,nl
%u(2) is the number of gazelles,ng
%In 2D u array, the columns time
function [t,u] = rk4mat(dudt,tmax,tmin,h,u0)

%Number of iterations
N=round((tmax-tmin)/h);

%Initial Values
t(1)=tmin;

%Sets the rows in first column equal to u0
%u=[u0(1);
%   u0(2)]
    
u(:,1)=u0;  


%Loop to solve ODE
for i=1:N
    
    %Find midpoint values
    %u(:,i) is both rows at time "i"
    %k1,k2,k3,k4 will be 2x1 arrays

    k1=h*dudt(t(i),u(:,i));
    k2=h*dudt(t(i)+.5*h,u(:,i)+0.5*k1);   
    k3=h*dudt(t(i)+.5*h,u(:,i)+0.5*k2);
    k4=h*dudt(t(i)+h,u(:,i)+k3);
    
    %Find average of midpoint values and update solution
    %u(:,i+1) sets the next column in my 2 row u array
    u(:,i+1)=u(:,i)+1/6*(k1+2*k2+2*k3+k4);
    t(i+1)=t(i)+h;

end

end