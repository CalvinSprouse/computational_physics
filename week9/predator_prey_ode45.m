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

[t,u]=ode45(@predatorprey,[tmin tmax],[nl0; ng0]);

%Extract number of lions, 1st column of u array
nl=u(:,1);
%Extract number of gazallels, 2nd column of u array
ng=u(:,2);

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


