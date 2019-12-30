% Author: Kshitij Singh.
% Title: Simulation of Thermostatically controlled loads.
% Description: Use the following code to simulate the behavior of N TCLs in
%              the following order.
%                      1. Clear and close all variable and processes.
%                      2. Parameters.
%                      3. Kuramoto's Temperature Equation for single TCL.
%                      4. Standard Kuramoto's Equation of N coupled oscillators.
%                      5. Ensemble of TCL using Heaviside function and standard continuous-time distributed average 'concensus' type protocol.
%                      6. Switiching function.
%                      7. Cumulative error at steady-state.
%                      8. Temperature of N TCLs.
%                      9. Phase plot on circle.
%                      10.Fourier transform comparision of Kuramoto Equation and Ensemble TCLs.
%                      11.Aggregate power consumption.
%                      12.Power v/s delay.
%                      13.Delay computer.



%% 1.Clear and close all variables and processes.
clear; clc;close all; % Clear all variables, clear command window, close all other processes.
%% 2. Parameters.
W=linspace(-26*0.9,-26*1.1);                              %Natural Frequency
k=0.1;                                                    %Coupling constant
P= 14e+3;                                                 %Power consumed by single TCL.
N=100;                                                    %Number of TCLs.
Pagg=N*P;                                                 %Total power of N TCLs.
deadband=1;                                               %Deadband for temperature.
Control_temp=20;                                          %Temperature to be maintained.
Ambient_temp=32;                                          %Ambient temperature.
R=2e-3;                                                   %Thermal Resistance.
C=10e+3;                                                  %Thermal Capacitance.
duty=0.43;                                                %Dutycycle for switching.
alpha_delay=[0;pi;-pi/2;-pi/2];                           %Delay alpha to be given to each TCL.
x0=rand(N,1);                                             %Initial Conditions.
multiplier=1.03;                                          %Multiplier for adjusting fourier transform
h=0.01;                                                   %Time step.
t_final=300;                                              %Time for simulation.
tspan= 1e-2:h:t_final;                                    %Time array for solving.


%% 3.Hybrid Equation for single TCL.

Kuramoto_tempdot = @(t,y,s) (-1/(R*C))*(y-Ambient_temp + s*P*R);  %Kuramoto temperature Equation.
Temperature_kuramoto = zeros(1,numel(tspan));                     %For storing temperature.
Temperature_kuramoto(1) = 20.1;                                   %Initial Temperature.
s_kuramoto=0;                                                     %Switching
store_switching=zeros(1,numel(tspan));
store_switching(1)=s_kuramoto;
for i = 2:numel(tspan) 
  
    k1 = h*Kuramoto_tempdot(tspan(i-1),Temperature_kuramoto(i-1),s_kuramoto);
    k2 = h*Kuramoto_tempdot(tspan(i-1)+h/2, Temperature_kuramoto(i-1)+k1/2, s_kuramoto);
    k3 = h*Kuramoto_tempdot(tspan(i-1)+h/2, Temperature_kuramoto(i-1)+k2/2, s_kuramoto);
    k4 = h*Kuramoto_tempdot(tspan(i-1)+h, Temperature_kuramoto(i-1)+k3, s_kuramoto);
    Temperature_kuramoto(i) = Temperature_kuramoto(i-1) + (k1+2*k2+2*k3+k4)/6;
    
  if Temperature_kuramoto(i)<(Control_temp-deadband*0.5)           %Switching based on temperature.
      s_kuramoto=0;
     store_switching(i)=s_kuramoto;
      
  elseif Temperature_kuramoto(i)>(Control_temp+deadband*0.5)   
      s_kuramoto=1;  
       store_switching(i)=s_kuramoto;
  end
  
end

 plotyy(tspan,Temperature_kuramoto,Temperature_kuramoto,store_switching);
% plot(store_switching,Temperature_kuramoto);
xlabel('Time in hours');
ylabel('Temperature in degree C');
title('Kuramoto Temperature v/s time');


%% 4. Standard Kuramoto's Equation of N coupled oscillators.(Synchronization)

Kuramoto_phidot  =@(t,x)[(W+k*(sin(x(2)-x(1))+sin(x(3)-x(1))+sin(x(4)-x(1)))); %Kuramoto's Equation
                (W+k*(sin(x(1)-x(2))+sin(x(3)-x(2))+sin(x(4)-x(2)))); 
                (W+k*(sin(x(2)-x(3))+sin(x(1)-x(3))+sin(x(4)-x(3)))); 
                (W+k*(sin(x(2)-x(4))+sin(x(3)-x(4))+sin(x(1)-x(4))))];
[t,phi_kuramoto]=ode45(Kuramoto_phidot, tspan, x0);                            %Solving using ODE45.
plot(t,phi_kuramoto);
title('Phi v/s time synchronization');
xlabel('Time(h)');
ylabel('Phi');
grid;
%% Ensemble of TCL using Heaviside function and standard continuous-time distributed average 'concensus' type protocol.
% 5. Phidot Equation.

k1=zeros(N,1);                                                                 %Parameters for Runge-Kutta Solution.
k2=zeros(N,1);
k3=zeros(N,1);
k4=zeros(N,1);
phi_heavi=zeros(N,numel(tspan));
phi_heavi(:,1)=asin(x0);                                                       %Substituting the initial conditions.
for i = 2:numel(tspan)                                                         %Runge-Kutta method to solve equation.
    
    k1 = h*Phidot_Heaviside(phi_heavi(:,i-1));
    k2 = h*Phidot_Heaviside(phi_heavi(:,i-1)+k1/2);
    k3 = h*Phidot_Heaviside(phi_heavi(:,i-1)+k2/2);
    k4 = h*Phidot_Heaviside(phi_heavi(:,i-1)+k3);
    phi_heavi(:,i) = phi_heavi(:,i-1) + (k1+2*k2+2*k3+k4)/6;
    
end
plot(tspan,phi_heavi);                                                          %Plot of phi v/s time.
xlabel('Time in hours');
ylabel('Phi');
title('Phi v/s Time De-synchronization');
grid;
%% 6. Switiching function.
ensemble_s = zeros(N,numel(tspan));                                                      %initializing s for storing values.
for a=1:N
        ensemble_s(a,:)= heaviside(sin(phi_heavi(a,:))-sin((pi-(2*pi*duty))/2));
end
plot(tspan,ensemble_s);                                                                  %Plotting switching function v/s time.
xlabel('Switiching function');
ylabel('Time in hours');
title('Switching function with Dutycycle 43%');
%% 7. Cumulative error at steady-state.
Phi_err=zeros(N,numel(tspan));                                                           %Initializing err for storing values.

for i=1:numel(tspan)
Phi_err(:,i)=sum(ones(N,1)*phi_heavi(:,i)'-phi_heavi(:,i)*ones(N,1)',2);
end

plot(tspan,Phi_err);
xlabel('Time in hours');
ylabel('Cumulative error(phi)');
title('Error v/s Time');
grid;

%% 8. Temperature of N TCLs.

Temperature = zeros(N,numel(tspan));
for b=1:N
for i=2:numel(tspan)
    a1 = TemperaturEquationEnsembleTCLs(tspan(i-1),Temperature(b,i-1), ensemble_s(b,i-1));
    a2 = TemperaturEquationEnsembleTCLs(tspan(i-1)+h/2, Temperature(b,i-1)+a1*h/2, ensemble_s(b,i-1));
    a3 = TemperaturEquationEnsembleTCLs(tspan(i-1)+h/2, Temperature(b,i-1)+a2*h/2, ensemble_s(b,i-1));
    a4 = TemperaturEquationEnsembleTCLs(tspan(i-1)+h, Temperature(b,i-1)+a3*h, ensemble_s(b,i-1));
    Temperature(b,i) = Temperature(b,i-1) + (a1+2*a2+2*a3+a4)*h/6;
end
end
plot(tspan,Temperature(1,:));
%% 9. Phase plot on circle.
%Circular plot showing phase separation.
phi_x=cos(phi_heavi);                                                                   %Cosine of phi.
phi_y=sin(phi_heavi);                                                                   %Sine of phi.

sp=linspace(0,2*pi,100);                                                                %Circle for reference.
cx=sin(sp);
cy=cos(sp);

figure
plot(phi_x(:,end),phi_y(:,end),'o',cx,cy);
axis([-1 1 -1 1]);
xlabel('cos(phase)');
    ylabel('sin(phase)');
title('Angular Separation');
axis square;
grid;

%% 10. Fourier transform comparision of Kuramoto Equation and Ensemble TCLs.
%Calculates the fundamental frequency of the signals.

calcFREQ(ensemble_s(1,:),h,tspan);

%% 11. Aggregate power consumption.
p=zeros(N,numel(tspan));
for n=2:numel(tspan)
p(:,n-1)=P*sum(ensemble_s(:,n-1));
end
plot(tspan,p);
xlabel('Time in hours');
ylabel('Power consumed');
title('Power v/s time');
grid;

%% 12. Varying delay and calculating aggregate power.

%Note: Only two TCLs are considered.
t=0:0.01:2*pi;                                 %Time.
alpha=0:0.05:2*pi;                             %Delay.
s_dp=zeros(2,numel(t));                        %For storing switching functions.
m=zeros(numel(alpha),numel(t));                %Summing the switching functions for rms calculation.                         
y=zeros(1,numel(alpha));                       %Storing Power.

for i=1:numel(alpha)
    s_dp(1,:)=heaviside(sin(t));
    s_dp(2,:)=heaviside(sin(t+alpha(i)));
    m(i,:)=(s_dp(1,:)+s_dp(2,:));
    y(1,i)=(14e+3)*rms(m(i,:));
end

plot(alpha,y);
title('Delay Calculator');
xlabel('Delay');
ylabel('Power(rms)');
grid;
%% 13. Delay computer.
%Input the value of alpha get the value of Power.
alpha12 = interp1(y,alpha,3.14,'spline');


%% 14.Test
Kuramoto_phidot  =@(t,x)[((k*(sum(x)-N*x(1)))); %Kuramoto's Equation
                ((k*(sum(x)-N*x(2)))); 
                 ((k*(sum(x)-N*x(3)))); 
                  ((k*(sum(x)-N*x(4)))); 
                   ((k*(sum(x)-N*x(5)))); 
                    ((k*(sum(x)-N*x(6)))); 
                     ((k*(sum(x)-N*x(7)))); 
                      ((k*(sum(x)-N*x(8))))
               ];
[t,phi_kuramoto]=ode45(Kuramoto_phidot, tspan, x0);                            %Solving using ODE45.
plot(t,phi_kuramoto);
title('Phi v/s time synchronization');
xlabel('Time(h)');
ylabel('Phi');
grid;
