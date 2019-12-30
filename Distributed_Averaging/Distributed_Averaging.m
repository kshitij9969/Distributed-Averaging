% Author: Kshitij Singh.
% Title: Simulation of Thermostatically controlled loads.
% Description: Use the following code to simulate the behavior of N TCLs in
%              the following order.

%% 1.Clear and close all variables and processes.
clear; clc;close all; % Clear all variables, clear command window, close all other processes.
%% 2.Parameters.
P= 14e+3;                                                 %Power consumed by single TCL.
N=1000;                                                   %Number of TCLs.
Pagg=N*P;                                                 %Total power of N TCLs.
Control_temp=20;                                          %Temperature to be maintained.
Ambient_temp=32;                                          %Ambient temperature.
R=2e-3;                                                   %Thermal Resistance.
C=10e+3;                                                  %Thermal Capacitance.
duty=linspace(0.5,0.5,N);                                 %Dutycycle for switching.
f_final=linspace(0.26,0.36,N);                              %Frequencies of the TCLs.
f_settle=mean(f_final);                                   %Mean frequency of the TCLs.
h=0.1
iter=4;                                                 %Steps taken to reach mean frequency.
t_s=0;                                                    %Starting time of transient part.
t_f=10;                                                   %Ending time of transient part.
t_total=t_s:h:t_f;                                        %Total transient time.
ensemble_vary_s=zeros(N,(numel(t_total)));                %Switching signal of transient part.
time=0;                                                   %Just a variable for resetting time.
time_index=0;                                             %Just a variable of algorithmic significance.
Power=zeros(1,(numel(t_total)));                          %For storing Power
m=1
%% 3. Delay
alpha=0;                                                  
for i=1:N-1                                               %Looping and storing delay.
    b=2*pi*i*m/N;
    alpha=[alpha;b];
end
alpha=alpha/iter;
%% 4. Generating signals.
for n=1:N                                                 %Caculating steps of frequency.
    steps=(abs(f_settle-f_final(n)))/iter;
    start_index=1;
   f_start=f_final(n);
   
 for i=1:iter                                             %Iterating for one TCL.

 t=0:h:1/(f_start);
 if (i)*alpha(n)>=alpha(n)
     delay=iter*alpha(n);
 else
     delay=(i)*alpha(n);
 end
 
 
 ensemble_vary_s(n,(start_index:start_index+numel(t)-1))=heaviside(sin(2*pi*f_start*t+delay)-sin((pi-(2*pi*duty(n)))/2));%Generating signals.
 start_index=start_index+numel(t);

if f_start>f_settle                                        %Conditions for stopping the stepping.
f_start=(f_start-steps);
end
if f_start<f_settle
    f_start=(f_start+steps);
end
if f_start==f_settle
   f_start=f_settle;
end



 end
end
%Power

for i=1:numel(t_total)


Power(i)=sum(P*(ensemble_vary_s(:,i)));


end
plot(t_total,Power)
hold on


t=10:0.1:200;                           %Steady State time.
ensemble_fake_s=zeros(N,numel(t));      %Storing steay state switching.
alpha=0;
f_final=linspace(0.2,0.3,N);
f_settle=mean(f_final);
duty=linspace(0.4,0.5,N);
for i=1:N-1
    b=2*pi*i*m/N;
    alpha=[alpha;b];
end

for n=1:N
ensemble_fake_s(n,:)=heaviside(sin(2*pi*f_settle*t+alpha(n))-sin((pi-(2*pi*duty(n)))/2));%Generating signals.
end

%Storing the steady state power.
Power=zeros(1,numel(t));
P=14e+3;
for i=1:numel(t)
    
    Power(i)=sum((P)*((ensemble_fake_s(:,(i)))));
 
end
title('Heterogenous TCLs Power vs Time');
xlabel('Time(hours)');
ylabel('Power(KW)');
plot(t,Power)

hold on
rms(Power)

