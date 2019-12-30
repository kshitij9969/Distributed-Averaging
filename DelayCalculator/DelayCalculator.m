% Author: Kshitij Singh

% Description: 
% This code generate a data set for delay calculator. Where N
% varies from 2 to 500, P varies from 100% to 40% and alpha varies from 0
% to 2*pi/N in 100 steps. All the parameters are obtained by observing a
% real TCL.

% Final dataset is stored in final_N_power


%%
clear;
clc;
close all;
%%
% Choosing power between 60% and 140% of 1.66kW
% Choosing dutycycle between 0.422 and 0.482
% The mean frequency is chosen as 0.5
clear all
N=500;
t = 0 : 0.01 : 2*pi;
a=(0.6*1.66);b=(1.4*1.66);
f =  a + (b-a).*rand(N,1);
final_N_power = zeros(N,101);
power = a + (b-a).*rand(N,1);
a = 0.422;
b = 0.482;
dutycycle = a + (b-a).*rand(N,1);
s_temp = zeros(1, numel(t));
for m = 2 : N
    h = 0.01*(2*pi/m);
    alpha = 0 : h : 100*h;
    for i = 1:length(alpha)
        s = heaviside(sin(2*pi*0.31*t));
        delay = (2*pi*0.31*t + alpha(i));
        for n = 1 : m-1
            s0 = sin((pi - dutycycle(n)*2)/2); %1/f = 1/0.5
            s_temp = (power(n))*heaviside(sin(delay)-s0);
            s = s + s_temp;
            delay = delay + alpha(i);
            clear s_temp
        end
           final_N_power(m,i) = rms(s(1,:));
    end
    clear s
end










