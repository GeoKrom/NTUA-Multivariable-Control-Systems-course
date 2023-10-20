%% Real System Simulation (Sine)
% Created for the Multivariable Systems Lab
% Athens 02/2023
% Authors: 
%   G. Kassavetakis AM 02121203
%   G. Krommydas    AM 02121208

clc
clear
close all

%% Parameters
xvmax = 2.79e-3;
s = tf('s');
% Electro-hydraulic Parameters 
% System (linear)with input Volt and output x_v normalised 
% System Equation (Laplace):
% x_v/u =(Kvp/xvmax)/(tau*Av*s^2 + Av*s + K3*Kvp)
Kvp = 1.06e-5;
tau = 0.0014;
Av = 1.964e-4;
K3 = 3579.13;

% Valve Parameters
Ps = 1.9e4;
Kv = 1.64e-2; % Kv = C_d*w*xvmax*sqrt(P_s/p)
K1 = 2*Kv/xvmax;

% Actuator Parameters
% Common Parameters for Eq1 and Eq2
A = 8.212e-3;
% Equation 1: Q_L = Ka*(dP_L/dt) + Cl*P_L + A*(dx_p/dt)
Ka = 7.56e-10;
Cl = 1.3e-8;
% Equation 2: A*P_L - Fc = mt*(d^2x_p/dt^2) + c*(dx_p/dt) + kx_p
%Fc = 0.1; %Not USED
mt = 500;
c = 4.33;
k = 698.8;

%% Linear Model parameters

Kp = 0; %Using 0 flow back due to Pressure diff
p = 0.0893; %Density from gamma/g=0.876/9.81=p
b3 = Ka*mt/A;
b2 = Ka*c/A + mt*(Cl + Kp)/A;
b1 = A*p + Ka*k/A + (Cl + Kp)*c/A;
% OR b1 = A + Ka*k/A + (Cl + Kp)*c/A;
b0 = k*(Kp+Cl)/A;

%% System Creation

tf_valve = tf(Kvp/xvmax,[tau*Av,Av,K3*Kvp]);
tf_actuator = tf(Kv,[b3,b2,b1,b0]);
tf_system = tf_valve*tf_actuator;
% tf_system = tf_system*10^5/10^5/(1/1.266e-06);
%% Input Creation

% Time vector
t = 0:0.01:5000;

% For Persistance of Excitation we need n/2=4 frequencies in order to obtain
% all the parameters
% f_train = [0.01, 0.02, 0.05]; %Frequency in Hz
% A_train = [0.1, 0.5, 0.25]; %Amp in Volt
f_train = [0.01, 0.02, 0.05, 0.08, 0.07]; %Frequency in Hz
A_train = [0.1, 0.5, 0.25, 0.1, 0.05]; %Amp in Volt

% Testing will be done using a sum of a step input and a sine wave (No need
% for the Persistance of Excitation condition here)
f_test = 0.015; %Frequency in Hz
A_test = 0.5; %Amp in Volt

% Input signals
u_train = A_train*sin(2*pi*f_train'*t);
u_test = ones(size(t)) + A_test*sin(2*pi*f_test*t);

%% System Response for sum of Sines Input

%Simulation for the training signal
y_train = lsim(tf_system,u_train,t);

%Simulation for the testing signal
y_test = lsim(tf_system,u_test,t);

% Plot of the 2 system responses
figure(1)
clf
subplot(2,2,1)
plot(t,u_train)
grid on
title('Input for the model training u_{train}')
ylabel('Volt')
xlabel('Time[s]')
subplot(2,2,2)
plot(t,y_train)
grid on
title('System Response for the Input u_{train}')
ylabel('Amplitude')
xlabel('Time[s]')
subplot(2,2,3)
plot(t,u_test)
grid on
title('Input for the model testing u_{test}')
ylabel('Volt')
xlabel('Time[s]')
subplot(2,2,4)
plot(t,y_test)
grid on
title('System Response for the Input u_{test}')
ylabel('Amplitude')
xlabel('Time[s]')

%Save to Files
names = {'Time','Input','Position'};
output1 = table(t',u_train',y_train,'VariableNames',names);
output2 = table(t',u_test',y_test,'VariableNames',names);
writetable(output1,'SineResultTrain.csv')
writetable(output2,'SineResultTest.csv')
