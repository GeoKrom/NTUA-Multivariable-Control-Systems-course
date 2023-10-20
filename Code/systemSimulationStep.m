%% Real System Simulation (Step)
% Created for the Multivariable Systems Lab
% Athens 02/2023
% Authors: 
%   G. Kassavetakis  AM 02121203
%   G. Krommydas     AM 02121208

clc
clear
close all

%% Parameters
xvmax = 2.79e-3;
s = tf('s');
% Electro-hydraulic Parameters 
% System (linear)with input Volt and output x_v normalised 
% System Equation (Laplace):
% x_v/u = (Kvp/xvmax)/(tau*Av*s^2 + Av*s + K3*Kvp)
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
b0 = k*(Kp + Cl)/A;

%% System Creation

tf_valve = tf(Kvp/xvmax, [tau*Av, Av, K3*Kvp]);
tf_actuator = tf(Kv, [b3, b2, b1, b0]);
tf_system = tf_valve*tf_actuator;

%% System step Response

% For training using the position measurement
[y_train, t_train] = step(tf_system);
u_train = ones(size(t_train));
stepplot(tf_system)
% For testing using the position measurement
[y_test, t_test] = step(tf_system*0.8);
u_test = 0.8*ones(size(t_train));
stepplot(tf_system*0.8)

% Plot of the 2 system responses
figure(1)
clf
subplot(2,2,1)
plot(t_train, u_train)
grid on
title('Input for the model training u_{train}')
ylabel('Volt')
xlabel('Time[s]')
xlim([0, t_train(end)])
subplot(2,2,2)
plot(t_train, y_train)
grid on
title('System Response for the Input u_{train}')
ylabel('Amplitude')
xlabel('Time[s]')
xlim([0, t_train(end)])
subplot(2,2,3)
plot(t_test, u_test)
grid on
title('Input for the model testing u_{test}')
ylabel('Volt')
xlabel('Time[s]')
xlim([0, t_test(end)])
subplot(2,2,4)
plot(t_test, y_test)
grid on
title('System Response for the Input u_{test}')
ylabel('Amplitude')
xlabel('Time[s]')
xlim([0, t_test(end)])

%Save to Files
names = {'Time', 'Input', 'Position'};
output1 = table(t_train, u_train, y_train, 'VariableNames', names);
output2 = table(t_test, u_test, y_test, 'VariableNames', names);
writetable(output1, 'StepResultTrain.csv')
writetable(output2, 'StepResultTest.csv')
