%% System Modeling (Ordinary Least Squares)
% Created for the Multivariable Systems Lab
% Athens 02/2023
% Authors: 
%   G. Kassavetakis AM 02121203
%   G. Krommydas    AM 02121208

clc
clear
close all

%% Reading the Training Result
Tbl = readtable('SineResultTrain.csv');
t = Tbl.Time;
y = Tbl.Position;
u = Tbl.Input;

%% Ordinary Least Squares Model Method

%Filter Creation for n=5
s = tf('s');
n = 5;
m = 4;
% L = [1.25, 0.625, 0.1563, 0.01953, 0.0009766];
% filter = tf(1, [1, L]);
filter = 1/(s+5)^n;
L = filter.Denominator{:};
L = L(2:end);
%Best Result
Hy = -[s^4;s^3;s^2;s;1]*filter;
Hu = [s^4;s^3;s^2;s;1]*filter;

% Backstepping Vector calculation
phi1 = lsim(Hy,y,t);
phi2 = lsim(Hu,u,t);
phi = [phi1, phi2]';

% Method Implementation
N = length(t);
S1 = double(zeros(n+m+1,n+m+1));
S2 = double(zeros(n+m+1,1));
for i = 1:N
   S1 = S1 + double(phi(:,i)*phi(:,i)');
   S2 = S2 + double(phi(:,i)*y(i));
end

%Theta Calculation
theta0 = double(pinv(1/N*S1)*(1/N)*S2);
a = theta0(1:n)'+L;
b = theta0((n+1):end)';

% Model's tf Creation
G = tf(b,[1, a]);

% Training Result Figure
y_m_train = lsim(G,u,t);
figure(1)
clf
plot(t,y,'b-')
hold on
plot(t,y_m_train,'g-.')
title('Model vs Training Measurement')
xlabel('Time [s]')
ylabel('Amplitude')
xlim([min(t),max(t)])
legend('Plant','Model','Location','northeast')

%% Model Testing

% Reading Test Data
Tbl = readtable('SineResultTest.csv');
t_test = Tbl.Time;
y_test = Tbl.Position;
u_test = Tbl.Input;
y_model = lsim(G,u_test,t_test);

% Testing Figure
figure(2)
clf
plot(t_test,y_test,'b-')
hold on
plot(t_test,y_model,'g-.')
grid minor
title('Model vs Testing Measurement')
xlabel('Time [s]')
ylabel('Amplitude')
xlim([min(t_test),max(t_test)])
legend('Plant','Model','Location','southeast')

figure(3)
clf
plot(t_test,y_test-y_model,'r-')
grid minor
title('Estimation Error on Testing Measurement')
xlabel('Time [s]')
ylabel('Amplitude')
xlim([min(t_test),max(t_test)])

%% Alternative no Zero Model
% 
% %Creating new model without zeros
% Den=G.Denominator{:};
% G_new=tf(dcgain(G)*Den(end),Den);
% 
% y_new_train=lsim(G_new,u,t);
% y_new_test=lsim(G_new,u_test,t_test);
% % Testing Figure
% figure(4)
% clf
% plot(t,y,'b-')
% hold on
% plot(t,y_new_train,'g-.')
% grid minor
% title('Model (no zeros) vs Training Measurement')
% xlabel('Time [s]')
% ylabel('Amplitude')
% xlim([min(t),max(t)])
% legend('Plant','Model','Location','southeast')
% 
% figure(5)
% clf
% plot(t_test,y_test,'b-')
% hold on
% plot(t_test,y_new_test,'g-.')
% grid minor
% title('Model (no zeros) on Testing Measurement')
% xlabel('Time [s]')
% ylabel('Amplitude')
% xlim([min(t_test),max(t_test)])
% legend('Plant','Model','Location','southeast')
% 
% figure(6)
% clf
% plot(t_test,y_test-y_model,'r-')
% hold on
% plot(t_test,y_test-y_new_test,'b-')
% grid minor
% title('Estimation Error on Testing Measurement')
% xlabel('Time [s]')
% ylabel('Amplitude')
% xlim([min(t_test),max(t_test)])
% legend('Model (plus zeros)','Model (no zeros)','Location','southeast')