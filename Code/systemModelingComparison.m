%% System Modeling (Comparing Models)
% Created for the Multivariable Systems Lab
% Athens 02/2023
% Authors: 
%   G.Kassavetakis AM 02121203
%   G.Krommydas    AM 02121208

clc
clear
close all

%% Training the Tree Parameter Model
Tbl = readtable('StepResultTrain.csv');
t = Tbl.Time;
y = Tbl.Position;
%u=Tbl.Input;

% Finding the Max Derivative point
dy = gradient(y)./gradient(t);
[dy_max, index] = max(dy);
if length(index) > 1
    index = index(1);
end
t_max = t(index);
y_line = y(index)+dy_max*(t-t_max);

% Finding the y>0 and 0.63Kss indexes
t_index = find(y_line<=0, 1, 'last' );
ind = find(y<=0.63*y(end), 1, 'last' );

% Display of the Resulting parameters of the model
K = y(end);
L = t(t_index);
Tav = t(ind);
T = Tav - L;

% Model's tf Creation
s = tf('s');
G_1 = K/(T*s + 1)*exp(-L*s);

%% Training the Ordinary Least Squares Model
Tbl=readtable('SineResultTrain.csv');
t = Tbl.Time;
y = Tbl.Position;
u = Tbl.Input;

%Filter Creation for n=5
s = tf('s');
n = 5;
m = 5;
% L = [1.25,0.625,0.1563,0.01953,0.0009766];
% filter = tf(1,[1, L]);
filter = 1/(s + 3)^n;
L = filter.Denominator{:};
L = L(2:end);
Hy = -[s^4;s^3;s^2;s;1]*filter;
Hu = [s^4;s^3;s^2;s;1]*filter;

% Backstepping Vector calculation
phi1 = lsim(Hy,y,t);
phi2 = lsim(Hu,u,t);
phi = [phi1,phi2]';

% Method Implementation
N = length(t);
S1 = double(zeros(n+m,n+m));
S2 = double(zeros(n+m,1));
for i = 1:N
   S1 = S1 + double(phi(:,i)*phi(:,i)');
   S2 = S2 + double(phi(:,i)*y(i));
end

%Theta Calculation
theta0 = double(pinv(1/N*S1)*(1/N)*S2);
a = theta0(1:n)'+L;
b = theta0((n+1):end)';

% Model's tf Creation
G_OLS = tf(b,[1, a]);

%% Comparison using the Testing Measurement

% Reading Test Data
Tbl = readtable('SineResultTest.csv');
t_test = Tbl.Time;
y_test = Tbl.Position;
u_test = Tbl.Input;

y_1 = lsim(G_1,u_test,t_test);
y_OLS = lsim(G_OLS,u_test,t_test);

figure(1)
clf
plot(t_test,y_test,'r-')
hold on
plot(t_test,y_OLS,'b--')
plot(t_test,y_1,'g-.')
grid minor
title('OLS Model vs 3 Parameter Model')
xlabel('Time [s]')
ylabel('Amplitude')
xlim([min(t_test),max(t_test)])
legend('Plant','OLS Model','First Order Model','Location','southeast')

figure(2)
clf
plot(t_test,y_test-y_OLS,'b-')
hold on
plot(t_test,y_test-y_1,'g-')
grid minor
title('Estimation Error Comparison')
xlabel('Time [s]')
ylabel('Amplitude')
xlim([min(t_test),max(t_test)])
legend('OLS Model Error','First Order Model Error','Location','southeast')
