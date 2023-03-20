%% System Modeling (Step)
% Created for the Multivariable Systems Lab
% Athens 02/2023
% Authors: 
%   G.Kassavetakis AM02121203
%   G.Krommidas    AM02121208

clc
clear
close all

%% Reading the Training Result
Tbl=readtable('StepResultTrain.csv');
t=Tbl.Time;
y=Tbl.Position;
u=Tbl.Input;
%% Step Response Model Method

% Finding the Max Derivative point
dy=gradient(y)./gradient(t);
[dy_max,index]=max(dy);
if length(index)>1
    index=index(1);
end
t_max=t(index);
y_line=y(index)+dy_max*(t-t_max);

% Finding the y>0 and 0.63Kss indexes
t_index=find(y_line<=0, 1, 'last' );
ind=find(y<=0.63*y(end), 1, 'last' );

% Figure Analysis
figure(1)
clf
plot(t,y,'b-')
hold on
plot(t,min(y_line,y(end)),'r--')
plot(t,y(end)*ones(size(t)),'r--')
xline(t(ind),'k--')
scatter(t(ind),y(ind),'kx','LineWidth',2)
grid minor
title('Step Response Analysis')
xlabel('Time [s]')
ylabel('Amplitude')
legend('Plant','Max Derivative Line','y=Kss','y=0.63*Kss','Location','southeast')
xlim([min(t),max(t)])

% Display of the Resulting parameters of the model
K=y(end);
L=t(t_index);
Tav=t(ind);
T=Tav-L;
disp(['The K Parameter of the step Response is: ',num2str(K)])
disp(['The L Parameter of the step Response is: ',num2str(L)])
disp(['The T Parameter of the step Response is: ',num2str(T)])

%% Model Testing

% Model's tf Creation
s=tf('s');
G=K/(T*s+1)*exp(-L*s);

% Reading Test Data
Tbl=readtable('StepResultTest.csv');
t_test=Tbl.Time;
y_test=Tbl.Position;
u_test=Tbl.Input;
y_model=lsim(G,u_test,t_test);

% Comparison Figure
figure(2)
clf
plot(t_test,y_test,'b-')
hold on
plot(t_test,y_model,'g-')
yline(y_test(end),'r--')
grid minor
title('Model vs Test')
xlabel('Time [s]')
ylabel('Amplitude')
xlim([min(t_test),max(t_test)])
legend('Plant','Model','y_{ss}','Location','southeast')
