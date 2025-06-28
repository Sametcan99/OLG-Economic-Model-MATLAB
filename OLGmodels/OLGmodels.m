%%                   Coding Lecture 1: OLG Models in MATALB                         
%                       Master's Degree in Applied Economics
%   Department of Economics and Management "Marco Fanno" - University of Padova
%                     Advanced Public Economics - Fall 2024 

%% OLG two-period model: general Set-up (Dimond, 1965; Samuelson, 1958) 

% For all the assumptions on utility and production fuction of the numerical example see slide 4/14 of Coding Lecture 1

clear;
clc;

%% 1. Define Model's Parameters (INPUT)

T = 20;               % Number of periods (generations)
n = 0.1;              % Population growth
L = 30;               % Lenght of one generation period (years)
r = 0.03;             % Annual interest rate
beta = 0.4;     % Discount factor for houdehold's utility in the two periods -> It derives from the consumption smoothing assumption!
alpha = 0.36;         % Capital share in production function (Empirical data tells us that labour share is around 64% in modern economies)
k_ss = ((beta/(1+beta))*((1-alpha)/(1+n)))^(1/(1-alpha)); % Steady state capital
k_0 = k_ss/3;         % Initial capital stock

%% 2. Define initialize variables of the model


k = zeros(T, 1);  % Build the time series of capital stock
k(1,1)=k_0;       % Imposing first observation of capital stock at time zero equal to initial capital stock


%% 2.1 Iterate capital dynamics according to optimal solutions of the model

for t = 1:T

    if t < T
        k(t+1) = (beta/(1+beta))*((1-alpha)/(1+n))* k(t).^alpha  %This is the optimal solution for capital dyanmics - Equation (1)
    end
end 

%% 2.2 Compute dynamics of the other model's variable accordin to general equilibrium's optimal solutions

    y = k .^alpha;                     % Compute output dynamics - Equation (2)

    w = (1 - alpha) * k .^alpha;       % Compute wage dynamics - Equation (3)

    r = alpha * k .^ (alpha-1);        % Compute interest rate dynamics - Equation (4)

    s = ((beta)/(1+beta)) * w;         % Compute savings for the young generation dynamics - Equation (5) 
  
    c1 = w - s;                        % Compute consumption for young generation dynamics - Equation (6) 

    c2 = (1+r(2:end,:)) .* s(1:19,:);  % Compute consumption for old generation - Equation (7)  
 

%% 3. Result's Plot

t1 = 1:T; %define time horizon for young generation
t2 = 2:T; %define time horizon for old generation

%Capital dynamics - Figure 1

figure(1); %figure number to save
plot(t1, k, '-ob', 'LineWidth', 1.5, 'MarkerSize', 6); %plot(time horizon, time series we want to plot, colour, options);
xlabel('t','FontSize',20);
ylabel('k_{t}','FontSize',20);
xlim([1, 20]);
ylim([0.0209 ,0.065]);
legend('Capital','Fontsize',20);
title('Diamond-Samuelson OLG Model - Dynamics of capital stock','FontSize',20);
grid on;

%Output dynamics - Figure 2

figure(2);
plot(t1, y, '-ok', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('t','FontSize',20);
ylabel('y_{t}','FontSize',20);
xlim([1, 20]);
legend('Output','Fontsize',20);
title('Diamond-Samuelson OLG Model - Dynamics of output','FontSize',20);
grid on;

%Wage dynamics - Figure 3

figure(3);
plot(t1, w, '-or', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('t','FontSize',20);
ylabel('w_{t}','FontSize',20);
xlim([1, 20]);
legend('Wage','Fontsize',20);
title('Diamond-Samuelson OLG Model - Dynamics of wage','FontSize',20);
grid on;

%Rate dynamics - Figure 4

figure(4);
plot(t1, r, '-or', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('t','FontSize',20);
ylabel('r_{t}','FontSize',20);
xlim([1, 20]);
legend('Interest rate','Fontsize',20);
title('Diamond-Samuelson OLG Model - Dynamics of interest rate','FontSize',20);
grid on;

%Savings dynamics - Figure 5

figure(5);
plot(t1, s, '-ob', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Time Period','FontSize',20);
ylabel('s_{t}','FontSize',20);
xlim([1, 20]);
legend('Savings','Fontsize',20);
title('Diamond-Samuelson OLG Model - Dynamics of savings','FontSize',20);
grid on;

% Consumption dynamics for young and old generation - Figure 6

figure(6);
plot(t1, c1, '-or', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(t2, c2, '-ob', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('t','FontSize',20);
ylabel('c^i_{t}    with i \in {1,2}','FontSize',20);
xlim([1, 20]);
legend('c^1_{t}','c^2_{t+1}','Fontsize',20);
title('Diamond-Samuelson OLG Model - Dynamics of consumption','FontSize',20);
grid on;

%% 4. Dynamics of capital in the numerical example - Plot of convergence 

% I define the function of capital's evolution over time - Equation (1)

K = 0:0.0001:0.1; %I define a reliable domain for the value of k_t; 
K1 = ((beta/(1+beta))*((1-alpha)/(1+n)))*K.^alpha; %I write down the specific optimal g(k_t) in the numerical example; 

% REMARK: Here I'm calling k_t K and k_{t+1} K1 to avoid the overwriting of the values and allow the user to run the whole code from the beginning

%I define the identity function (Y=X) which in our case stays for capital stability (k_t = k_{t+1}) 
X = linspace(0, 0.1, 20); 
Y = X;

% Calculate the intersection point (borgig coding)
intersection_x = fminbnd(@(x) abs(interp1(K, K1, x) - interp1(X, Y, x)), 0, 0.08);
intersection_y = interp1(K, K1, intersection_x);


%Here I simply write down the equation as string to put them in plots close to relative curve

equation1 = '$$k_{t+1} = k_{t}$$';
equation2 = '$$k_{t+1} = \frac{\beta}{1+\beta} \frac{1-\alpha}{1+n} k_{t}^\alpha$$';
x1_pos = 0.042;
y1_pos = 0.04;
x2_pos = 0.025;
y2_pos = 0.07;

%I plot both the functions to identify the steady state solution for capital

figure(7);
plot(K,K1,'r','Linewidth',2);
hold on; %this command is needed to obtain two different plots in the same figure
plot(Y,X,'b','Linewidth',2);
hold on;
plot(intersection_x, intersection_y,'-mo', 'MarkerSize', 10);
text(x1_pos, y1_pos, equation1, 'FontSize', 20,'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
text(x2_pos, y2_pos, equation2, 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
xlabel('k_{t}', 'FontSize',20);
ylabel('k_{t+1}','FontSize',20);
xlim([0, 0.08]);
ylim([0 ,0.08]);
title('Stability in OLG model - Capital dynamics','FontSize',20);
grid on;
%% Steady-State Computation for OLG Model

% Parameters
beta = 0.4;          % Discount factor
alpha = 0.36;        % Capital share in production
n = 0.1;             % Population growth rate

% Steady-State Capital (k_ss)
k_ss = ((beta / (1 + beta)) * ((1 - alpha) / (1 + n)))^(1 / (1 - alpha));

% Steady-State Wage (w_ss) and Interest Rate (r_ss)
w_ss = (1 - alpha) * k_ss^alpha;            % Wage at steady state
r_ss = alpha * k_ss^(alpha - 1);            % Interest rate at steady state

% Steady-State Savings (s_ss)
s_ss = (beta / (1 + beta)) * w_ss;

% Steady-State Consumption
c1_ss = w_ss - s_ss;                        % Consumption for young generation
c2_ss = (1 + r_ss) * s_ss;                  % Consumption for old generation

% Steady-State Utility (U_ss)
U_ss = log(c1_ss) + beta * log(c2_ss);      % Utility at steady state

% Display Results
fprintf('Steady-State Results:\n');
fprintf('k_ss = %.5f\n', k_ss);
fprintf('w_ss = %.5f\n', w_ss);
fprintf('r_ss = %.5f\n', r_ss);
fprintf('s_ss = %.5f\n', s_ss);
fprintf('c1_ss = %.5f\n', c1_ss);
fprintf('c2_ss = %.5f\n', c2_ss);
fprintf('U_ss = %.5f\n', U_ss);
%% EXTENSION 5. Two-period OLG with technological progress
% What will change if we introduce dynamic technological progress? 

%% 5.1. Define Model's Parameters (INPUT)

gamma = 0.007                %exogenous grwoth rate of techonolgical progress - empirically driven from TFP's growth data
delta = 0.05                  %depreciation rate of capital (typically range between 5% and 10%
k_ss_tec = ((beta/(1+beta))*((1-alpha)/((1+n)*(1+gamma))))^(1/(1-alpha)); % steady state capital
k_tec_0 = k_ss_tec/3;         % Initial capital stock

%% 5.2. Build and iterate the dynamics of technological progress assuming it is growing only for 10 generations


A = zeros(T, 1);  %build the time series for technological progress
A_0 = 1.2;      %assume a value for the initial observation 
A(1,1)=A_0;      %impose the initial observation of the vector equal to the assumption




J = 15;  %horizon until which A is assumed to grow (15 periods), afterwards is assumed fixed 

for t = 1:T

    if t <= J
        A(t+1) = (1+gamma)*A(t) 
    else
         A(t) = A(t-1) 
    end
end
   


%% %5.3. Iterate capital dynamics according to optimal solutions of the model

k_tec = zeros(T, 1);  % Build the time series of capital stock
k_tec(1,1)=k_tec_0;   % Imposing first observation of capital stock at time zero equal to initial capital stock

for t = 1:T

    if t < T
        k_tec(t+1) = (beta/(1+beta))*((1-alpha)/((1+n)*(1+gamma)))* k_tec(t).^alpha  %This is the optimal solution for capital dyanmics - Equation (8)
    end
end 

%% 5.4. Compute dynamics of the other model's variable accordin to general equilibrium's optimal solutions

    y_tec = k_tec .^alpha;                     % Compute output dynamics - Equation (9)

    w_tec = (1 - alpha) * k_tec .^alpha;       % Compute wage dynamics - Equation (10)

    r_tec = alpha * k_tec .^ (alpha-1) - delta;        % Compute interest rate dynamics - Equation (11)

    s_tec = ((beta)/(1+beta))* w_tec .* A;         % Compute savings for the young generation dynamics - Equation (12) 
  
    c1_tec = w_tec - s_tec;                        % Compute consumption for young generation dynamics - Equation (13) 

    c2_tec = (1+r_tec(2:end,:)) .* s_tec(1:19,:);  % Compute consumption for old generation - Equation (14)  

%% 5.5. Comparing Result's Plot

t1 = 1:T;
t2 = 2:T;

%Capital dynamics - Figure 1

figure(8);
plot(t1, k, '-ob', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(t1, k_tec, '-or', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('t','FontSize',20);
ylabel('k_{t}','FontSize',20);
xlim([1, 20]);
ylim([0.0209 ,0.065]);
legend('Capital','Capital with technological progress','Fontsize',20);
title('Diamond-Samuelson OLG Model - Dynamics of capital stock','FontSize',20);
grid on;

%Output dynamics - Figure 2

figure(9);
plot(t1, y, '-ok', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(t1, y_tec, '-or', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('t','FontSize',20);
ylabel('y_{t}','FontSize',20);
xlim([1, 20]);
legend('Output','Output with technological progress','Fontsize',20);
title('Diamond-Samuelson OLG Model - Dynamics of output','FontSize',20);
grid on;

%Wage dynamics - Figure 3

figure(10);
plot(t1, w, '-or', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(t1, w_tec, '-ob', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('t','FontSize',20);
ylabel('w_{t}','FontSize',20);
xlim([1, 20]);
legend('Wage','Wage with technological progress','Fontsize',20);
title('Diamond-Samuelson OLG Model - Dynamics of wage','FontSize',20);
grid on;

%Rate dynamics - Figure 4

figure(11);
plot(t1, r, '-or', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on; 
plot(t1, r_tec, '-ob', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('t','FontSize',20);
ylabel('r_{t}','FontSize',20);
xlim([1, 20]);
legend('Interest rate','Interest rate with technological progress','Fontsize',20);
title('Diamond-Samuelson OLG Model - Dynamics of interest rate','FontSize',20);
grid on;

%Savings dynamics - Figure 5

figure(12);
plot(t1, s, '-or', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(t1, s_tec, '-ob', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Time Period','FontSize',20);
ylabel('s_{t}','FontSize',20);
xlim([1, 20]);
legend('Savings','Savings with technological progress','Fontsize',20);
title('Diamond-Samuelson OLG Model - Dynamics of savings','FontSize',20);
grid on;

% Consumption dynamics for young and old generation - Figure 6

figure(13);
plot(t1, c1, '-or', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(t2, c2, '-ob', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on; 
plot(t1, c1_tec, '-ok', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(t2, c2_tec, '-oy', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('t','FontSize',20);
ylabel('c^i_{t}    with i \in {1,2}','FontSize',20);
xlim([1, 20]);
legend('c^1_{t}','c^2_{t+1}','c^1_{t} with technological progress','c^2_{t+1} with technological progress','Fontsize',12);
title('Diamond-Samuelson OLG Model - Dynamics of consumption','FontSize',20);
grid on;
