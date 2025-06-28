%%                   Coding Lecture 3: OLG Models with Pensions in MATLAB                          
%                       Master in Economics and Finance
%   Department of Economics and Management "Marco Fanno" - University of Padova
%                     Advanced Public Economics - Fall 2023 

%% OLG two-period model with PAYG public pension system and inelastic labour supply
% REMARK: Put the file interest.m, wage.m, ksteady.m and kdyn.m in the working directory before running the code

clear;
clc;
addpath('C:\Users\Samet\Desktop\MATLAB.code'); % Change path as necessary

%% 1. Define Model's Parameters (INPUT)

t = 20;          % Number of transaction periods (generations)
n = 0.1;         % Population growth
beta = 0.40;     % Discount factor for household's utility in the two periods
alpha = 0.36;    % Capital share in production function
ls = 0.3;        % Steady state labour supply
nu0 = 257.15;    % Disutility parameters
nu1 = 3.33;      % 1/v1 in the utility function
tau = 0.5;       % Replacement rate for initial steady state analysis

% Save all model parameters in the same structure
par.t = t;
par.n = n;
par.beta = beta;
par.alpha = alpha;
par.ls = ls;
par.nu0 = nu0;
par.nu1 = nu1;
par.tau = tau;

% Build a vector of potential values for K to iterate the non-linear function and find solutions
kmin = 0.001; % Define the minimum
kmax = 5.0;   % Define the maximum
nk = 1000;    % Define the number of observations
k = linspace(kmin, kmax, nk); % Build the vector

% Initialize vector of solutions
yk = zeros(nk, 1);

for i = 1:nk
    par.i = i;
    yk(i) = ksteady(k(i), par); % Use built-in function ksteady.m
end

% Retain only positive solutions
yk = abs(yk);

[M, I] = min(yk); % Find the minimum solution
kinitial = k(I);  % Initial value for steady state capital

% Steady state capital and associated values
kss = fsolve(@(x) ksteady(x, par), kinitial);
yss = kss^alpha * ls^(1 - alpha);

c1ss = 1 / (1 + beta) * (wage(kss, ls, par) * ls + (n - interest(kss, ls, par)) / (1 + interest(kss, ls, par)) * tau * wage(kss, ls, par) * ls);
c2ss = beta * c1ss * (1 + interest(kss, ls, par));
utilss = log(c1ss) + beta * log(c2ss) - nu0 * ls^(1 + nu1) / (1 + nu1);

%% Scenario 2: Reduce Replacement Rate (tau = 0.2)

tau = 0.2;
par.tau = tau;

for i = 1:nk
    par.i = i;
    yk(i) = ksteady(k(i), par);
end

yk = abs(yk);
[M, I] = min(yk);
kinitial = k(I);

kssd = fsolve(@(x) ksteady(x, par), kinitial);
yssd = kssd^alpha * ls^(1 - alpha);

c1ssd = 1 / (1 + beta) * (wage(kssd, ls, par) * ls + (n - interest(kssd, ls, par)) / (1 + interest(kssd, ls, par)) * tau * wage(kssd, ls, par) * ls);
c2ssd = beta * c1ssd * (1 + interest(kssd, ls, par));
utilssd = log(c1ssd) + beta * log(c2ssd) - nu0 * ls^(1 + nu1) / (1 + nu1);

% Welfare impact (consumption equivalent change)
deltak = kssd - kss;
cec = exp((utilssd - utilss) / (1 + beta)) - 1;

%% Scenario 3: Dynamic Analysis (tau changes from 0.2 to 0.5)

kt = zeros(t + 1, 1);
utilt = kt;
periods = linspace(0, t, t + 1);
kt(1) = kssd;
utilt(1) = utilssd;
cect = zeros(t + 1, 1);
cect(1) = exp((utilt(1) - utilss) / (1 + beta)) - 1;

tau0 = 0.2;
par.tau0 = tau0;
tau1 = 0.5;
par.tau1 = tau1;

for i = 1:t
    if i == 2
        tau1 = 0.5; % Increase replacement rate for future generations
        par.tau1 = tau1;
    end
    k0 = kt(i);
    par.k0 = k0;
    k1 = fsolve(@(x) kdyn(x, par), k0);
    d0 = tau0 * wage(k0, ls, par) * ls;
    d1 = tau1 * wage(k1, ls, par) * ls;

    c1 = 1 / (1 + beta) * (wage(k0, ls, par) * ls - d0 + d1 * (1 + n) / (1 + interest(k1, ls, par)));
    c2 = beta * c1 * (1 + interest(k1, ls, par));
    utilt(i + 1) = log(c1) + beta * log(c2) - nu0 * ls^(1 + nu1) / (1 + nu1);
    cect(i + 1) = exp((utilt(i + 1) - utilss) / (1 + beta)) - 1;
    kt(i + 1) = k1;
end

%% Plot Results

figure(1);
plot(periods, kt, 'b', 'LineWidth', 2.5);
xlabel('Period t', 'FontSize', 18);
ylabel('Capital K_t', 'FontSize', 18);
title('Capital Dynamics: PAYG Pension System', 'FontSize', 20);

grid on;

figure(2);
plot(periods, cect * 100, 'r', 'LineWidth', 2.5);
xlabel('Period t', 'FontSize', 18);
ylabel('Welfare Change (%)', 'FontSize', 18);
title('Welfare Dynamics: PAYG Pension System', 'FontSize', 20);

grid on;
disp(['Steady-State Capital: ', num2str(kss)]);
disp(['Steady-State Utility: ', num2str(utilss)]);

tau = 0.2;
par.tau = tau;
kssd = fsolve(@(x) ksteady(x, par), kinitial);
c1ssd = 1 / (1 + beta) * (wage(kssd, ls, par) * ls + ...
        (n - interest(kssd, ls, par)) / (1 + interest(kssd, ls, par)) * ...
        tau * wage(kssd, ls, par) * ls);
c2ssd = beta * c1ssd * (1 + interest(kssd, ls, par));
utilssd = log(c1ssd) + beta * log(c2ssd) - nu0 * ls^(1 + nu1) / (1 + nu1);
cec = exp((utilssd - utilss) / (1 + beta)) - 1;
disp(['New Steady-State Capital (tau = 0.2): ', num2str(kssd)]);
disp(['New Steady-State Utility (tau = 0.2): ', num2str(utilssd)]);
disp(['Welfare Change (CEC): ', num2str(cec * 100), '%']);


