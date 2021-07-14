% Tested only on Matlab-Online
syms v T;

T_c = 545.5;
P_c = 4.83 * 10^6;
R = 8.314;
m = 0.480 + 1.574*0.338 - 0.176*(0.338^2);
% a = 0.42748*((1 + m*(1 - (T/T_c)^0.5))^2) * R^2 * T_c^2 / P_c;
b = 0.08664*R * T_c / P_c;

p = (R * T) / (v - b) - (0.42748 * ((1 + m*(1 - (T/T_c)^0.5))^2) * R^2 * T_c^2 / P_c)/(v * (v + b));
k = diff(p, v); % First derivative wrt v
k2 = diff(p, v, 2); % Second derivative wrt v

% Condition for critical point is k=0 and k2=0 simultaneously
[solt, solv] = solve(k == 0, k2 == 0);
V = vpa(solv);
T = vpa(solt);

disp(V(1)) % Gives the real root for critical volume
disp(T(1)) % Gives the real root for critical temperature