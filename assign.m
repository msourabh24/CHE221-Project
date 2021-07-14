% Assignment - CHE221

% System: Acetonitrile
% EOS: Soave-Redlich-Kwong

% Parameters: P_c, T_c, R, omega
P_c = 4.83 * 10^6;
T_c = 545.5;
R = 8.314;
omega = 0.338;

% Define a,b,alpha
alpha = @(T) (1 + (0.480 + 1.574*omega - 0.176*omega^2)*(1 - (T./T_c).^0.5)).^2;
a = @(T) (0.42748 * alpha(T) * R^2 * T_c^2) / P_c;
b = 0.08664 * R * T_c / P_c;

% Expression for P_sat
C = [58.308 -5385.6 -5.4954 5.3634*10^-6 2];
P_sat = @(T) exp(C(1) + (C(2) / T) + C(3) * log(T) + C(4) * (T^C(5)));

% Temperature range from T = 300K to 600K with increments of 20K
T = 300:10:600;
A = zeros(31);
P_s = zeros(25);
% Calculate guess P_sat
for i = 1:31
    A(i) = a(T(i));
    if (i < 26)
        P_s(i) = P_sat(T(i));
    end
end

% Define P, Z and mu
P = @(V,T) R.*T./(V-b) - a(T)./(V.*(V+b));
Z = @(V,T) (V/(V-b)) - (a(T) / (R*T*(V+b)));
mu = @(V,T) log(V/(V-b)) + a(T)*log(V/(V+b))/(b*R*T) + Z(V,T) - 1 - log(Z(V,T));

tol = 0.01;
del_mu = zeros(25);
v = zeros(25,3);
for i = 1:25
    % Calculate initial values of V liquid and vapor
    v(i,:) = roots([P_s(i) -R*T(i) -(P_s(i)*b*b + R*T(i)*b - A(i)) -A(i)*b]);

    % Calculate difference in mu l and mu v
    del_mu(i) = mu(v(i,3),T(i)) - mu(v(i,1),T(i));
    % disp(del_mu(i));

    if (del_mu(i) > 0)
        k = 1;
    else
        k = -1;
    end

    while (abs(del_mu(i)) > tol)
        P_s(i) = P_s(i) + k * 10^2; % Update the P_sat value by incrementing or decrementing by 10Pa
        v(i,:) = roots([P_s(i) -R*T(i) -(P_s(i)*b*b + R*T(i)*b - A(i)) -A(i)*b]);
        del_mu(i) = mu(v(i,3),T(i)) - mu(v(i,1),T(i));
        % disp(del_mu(i));
    end
end

% Derived from EOS from critical.m
v_c = sqrt(3*a(T_c)*b / (R*T_c));
p1 = P(v_c,T_c);

hold on;

% Plot the dome for vapor-liquid equilibria
xx = [v(:,3)' v_c flip(v(:,1)')];
yy = [P_s(:,1)' p1 flip(P_s(:,1)')];

plot(xx(1:3:25),yy(1:3:25),'o','markerfacecolor','black');
plot(xx(26),yy(26),'o','markerfacecolor','black');
plot(xx(27:3:51),yy(27:3:51),'o','markerfacecolor','black');

plot(xx(1:25),yy(1:25),'linestyle','-.','linewidth',2,'color','black','Marker','none');
plot(xx(27:51),yy(27:51),'linestyle','-.','linewidth',2,'color','black','Marker','none');
xspline = linspace(xx(25),xx(27),200);
yspline = interp1(xx,yy,xspline,'spline');
plot(xspline,yspline,'linestyle','-.','linewidth',2,'color','black');

% Plot isotherms for T = [300 330 360 390 420 450 480 510 540]
for i = 1:3:25
    p = @(V) (R.*T(i)./(V-b) - (A(i)./(V.*(V+b))));
    fplot(p,[b+10^-5 3*10^-3],'linewidth',1.5);
end

% Plot at T_c = 545.5
p = @(V) P(V,545.5);
fplot(p,[b+10^-5 3*10^-3],'linewidth',1.5);

% Plot isotherms for T = [550 580]
for i = 26:3:31
    p = @(V) (R.*T(i)./(V-b) - (A(i)./(V.*(V+b))));
    fplot(p,[b+10^-5 3*10^-3],'linewidth',1.5);
end

% Scale the plot to display the curves properly
grid on;
xlim([0 3*10^-3]);
ylim([0, 70 * 10^5]);

xlabel("Volume (in m^3/mol)",'fontsize',15);
ylabel("Pressue (in Pa)",'fontsize',15);

hold off;