close all; clear all; clc;

global C_La C_D0 eta c rho T
% table 1. Lift, Drag coefficients and eta vs. Mach number
Ma_in = 0.0:0.4:3.2;
C_La_out = [2.24, 2.325, 2.35, 2.29, 2.16, 1.95, 1.70, 1.435, 1.25];
C_D0_out = [0.0065, 0.0055, 0.006, 0.0118, 0.0111, 0.0086, 0.0074, 0.0069, 0.0068];
C_La = griddedInterpolant(Ma_in, C_La_out, 'spline', 'nearest');
C_D0 = griddedInterpolant(Ma_in, C_D0_out, 'spline', 'nearest');
Ma_in = [0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0];
eta_out = [0.54, 0.75, 0.79, 0.78, 0.89, 0.93, 0.97, 1.00]; 
eta = griddedInterpolant(Ma_in, eta_out, 'spline', 'nearest');
% table 2. Speed of sound, Air density vs. Altitude
h_in = [0, 5000, 10000, 15000, 20000, 25000, 30000, 36090, 40000, 45000, 50000, 55000, 60000, 70000, 80000, 82020, 90000, 100000];
c_out = [1116, 1097, 1077, 1057, 1037, 1016, 994.7, 968.1, 968.1, 968.1, 968.1, 968.1, 968.1, 968.1, 968.1, 968.1, 984.2, 1004];
rho_out = [2377, 2048, 1755, 1496, 1266, 1065, 889.3, 706.1, 585.1, 460.1, 361.8, 284.5, 223.8, 138.4, 85.56, 77.64, 51.51, 31.38] * 1e-6;
c = griddedInterpolant(h_in, c_out, 'spline', 'nearest');
rho = griddedInterpolant(h_in, rho_out, 'spline', 'nearest');
% table 3. Thrust
Ma_in = (0.0:0.4:3.2).' * ones([1,12]);
h_in = ones([9,1]) * [0, 5000, 15000, 25000, 35000, 45000, 55000, 65000, 75000, 85000, 95000, 105000];
T_out = [23.3,20.6,15.4,9.9,5.8,2.9,1.3,0.7,0.3,0.1,0.1,0.0;
    22.8,19.8,14.4,9.9,6.2,3.4,1.7,1.0,0.5,0.3,0.1,0.1;
    24.5,22.0,16.5,12.0,7.9,4.9,2.8,1.6,0.9,0.5,0.3,0.2;
    29.4,27.3,21.0,15.8,11.4,7.2,3.8,2.7,1.6,0.9,0.6,0.4;
    29.7,29.0,27.5,21.6,15.7,10.5,6.5,3.8,2.3,1.4,0.8,0.5;
    29.9, 29.4,28.4,26.6,21.2,14.0,8.7,5.1,3.3,1.9,1.0,0.5;
    29.9,29.3,28.4,27.1,25.6,17.2,10.7,6.5,4.1,2.3,1.2,0.5;
    29.8,29.1,28.2,26.8,25.6,20.0,12.2,7.6,4.7,2.8,1.4,0.5;
    29.7,28.9,27.5,26.1,24.9,20.3,13.0,8.0,4.9,2.8,1.4,0.5] * 1e3;
T = griddedInterpolant(Ma_in, h_in, T_out, 'spline', 'nearest');

global alt epsilon S mg
epsilon = 3*pi/180; S = 500; mg = 34200;
altitudes = 0:2500:47500;
diff = 1e-6; % finite difference
%% Problem statement
% Find gamma that maximizes the rate of climb vsin(gamma) for different 
% altitudes, under the constraint f(x, u) = 0 (x=[v,alpha], u=[gamma]):
% V_1(v,gamma,alpha) = 0 & V_2(v,gamma,alpha) = 0

%% (a) steepest descent
learning_rate_x = 1e-2; learning_rate_u = 1e-5;

x0 = [1000; 5*pi/180]; u0 = 15*pi/180;
x = x0; u = u0; L = -Inf;
x_save = zeros(length(altitudes),2);
u_save = zeros(length(altitudes),1);
L_save = zeros(length(altitudes),1);
total_x_iters = 0;
for i=1:length(altitudes)
    alt = altitudes(i);
    iter = 0;
    while iter < 100
        % x update
        x_iter = 0;
        f = f_val(x, u);
        while x_iter < 10000
            f_x = [0, 0; 0, 0];
            f_x(:,1) = (f_val(x + diff*[1; 0], u) - f)/diff;
            f_x(:,2) = (f_val(x + diff*[0; 1], u) - f)/diff;
            x = x - learning_rate_x*(f_x\f);
            f = f_val(x, u);
            x_iter = x_iter + 1;
            if norm(f) < mg*0.001
                fprintf('x_iter: %d \n', x_iter);
                break;
            end
        end
        total_x_iters = total_x_iters + x_iter;
        % maximum L that satisfies the constraint
        L_prev = L;
        L = L_val(x, u);
        iter = iter + 1;
        if abs(L - L_prev) < 1
            fprintf('altitude: %d, velocity: %f, rate-of-climb: %f \n', alt, x(1), L);
            break;
        end
        % u update
        f_x = [0, 0; 0, 0];
        f_x(:,1) = (f_val(x + diff*[1; 0], u) - f)/diff;
        f_x(:,2) = (f_val(x + diff*[0; 1], u) - f)/diff;
        lambda_T = -[sin(u), 0]/f_x;
        f_u = (f_val(x, u + diff) - f)/diff;
        L_u = x(1)*cos(u);
        u = u + learning_rate_u*(L_u + lambda_T*f_u);
    end
    x_save(i,:) = x.'; 
    u_save(i) = u;
    L_save(i) = L;
end

figure(1);
plot(altitudes, x_save(:,1));
xlabel('altitude [ft]');
ylabel('velocity [ft/s]');
figure(2);
plot(altitudes, L_save);
xlabel('altitude [ft]');
ylabel('rate of climb [ft/s]');
fprintf('total x iterations: %d', total_x_iters);

%% (b) newton-raphson method
learning_rate_x = 1e-2; learning_rate_u = 1e-2;

x0 = [1000; 5*pi/180]; u0 = 15*pi/180;
x = x0; u = u0; L = -Inf;
x_save = zeros(length(altitudes),2);
u_save = zeros(length(altitudes),1);
L_save = zeros(length(altitudes),1);
for i=1:length(altitudes)
    alt = altitudes(i);
    iter = 0;
    while iter < 100
        % x update
        x_iter = 0;
        f = f_val(x, u);
        while x_iter < 10000
            f_x = [0, 0; 0, 0];
            f_x(:,1) = (f_val(x + diff*[1; 0], u) - f)/diff;
            f_x(:,2) = (f_val(x + diff*[0; 1], u) - f)/diff;
            x = x - learning_rate_x*(f_x\f);
            f = f_val(x, u);
            x_iter = x_iter + 1;
            if norm(f) < mg*0.001
                fprintf('x_iter: %d \n', x_iter);
                break;
            end
        end
        % maximum L that satisfies the constraint
        L_prev = L;
        L = L_val(x, u);
        iter = iter + 1;
        if abs(L - L_prev) < 1
            fprintf('altitude: %d, velocity: %f, rate-of-climb: %f \n', alt, x(1), L);
            break;
        end
        % u update
        L_u = x(1)*cos(u);
        lambda_T = -[sin(u), 0]/f_x;
        f_u = (f_val(x, u + diff) - f)/diff;
        H_u = L_u + lambda_T*f_u;
        L_uu = -L_val(x, u);
        du = -inv(L_uu)*H_u;
        u = u + learning_rate_u*du;
    end
    x_save(i,:) = x.'; 
    u_save(i) = u;
    L_save(i) = L;
end

figure(1);
plot(altitudes, x_save(:,1));
xlabel('altitude [ft]');
ylabel('velocity [ft/s]');
figure(2);
plot(altitudes, L_save);
xlabel('altitude [ft]');
ylabel('rate of climb [ft/s]');

%% (c) conjugate gradient method
learning_rate_x = 1e-2; learning_rate_u = 1e-5;

x0 = [1000; 5*pi/180]; u0 = 15*pi/180;
x = x0; u = u0; L = -Inf;
x_save = zeros(length(altitudes),2);
u_save = zeros(length(altitudes),1);
L_save = zeros(length(altitudes),1);
total_x_iters = 0;
for i=1:length(altitudes)
    alt = altitudes(i);
    iter = 0;
    while iter < 100
        % x update
        x_iter = 0;
        f = f_val(x, u);
        f_x = [0, 0; 0, 0];
        f_x(:,1) = (f_val(x + diff*[1; 0], u) - f)/diff;
        f_x(:,2) = (f_val(x + diff*[0; 1], u) - f)/diff;
        while x_iter < 10000
            g_k = -f_x\f;
            if x_iter == 0
                d_k = g_k;
            else
                d_k = g_k + beta_k*d_k;
            end
            x = x + learning_rate_x*d_k;
            f = f_val(x, u);
            f_x = [0, 0; 0, 0];
            f_x(:,1) = (f_val(x + diff*[1; 0], u) - f)/diff;
            f_x(:,2) = (f_val(x + diff*[0; 1], u) - f)/diff;
            g_k_ = g_k; g_k = -f_x\f;
            delg_k = g_k - g_k_;
            beta_k = (delg_k.'*g_k)/(delg_k.'*d_k);
            x_iter = x_iter + 1;
            if norm(f) < mg*0.001
                fprintf('x_iter: %d \n', x_iter);
                break;
            end
        end
        total_x_iters = total_x_iters + x_iter;
        % maximum L that satisfies the constraint
        L_prev = L;
        L = L_val(x, u);
        iter = iter + 1;
        if abs(L - L_prev) < 1
            fprintf('altitude: %d, velocity: %f, rate-of-climb: %f \n', alt, x(1), L);
            break;
        end
        % u update
        f_x = [0, 0; 0, 0];
        f_x(:,1) = (f_val(x + diff*[1; 0], u) - f)/diff;
        f_x(:,2) = (f_val(x + diff*[0; 1], u) - f)/diff;
        lambda_T = -[sin(u), 0]/f_x;
        f_u = (f_val(x, u + diff) - f)/diff;
        L_u = x(1)*cos(u);
        u = u + learning_rate_u*(L_u + lambda_T*f_u);
    end
    x_save(i,:) = x.'; 
    u_save(i) = u;
    L_save(i) = L;
end

figure(1);
plot(altitudes, x_save(:,1));
xlabel('altitude [ft]');
ylabel('velocity [ft/s]');
figure(2);
plot(altitudes, L_save);
xlabel('altitude [ft]');
ylabel('rate of climb [ft/s]');
fprintf('total x iterations: %d', total_x_iters);

%% functions
% L(x,u)
function L = L_val(x, u)
    L = x(1)*sin(u);
end
% f(x,u)
function f = f_val(x, u)
    global C_La C_D0 eta c rho T
    global alt epsilon S mg
    Ma = x(1)/c(alt);
    L = 0.5*rho(alt)*x(1)^2*C_La(Ma)*x(2)*S;
    D = 0.5*rho(alt)*x(1)^2*(C_D0(Ma) + eta(Ma)*C_La(Ma)*x(2)^2)*S;
    f1 = T(x(1),alt)*cos(x(2)+epsilon) - D - mg*sin(u);
    f2 = T(x(1),alt)*sin(x(2)+epsilon) + L - mg*cos(u);
    f = [f1; f2];
end