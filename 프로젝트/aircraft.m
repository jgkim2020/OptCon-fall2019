clear all; close all; clc;

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

%% plot coefficients & max thrust
figure();
aoa = -15:30;
Ma_plot = [0.7 1.2];
for i = 1:length(Ma_plot)
    Cl_plot = C_La(Ma_plot(i))*(aoa/180*pi);
    Cd_plot = C_D0(Ma_plot(i)) + eta(Ma_plot(i))*C_La(Ma_plot(i))*(aoa/180*pi).^2;
    ax = gca;
    ax.ColorOrderIndex = i;
    plot(aoa, Cl_plot, '-');
    hold on;
    ax.ColorOrderIndex = i;
    plot(aoa, Cd_plot, '--');
    hold on;
end
xlabel('\alpha (deg)');
ylabel('$C_L, C_D$', 'interpreter', 'latex');

figure();
[Ma_plot, h_plot] = ndgrid(0:0.1:3, 0:5000:100000);
surf(Ma_plot, h_plot/1000, T(Ma_plot, h_plot)/1000);
xlabel(sprintf('$Ma$'), 'interpreter', 'latex');
ylabel(sprintf('$h$ (1000 ft)'), 'interpreter', 'latex');
zlabel(sprintf('$Thrust$ (1000 lbf)'), 'interpreter', 'latex');
