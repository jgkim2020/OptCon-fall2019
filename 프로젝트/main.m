clear all; close all; clc;

%% Tables
% i) lift, drag coefficients, eta vs. Mach number
Ma_in = 0.0:0.4:3.2;
C_La_out = [2.24, 2.325, 2.35, 2.29, 2.16, 1.95, 1.70, 1.435, 1.25];
C_D0_out = [0.0065, 0.0055, 0.006, 0.0118, 0.0111, 0.0086, 0.0074, 0.0069, 0.0068];
C_La = griddedInterpolant(Ma_in, C_La_out, 'spline', 'nearest');
C_D0 = griddedInterpolant(Ma_in, C_D0_out, 'spline', 'nearest');
Ma_in = [0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0];
eta_out = [0.54, 0.75, 0.79, 0.78, 0.89, 0.93, 0.97, 1.00]; 
eta = griddedInterpolant(Ma_in, eta_out, 'spline', 'nearest');
% ii) speed of sound, air density vs. altitude
h_in = [0, 5000, 10000, 15000, 20000, 25000, 30000, 36090, 40000, 45000, ...
    50000, 55000, 60000, 70000, 80000, 82020, 90000, 100000];
a_out = [1116, 1097, 1077, 1057, 1037, 1016, 994.7, 968.1, 968.1, 968.1, ...
    968.1, 968.1, 968.1, 968.1, 968.1, 968.1, 984.2, 1004];
rho_out = [2377, 2048, 1755, 1496, 1266, 1065, 889.3, 706.1, 585.1, 460.1, ...
    361.8, 284.5, 223.8, 138.4, 85.56, 77.64, 51.51, 31.38] * 1e-6;
rho_out = rho_out*32.174; % slug -> lb
a = griddedInterpolant(h_in, a_out, 'spline', 'nearest');
rho = griddedInterpolant(h_in, rho_out, 'spline', 'nearest');
% iii) maximum thrust vs. Mach number, altitude
Ma_in = (0.0:0.4:3.2).' * ones([1,12]);
h_in = ones([9,1]) * [0, 5000, 15000, 25000, 35000, 45000, 55000, 65000, ...
    75000, 85000, 95000, 105000];
T_out = [23.3,20.6,15.4,9.9,5.8,2.9,1.3,0.7,0.3,0.1,0.1,0.0;
    22.8,19.8,14.4,9.9,6.2,3.4,1.7,1.0,0.5,0.3,0.1,0.1;
    24.5,22.0,16.5,12.0,7.9,4.9,2.8,1.6,0.9,0.5,0.3,0.2;
    29.4,27.3,21.0,15.8,11.4,7.2,3.8,2.7,1.6,0.9,0.6,0.4;
    29.7,29.0,27.5,21.6,15.7,10.5,6.5,3.8,2.3,1.4,0.8,0.5;
    29.9, 29.4,28.4,26.6,21.2,14.0,8.7,5.1,3.3,1.9,1.0,0.5;
    29.9,29.3,28.4,27.1,25.6,17.2,10.7,6.5,4.1,2.3,1.2,0.5;
    29.8,29.1,28.2,26.8,25.6,20.0,12.2,7.6,4.7,2.8,1.4,0.5;
    29.7,28.9,27.5,26.1,24.9,20.3,13.0,8.0,4.9,2.8,1.4,0.5] * 1e3;
T_out = T_out*32.174; % lbf -> poundal
Tmax = griddedInterpolant(Ma_in, h_in, T_out, 'spline', 'nearest');

%% Setup
% [L, T, M] units = [ft, s, lb]
% state: [V, gamma, x, h] 
% control: [alpha, T]
% time scale: w

% constants
g = 32.174;
m = 34200;
S = 500;
epsilon = 3*pi/180;

% parameters
tau0 = 0; tauf = 60; % normalized time
dtau = 0.1;
taus = tau0:dtau:tauf;
Vtanker = 800; % KC-46 cruise speed: 777, max speed: 836
htanker = 36000; % KC-46 service ceiling: 40100
xtanker0 = 10000; % KC-46 initial x position
Voffset = 200; % fighter jet initial V offset
hoffset = 4000; % fighter jet initial h offset

% hyperparameters
c = 1; % c in [0: minimum fuel prob. ~ 1: minimum time prob.]
del = 1e-6; % step size for finite difference
maxIter = 1e3;
Ku = 1e-3; Ku_min = 1e-5; Ku_max = 5e-3;
Kw = 5e-5; Kw_min = 5e-7; Kw_max = 1e-4;
err_thresh = 0;

% terminal cost
h = @(s,w) (s(1)-Vtanker)^2/Vtanker + s(2)^2/0.002 + ...
    (s(3)-xtanker0-w*Vtanker*tauf)^2/xtanker0 + (s(4)-htanker)^2/htanker;
dhds = @(s,w) partialderiv(h, {s,w}, 1, del);
dhdw = @(s,w) partialderiv(h, {s,w}, 2, del);

% integration cost
L = @(s,u,w) w*(c + (1-c)*abs(u(2))/Tmax(s(1)/a(s(4)),s(4)));

% lift, drag
Lift = @(s,u) 0.5*rho(s(4))*s(1)^2*S*C_La(s(1)/a(s(4)))*u(1);
Drag = @(s,u) 0.5*rho(s(4))*s(1)^2*S*(C_D0(s(1)/a(s(4))) + ...
    eta(s(1)/a(s(4)))*C_La(s(1)/a(s(4)))*u(1)^2);

% final time normalized state dynamics
dsdtau = @(s,u,w) w*[(u(2)*cos(u(1)+epsilon) - Drag(s,u))/m - g*sin(s(2));
    (u(2)*sin(u(1)+epsilon) + Lift(s,u))/(m*s(1)) - g*cos(s(2))/s(1);
    s(1)*cos(s(2));
    s(1)*sin(s(2))];

% final time normalized Hamiltonian
H = @(s,u,lambda,w) L(s,u,w) + lambda.'*dsdtau(s,u,w);
dHds = @(s,u,lambda,w) partialderiv(H, {s,u,lambda,w}, 1, del);
dHdu = @(s,u,lambda,w) partialderiv(H, {s,u,lambda,w}, 2, del);
dHdlambda = @(s,u,lambda,w) partialderiv(H, {s,u,lambda,w}, 3, del);
dHdw = @(s,u,lambda,w) partialderiv(H, {s,u,lambda,w}, 4, del);

% temp
dsdtau_test = @(s,u,w) dHdlambda(s,u,[1;1;1;1],w).';

% final time normalized costate dynamics
dlambdadtau = @(s,u,lambda,w) -dHds(s,u,lambda,w).';

% boundary conditions
s0 = [Vtanker+Voffset; 0; 0; htanker+hoffset];
% sf = [Vtanker; 0; xtanker+Vtanker*w*tauf; htanker]
lambdaf = @(s,w) dhds(s,w).';

% control bounds
u_min = @(s) [-25*pi/180; 0];
u_max = @(s) [25*pi/180; Tmax(s(1)/a(s(4)),s(4))];

%% Trajectory optimization
% initial guess
u_init = [-1*pi/180*ones(1,length(taus)); 0*g*ones(1,length(taus))];
w_init = 1;
u = u_init; w = w_init;

% history
J_his = zeros(1,maxIter);
h_his = zeros(1,maxIter);
s_his = zeros(4,length(taus),maxIter);
u_his = zeros(2,length(taus),maxIter);
lambda_his = zeros(4,length(taus),maxIter);
w_his = zeros(1,maxIter);
Tmax_his = zeros(1,length(taus),maxIter);

deltau = Inf;
deltaw = Inf;
for iter = 1:maxIter
    % forward pass
    [~, s] = ode45(@(tau,s) statedynmcs(tau,s,taus,u,w,dsdtau_test), taus, s0);
    s = s.';
    ub = zeros(size(u));
    for i = 1:size(u,2)
        ub(:,i) = u_max(s(:,i));
    end
    % backward pass
    sf = s(:,end);
    [~, lambda] = ode45(@(tau,lambda) costatedynmcs(tau,lambda,taus,s,u,w,dlambdadtau), ...
        taus(end:-1:1), lambdaf(sf,w));
    lambda = lambda(end:-1:1,:).';
    
    % log
    J_his(iter) = perfidx(h,L,s,u,w,dtau);
    h_his(iter) = h(s(:,end),w);
    s_his(:,:,iter) = s;
    u_his(:,:,iter) = u;
    lambda_his(:,:,iter) = lambda;
    w_his(iter) = w;
    Tmax_his(:,:,iter) = ub(2,:);
    
    % termination condition
    fprintf('iter: %d, J: %.4f, h: %.4f, J-h: %.4f, deltau: %.4f, Ku: %f, deltaw: %.4f, Kw: %f, w: %.2f\n', ...
        iter, J_his(iter), h_his(iter), J_his(iter)-h_his(iter), deltau, Ku, ...
        deltaw, Kw, w_his(iter));
    if deltau < err_thresh || iter == maxIter
        break;
    end
    
    % adjust Ku, Kw
    if iter > 1
        if J_his(iter) < J_his(iter-1)
            Ku = min(Ku_max, Ku*1.01);
            Kw = min(Kw_max, Kw*1.01);
        else
            Ku = max(Ku_min, Ku/1.5);
            Kw = max(Kw_min, Kw/1.5);
        end
    end
    
    % input update
    ustep = zeros(size(u));
    lb = zeros(size(u)); ub = zeros(size(u));
    dHdws = zeros(1,length(taus));
    for i = 1:size(u,2)
        ustep(:,i) = [(pi/180)^2; (1.0*m*g)^2].*dHdu(s(:,i),u(:,i),lambda(:,i),w).';
        lb(:,i) = u_min(s(:,i));
        ub(:,i) = u_max(s(:,i));
        dHdws(i) = dHdw(s(:,i),u(:,i),lambda(:,i),w);
    end
    u_old = u;
    u = u - Ku*ustep;
    u = min(max(u,lb), ub); % satisfy control bounds
    deltau = mean(vecnorm(u-u_old,2,1)); % time-averaged deltau norm
    % time scale update
    wstep = sum(dHdws*dtau) + dhdw(s(:,end),w);
    w_old = w;
    w = w - Kw/tauf*wstep;
    w = max(w, 0.1); % prevent w from becoming 0
    deltaw = mean(vecnorm(w-w_old,2,1)); % time-averaged deltaw norm
end

% figure 1. trajectory
figure(1);
plot(s_his(3,:,end)*0.001, s_his(4,:,end)*0.001, 'linewidth', 3); hold on; grid on;
plot((xtanker0 + w*Vtanker*taus)*0.001, htanker*ones(size(taus))*0.001, 'linewidth', 3);
set(gca,'ColorOrderIndex',1);
plot(s_his(3,end,end)*0.001, s_his(4,end,end)*0.001, 'x', 'linewidth', 3, 'markersize', 15);
plot((xtanker0 + w*Vtanker*tauf)*0.001, htanker*0.001, 'x', 'linewidth', 3, 'markersize', 15);
xlabel('Downrange, $x$ (1000ft)', 'interpreter', 'latex', 'fontsize', 14);
ylabel('Altitude above sea level, $h$ (1000ft)', 'interpreter', 'latex', 'fontsize', 14);
axis equal;
legend({'fighter jet', 'tanker'}, 'location', 'northeast', 'autoupdate', 'off');
saveas(gcf, 'figure1.jpg');

% figure 2. s, u vs. t (1/2)
figure(2);
[ax, h1, h2] = plotyy(w*taus, s_his(3:4,:,end)/1000, w*taus, s_his(1,:,end)); hold on; grid on;
set(h1, 'linewidth', 3); set(h2, 'linewidth', 3, 'linestyle', ':'); set(ax(2), 'ycolor', 'k');
xlabel('Time (s)', 'interpreter', 'latex', 'fontsize', 14);
ylabel(ax(1), '$x(t)$, $h(t)$ (1000ft)', 'interpreter', 'latex', 'fontsize', 14);
ylabel(ax(2), '$V(t)$ (ft/s)', 'interpreter', 'latex', 'fontsize', 14);
outerPositionDifference = get(ax(2), 'OuterPosition') - get(ax(1), 'OuterPosition');
rightEdgeShift = outerPositionDifference(1) + outerPositionDifference(3);
set(ax(1), 'OuterPosition', [0, 0, 1-rightEdgeShift, 1]);
legend({'$x(t)$', '$h(t)$', '$V(t)$'}, 'location', 'northeast', 'interpreter', 'latex', 'fontsize', 12);
yrange = get(gca, 'ylim');
line([w*tauf, w*tauf], yrange, 'linewidth', 1, 'color', 'k');
text(w*tauf*0.99, yrange(1)+(yrange(2)-yrange(1))*0.25, sprintf('t_f = %.1f', w*tauf), 'horizontalalignment', 'right');
saveas(gcf, 'figure2.jpg');

% figure 3. s, u vs. t (2/2)
figure(3);
[ax, h1, h2] = plotyy(w*taus, [s_his(2,:,end); u_his(1,:,end)]*180/pi, w*taus, [u_his(2,:,end); Tmax_his(:,:,end)]/g); hold on; grid on;
set(h1, 'linewidth', 3); set(h2(1), 'linewidth', 3, 'linestyle', ':'); set(h2(2), 'linewidth', 3, 'linestyle', '--'); set(ax(2), 'ycolor', 'k');
xlabel('Time (s)', 'interpreter', 'latex', 'fontsize', 14);
ylabel(ax(1), '$\gamma(t)$, $\alpha(t)$ (deg)', 'interpreter', 'latex', 'fontsize', 14);
ylabel(ax(2), '$T(t)$, $T_max(t)$ (lbf)', 'interpreter', 'latex', 'fontsize', 14);
outerPositionDifference = get(ax(2), 'OuterPosition') - get(ax(1), 'OuterPosition');
rightEdgeShift = outerPositionDifference(1) + outerPositionDifference(3);
set(ax(1), 'OuterPosition', [0, 0, 1-rightEdgeShift, 1]);
legend({'$\gamma(t)$', '$\alpha(t)$', '$T(t)$', '$T_{max}(t)$'}, 'location', 'northeast', 'interpreter', 'latex', 'fontsize', 12);
yrange = get(gca, 'ylim');
line([w*tauf, w*tauf], yrange, 'linewidth', 1, 'color', 'k');
text(w*tauf*0.99, yrange(1)+(yrange(2)-yrange(1))*0.25, sprintf('t_f = %.1f', w*tauf), 'horizontalalignment', 'right');
saveas(gcf, 'figure3.jpg');

% figure 4. lambda vs. t
figure(4);
y = lambda_his(:,:,end);
h1 = plot(w*taus, sign(y).*log10(1+abs(y)), 'linewidth', 3); hold on; grid on;
xlabel('Time (s)', 'interpreter', 'latex', 'fontsize', 14);
ylabel('Costate', 'interpreter', 'latex', 'fontsize', 14);
yticks = get(gca, 'ytick');
set(gca, 'yticklabel', sign(yticks).*(10.^abs(yticks)-1));
legend({'$\lambda_V(t)$', '$\lambda_\gamma(t)$', '$\lambda_x(t)$', '$\lambda_h(t)$'}, 'location', 'northeast', 'interpreter', 'latex', 'fontsize', 12, 'autoupdate', 'off');
yrange = get(gca, 'ylim');
line([w*tauf, w*tauf], yrange, 'linewidth', 1, 'color', 'k');
text(w*tauf*0.99, yrange(1)+(yrange(2)-yrange(1))*0.25, sprintf('t_f = %.1f', w*tauf), 'horizontalalignment', 'right');
saveas(gcf, 'figure4.jpg');

% figure 5. lift, drag, acc vs. t
figure(5);
lift_his = zeros(size(taus)); drag_his = zeros(size(taus));
at_his = zeros(size(taus)); an_his = zeros(size(taus)); a_his = zeros(size(taus));
for i = 1:length(taus)
    lift_his(i) = Lift(s_his(:,i,end), u_his(:,i,end));
    drag_his(i) = Drag(s_his(:,i,end), u_his(:,i,end));
    ret = dsdtau(s_his(:,i,end), u_his(:,i,end), 1);
    at_his(i) = ret(1); an_his(i) = s_his(1,i,end)*ret(2);
    a_his(i) = norm([at_his(i), an_his(i)]);
end
[ax, h1, h2] = plotyy(w*taus, [lift_his; drag_his]/g, w*taus, [at_his; an_his; a_his]/g); hold on; grid on;
set(h1, 'linewidth', 3); set(h2, 'linestyle', ':', 'linewidth', 2); set(ax(2), 'ycolor', 'k');
xlabel('Time (s)', 'interpreter', 'latex', 'fontsize', 14);
ylabel(ax(1), 'Lift, Drag (lbf)', 'interpreter', 'latex', 'fontsize', 14);
ylabel(ax(2), 'Acceleration in g', 'interpreter', 'latex', 'fontsize', 14);
outerPositionDifference = get(ax(2), 'OuterPosition') - get(ax(1), 'OuterPosition');
rightEdgeShift = outerPositionDifference(1) + outerPositionDifference(3);
set(ax(1), 'OuterPosition', [0, 0, 1-rightEdgeShift, 1]);
legend({'Lift', 'Drag', '$a_t(t)$', '$a_n(t)$', '$a(t)$'}, 'location', 'northeast', 'interpreter', 'latex', 'fontsize', 12);
yrange = get(gca, 'ylim');
line([w*tauf, w*tauf], yrange, 'linewidth', 1, 'color', 'k');
text(w*tauf*0.99, yrange(1)+(yrange(2)-yrange(1))*0.25, sprintf('t_f = %.1f', w*tauf), 'horizontalalignment', 'right');
saveas(gcf, 'figure5.jpg');

% figure 6. J, w vs. iter
figure(6);
[ax, h1, h2] = plotyy(1:iter, [J_his; h_his; J_his-h_his], 1:iter, w_his); hold on; grid on;
set(h1, 'linewidth', 3); set(h2, 'linestyle', ':', 'linewidth', 2); set(ax(2), 'ycolor', 'k');
xlabel('Iteration', 'interpreter', 'latex', 'fontsize', 14);
ylabel(ax(1), 'Cost', 'interpreter', 'latex', 'fontsize', 14);
ylabel(ax(2), 'Time scale, $\omega$', 'interpreter', 'latex', 'fontsize', 14);
outerPositionDifference = get(ax(2), 'OuterPosition') - get(ax(1), 'OuterPosition');
rightEdgeShift = outerPositionDifference(1) + outerPositionDifference(3);
set(ax(1), 'OuterPosition', [0, 0, 1-rightEdgeShift, 1]);
legend({'$J$', '$h(s(t_f),\omega)$', '$J$ - $h(s(t_f),\omega)$', '$\omega$'}, 'location', 'best', 'interpreter', 'latex', 'fontsize', 12);
saveas(gcf, 'figure6.jpg');

close all;

%% functions
function ret = perfidx(h, L, ss, us, w, dtau)
    Ls = zeros(1,size(us,2));
    for i = 1:size(us,2)
        Ls(i) = L(ss(:,i),us(:,i),w);
    end
    ret = h(ss(:,end),w) + sum(Ls*dtau);
end

function ret = statedynmcs(tau, s, taus, us, w, dsdtau)
    utau = zeros(size(us,1),1);
    for i = 1:size(us,1)
        utau(i) = interp1(taus, us(i,:), tau);
    end
    ret = dsdtau(s, utau, w);
end

function ret = costatedynmcs(tau, lambda, taus, ss, us, w, dlambdadtau)
    stau = zeros(size(ss,1),1); 
    utau = zeros(size(us,1),1);
    for i = 1:size(ss,1)
        stau(i) = interp1(taus, ss(i,:), tau);
    end
    for i = 1:size(us,1)
        utau(i) = interp1(taus, us(i,:), tau);
    end
    ret = dlambdadtau(stau, utau, lambda, w);
end

function ret = partialderiv(func, inputs, index, step)
    f0 = func(inputs{:});
    ret = zeros(length(f0), length(inputs{index}));
    for i = 1:length(inputs{index})
        inputs_perturb = inputs;
        inputs_perturb{index}(i) = inputs_perturb{index}(i) + step;
        ret(:,i) = (func(inputs_perturb{:}) - f0)/step;
    end
end