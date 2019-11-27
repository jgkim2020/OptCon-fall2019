clear all; close all; clc;

epsilon = 1e-3;
maxIter = 10;
res = Inf;

V = 1; h = 1; T = 2.5;
p = 2; q = 2;
y1_init = [0; 0]; y2_term = [1; 0]; % boundary conditions

% symbolic
y = sym('y', [4 1]);
dydt = [V*y(3)/sqrt(y(3)^2 + y(4)^2) + V*y(2)/h;
    V*y(4)/sqrt(y(3)^2 + y(4)^2);
    0;
    -y(3)*V/h];
psi = [y(3) - y2_term(1);
    y(4) - y2_term(2)];
dfdy = jacobian(dydt, y);
psiy = jacobian(psi, y);
% function-ize
dydt = matlabFunction(dydt, 'vars', {y});
psi = matlabFunction(psi, 'vars', {y});
dfdy = matlabFunction(dfdy, 'vars', {y});
psiy = matlabFunction(psiy, 'vars', {y});

% true initial
y_init = [y1_init; 1; 2.5]; % true initial satisfying terminal condition

%% (a) analytical solution
t = 0:0.01:T;
x = zeros(size(t)); y = zeros(size(t)); theta = zeros(size(t));
for i = 1:length(t)
    a = sqrt(V^2*T^2 + h^2); bt = sqrt(V^2*(T-t(i))^2 + h^2);
    x(i) = V/h*((t(i)-T/2)*a + (T-t(i))/2*bt + h^2/2/V*log((bt-V*(T-t(i)))/(a-V*T)));
    y(i) = a - bt;
    theta(i) = atan2(V/h*(T-t(i)), 1);
end
plotxu(t, x, y, theta, 'analytic');

%% (b) i) state transition approach
fprintf('\ni) state transition approach \n');
y_init = [y1_init; 2; 1]; % best guess initial for nominal trajectory
nominaly_init = y_init;
for iter = 1:maxIter
    % rollout nominal y with initial guess
    [nominalt, nominaly] = ode45(@(t,y) dydt(y), [0 T], nominaly_init);
    if iter == 1
        figure();
        plot(nominalt, nominaly, '--');
        xlabel(sprintf('$t$'), 'interpreter', 'latex');
        hold on;
    end
    % find C (C1, C2) and gamma
    C = psiy(nominaly(end,:).');
    gamma = -psi(nominaly(end,:).');
    C1 = C(:,1:p); C2 = C(:,1+p:p+q);
    % find PHIT
    PHI0 = eye(p+q);
    [nominalt, PHI] = ode45(@(t,PHI) deltalindyns(t, PHI, nominalt, nominaly, dfdy), nominalt, PHI0(:));
    PHIT = reshape(PHI(end,:), size(PHI0));
    PHI11 = PHIT(1:p,1:p); PHI12 = PHIT(1:p,1+p:p+q);
    PHI21 = PHIT(1+p:p+q,1:p); PHI22 = PHIT(1+p:p+q,1+p:p+q);
    % find dely2_init
    dely1_init = zeros(p,1);
    dely2_init = (C1*PHI12 + C2*PHI22)\(gamma - (C1*PHI11 + C2*PHI21)*dely1_init);
    % rollout dely
    dely_init = [dely1_init; dely2_init];
    [nominalt, dely] = ode45(@(t,dely) deltalindyns(t, dely, nominalt, nominaly, dfdy), nominalt, dely_init);
    
    err = max(max(abs(dely)));
    fprintf('iteration %d, error %f, num intervals %d \n', iter-1, err, length(nominalt)-1);
    if err < epsilon
        break;
    else
        % update initial guess
        nominaly_init = nominaly_init + dely_init;
    end
end
ax = gca;
ax.ColorOrderIndex = 1;
plot(nominalt, nominaly);
plotxu(nominalt, nominaly(:,1), nominaly(:,2), atan2(nominaly(:,4),nominaly(:,3)), 'state transition');

%% (b) ii) particular solution approach
fprintf('\nii) particular solution approach \n');
y_init = [y1_init; 2; 1]; % best guess initial for nominal trajectory
nominaly_init = y_init;
% rollout nominal y with initial guess
[nominalt, nominaly] = ode45(@(t,y) dydt(y), [0 T], nominaly_init);
for iter = 1:maxIter
    if iter == 1
        figure();
        plot(nominalt, nominaly, '--');
        xlabel(sprintf('$t$'), 'interpreter', 'latex');
        hold on;
    end
    % find C (C1, C2) and gamma
    C = psiy(nominaly(end,:).');
    gamma = C*nominaly(end,:).' - psi(nominaly(end,:).');
    % find y^j(T) (j = 1, ..., q+1) and construct M
    M = ones(q+1,q+1);
    Ytilde = zeros(length(nominalt),p+q,q+1);
    for j = 1:q+1
        ytildej_init = zeros(p+q,1);
        ytildej_init(1:p,1) = y1_init;
        if j < q+1
            ytildej_init(p+j,1) = 1;
        end
        [nominalt, ytildej] = ode45(@(t,y) tildelindyns(t, y, nominalt, nominaly, dydt, dfdy), nominalt, ytildej_init);
        M(2:q+1,j) = C*ytildej(end,:).';
        Ytilde(:,:,j) = ytildej;
    end
    % find coefficient k
    k = M\[1; gamma];
    ytilde = zeros(size(ytildej));
    for j = 1:q+1
        ytilde = ytilde + k(j)*Ytilde(:,:,j);
    end
    
    err = max(max(abs(ytilde - nominaly)));
    fprintf('iteration %d, error %f, num intervals %d \n', iter-1, err, length(nominalt)-1);
    if err < epsilon
        break;
    else
        nominaly = ytilde;
    end
end
ax = gca;
ax.ColorOrderIndex = 1;
plot(nominalt, nominaly);
plotxu(nominalt, nominaly(:,1), nominaly(:,2), atan2(nominaly(:,4),nominaly(:,3)), 'particular solution');

%% (b) iii) method of differential correction
fprintf('\niii) differential correction \n')
dx0 = 1e-6;

y_init = [y1_init; 2; 1]; % best guess initial for nominal trajectory
nominaly_init = y_init;
for iter = 1:maxIter
    % rollout nominal y with initial guess
    [nominalt, nominaly] = ode45(@(t,y) dydt(y), [0 T], nominaly_init);
    if iter == 1
        figure();
        plot(nominalt, nominaly, '--');
        xlabel(sprintf('$t$'), 'interpreter', 'latex');
        hold on;
    end
    % find dy(T)/dy(0)
    nominalyT = nominaly(end,:).';
    M = zeros(p+q,q);
    for j = 1:q
        perturb = zeros(size(nominaly_init));
        perturb(p+j) = dx0;
        [~, nominalyj] = ode45(@(t,y) dydt(y), [0 T], nominaly_init + perturb);
        nominalyjT = nominalyj(end,:).';
        M(:,j) = (nominalyjT - nominalyT)/dx0;
    end
    % initial value update such that terminal condition holds
    dely2_init = M(1+p:p+q,:)\(y2_term - nominalyT(1+p:p+q));
    
    err = max(abs(dely2_init));
    fprintf('iteration %d, error %f, num intervals %d \n', iter-1, err, length(nominalt)-1);
    if err < epsilon
        break;
    else
        % update initial guess
        nominaly_init(1+p:p+q) = nominaly_init(1+p:p+q) + dely2_init;
    end
end
ax = gca;
ax.ColorOrderIndex = 1;
plot(nominalt, nominaly);
plotxu(nominalt, nominaly(:,1), nominaly(:,2), atan2(nominaly(:,4),nominaly(:,3)), 'differential correction');

%% functions
function ret = deltalindyns(t, x, nominalt, nominaly, dfdy)
% linearized dynamics for delta x
    dimy = size(nominaly,2);
    yt = zeros(dimy,1);
    for i = 1:dimy
       yt(i) = interp1(nominalt, nominaly(:,i), t);
    end
    At = dfdy(yt);
    if numel(x) == dimy^2
        x = reshape(x, size(At));
        ret = At*x;
        ret = ret(:);
    elseif numel(x) == dimy
        ret = At*x;
    else
        error('Invalid x dimension!');
    end
end

function ret = tildelindyns(t, y, nominalt, nominaly, dydt, dfdy)
% linearized dynamics for tilde x
    dimy = size(nominaly,2);
    yt = zeros(dimy,1);
    for i = 1:dimy
       yt(i) = interp1(nominalt, nominaly(:,i), t);
    end
    At = dfdy(yt);
    Bt = dydt(yt) - At*yt;
    ret = At*y + Bt;
end

function plotxu(t, x, y, theta, titlestr)
% plot state and control histories
    figure();
    plot(t, x, t, y);
    xlabel(sprintf('$t$'), 'interpreter', 'latex');
    ylabel(sprintf('$state$'), 'interpreter', 'latex');
    title(titlestr);
    h = legend('$x$', '$y$');
    set(h, 'Interpreter', 'latex', 'fontSize', 12);
    figure();
    plot(t, theta);
    xlabel(sprintf('$t$'), 'interpreter', 'latex');
    ylabel(sprintf('$control$'), 'interpreter', 'latex');
    title(titlestr);
    h = legend('$\theta$');
    set(h, 'Interpreter', 'latex', 'fontSize', 12);
end