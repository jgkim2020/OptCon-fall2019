clear all; close all; clc;

% TPBVP for linear problem
% ydot(t) = A(t)y(t)
% y(0) = [alpha_1, ..., alpha_p, ...] (p initial condition)
% Cy(tf) = gamma (q terminal condition)
% y = [x1, x2, lambda1, lambda2]
% y(0) = [1, 0, unspecified, unspecified]
% y(1) = [0, 0, unspecified, unspecified]
A = [0, 1, 0, 0;
    0, 0, 0, -1;
    0, 0, 0, 0;
    0, 0, -1, 0];
C = [1, 0, 0, 0;
    0, 1, 0, 0];
y1_init = [1; 0]; gamma = [0; 0];
p = 2; q = length(A) - p;

%% i) state transition approach
PHI = expm(A);
C1 = C(:,1:p); C2 = C(:,1+p:p+q);
PHI11 = PHI(1:p,1:p); PHI12 = PHI(1:p,1+p:p+q);
PHI21 = PHI(1+p:p+q,1:p); PHI22 = PHI(1+p:p+q,1+p:p+q);
y2_init = (C1*PHI12 + C2*PHI22)\(gamma - (C1*PHI11 + C2*PHI21)*y1_init);

[t,y] = ode45(@(t,y) linsys(y,A), [0 1], [y1_init; y2_init]);
x1 = y(:,1); x2 = y(:,2); u = -y(:,4);
figure();
plot(t, x1, t, x2);
xlabel(sprintf('$t$'), 'interpreter', 'latex');
ylabel(sprintf('$state$'), 'interpreter', 'latex');
title('state transition');
h = legend('$x_1$', '$x_2$');
set(h, 'Interpreter', 'latex', 'fontSize', 12);
figure();
plot(t, u);
xlabel(sprintf('$t$'), 'interpreter', 'latex');
ylabel(sprintf('$control$'), 'interpreter', 'latex');
title('state transition');
h = legend('$u$');
set(h, 'Interpreter', 'latex', 'fontSize', 12);

%% ii) particular solution approach
M = [ones(1,q+1);
    C*PHI*[y1_init*ones(1,q+1); eye(q), zeros(q,1)]];
k = M\[1; gamma];
y_init = [y1_init*ones(1,q+1); eye(q), zeros(q,1)]*k;

[t,y] = ode45(@(t,y) linsys(y,A), [0 1], y_init);
x1 = y(:,1); x2 = y(:,2); u = -y(:,4);
figure();
plot(t, x1, t, x2);
xlabel(sprintf('$t$'), 'interpreter', 'latex');
ylabel(sprintf('$state$'), 'interpreter', 'latex');
title('particular solution');
h = legend('$x_1$', '$x_2$');
set(h, 'Interpreter', 'latex', 'fontSize', 12);
figure();
plot(t, u);
xlabel(sprintf('$t$'), 'interpreter', 'latex');
ylabel(sprintf('$control$'), 'interpreter', 'latex');
title('particular solution');
h = legend('$u$');
set(h, 'Interpreter', 'latex', 'fontSize', 12);

%% functions
function dydt = linsys(y, A)
    dydt = A*y;
end