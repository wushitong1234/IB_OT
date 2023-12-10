%% construct variables
% x: M, t: N, y: K
M = 2;          N = 10;             K = 2;
e = 0.15;       p = 1 - e;          s = [p, e; e, p];
p = [0.5; 0.5];
pxy = 0.5 * s;
q = [0.5; 0.5];


%% theoretical bound of information curve
IXY = sum(sum(pxy .* log(pxy))) - sum(pxy * log(p)) - sum(log(q') * pxy);     card_T = log(2);
x = [0, card_T];        y = [IXY, IXY];         plot(x, y, 'b--', 'HandleVisibility', 'off');           hold on;
x = [log(2), log(2)];   y = [0, IXY];           plot(x, y, 'b--', 'HandleVisibility', 'off');           hold on;

%% theoretical information curve
u_target = 0.01 : 0.01 : 0.99;
x_axis = log(2) + u_target .* log(u_target) + (1 - u_target) .* log(1 - u_target);
x_axis = [x_axis, log(2)];          % u = 1.0
v_target = e + (1 - 2 * e) * u_target;
v_target = [v_target, 1 - e];
y_axis = log(2) + v_target .* log(v_target) + (1 - v_target) .* log(1 - v_target);
plot(x_axis, y_axis, 'b', 'LineWidth', 1.5);       hold on;

%% approximated information curve
u = 0.1 : 0.1 : 0.4;
v = e + (1 - 2 * e) * u;
ITX_target = log(2) + u .* log(u) + (1 - u) .* log(1 - u);          % R
ITY_target = log(2) + v .* log(v) + (1 - v) .* log(1 - v);          % I

ITX_target = 0 : 0.02 : 0.68;
[~, trials] = size(ITX_target);
time_list = zeros(1, trials);
iteration_list = zeros(1, trials);
lambda_list = zeros(1, trials);
max_loops = 10;

ITX_plot = zeros(trials, 1);
ITY_plot = zeros(trials, 1);

points_cba = zeros(trials, 2);
for trial = 1 : trials
    for loop = 1 : max_loops
        tic;
        [w, r, z, d, ITY, iter, lambda_ret, ITX_ret, ITY_ret, ITX_list, ITY_list] = ...
            CBA_IB_IR(M, N, K, p, s, q, ITX_target(trial), 0, 10, 1e4, 1e-2, 1e-6);
%         [w, r, z, d, iter, lambda_ret, ITX_ret, ITY_ret] = ...
%             CBA_IB_RI(M, N, K, p, s, q, ITY_target(trial), 10, 1e4, 1e-2, 1e-6);
        time = toc;
        time_list(trial) = time_list(trial) + time / max_loops;
        iteration_list(trial) = iteration_list(trial) + iter / max_loops;
        lambda_list(trial) = lambda_list(trial) + lambda_ret / max_loops;
    end
%     ITX = sum(p' * (w .* log(w ./ r)));
%     ITY = -sum(p' * (w .* d)) - sum(q .* log(q));
    ITX = ITX_ret;       ITY = ITY_ret;
    % fprintf("ITX_pred = %f, ITX_ret = %f; ITY_pred = %f, ITY_ret = %f; lambda_ret = %f; number of iteration = %f, time = %f sec\n", ...
    %      ITX_target(trial), ITX, ITY_target(trial), ITY, lambda_list(trial), iteration_list(trial), time_list(trial));
    plot(ITX_ret, ITY_ret, 'ro');
    ITX_plot(trial, 1) = ITX_ret;
    ITY_plot(trial, 1) = ITY_ret;
    hold on;
end

