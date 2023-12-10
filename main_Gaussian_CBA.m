%% construct variables
% x = -8 + 0.08 : 0.16 : 8 - 0.08;
% y = -8 + 0.08 : 0.16 : 8 - 0.08;
x = -10 + 0.1 : 0.2 : 10 - 0.1;
y = -10 + 0.1 : 0.2 : 10 - 0.1;
% x: M, t: N, y: K
M = 100;        N = 100;            K = 100;
% s = p(y|x), K * M
s = zeros(K, M);
for k = 1 : K
    for i = 1 : M
        s(k, i) = exp(-(x(i) - y(k)) * (x(i) - y(k)) / 2) / sqrt(2 * pi);
    end
end
for i = 1 : M
    s(:, i) = s(:, i) / sum(s(:, i));
end
% px = p(x), 1 * M
px = zeros(1, M);
for i = 1 : M
    px(i) = exp(-x(i) * x(i) / 2) / sqrt(2 * pi);
end
px = px / sum(px);
% pxy = p(x, y), K * M
% pxy = s .* px;
pxy = zeros(K, M);
for i = 1 : M
   pxy(:, i) = s(:, i) * px(i); 
end
% py = p(y), 1 * K
py = ones(1, M) * pxy';

p = px';        q = py';
% save("parameters_gaussian", "pxy", "s", "p", "q");

%% theoretical bound of information curve
IXY = sum(sum(pxy .* log(pxy))) - sum(pxy * log(px')) - sum(log(py) * pxy);
% card_T = log2(N);
x = [0, 2.5];        y = [IXY, IXY];         plot(x, y, 'b--', 'HandleVisibility', 'off');                 hold on;
% x = [card_T, card_T];   y = [0, IXY];           line(x, y);             hold on;

%% theoretical information curve
x_axis = 0 : 0.1 : 2;
y_axis = -0.5 * log( (1 + exp(-2 * x_axis)) / 2 );
plot(x_axis, y_axis, 'c', 'LineWidth', 5.0);       hold on;

%% approximated information curve
% initialization for IR case
ITX_target = 0.05 : 0.05 : 2.3;
ITY_target = -0.5 * log( (1 + exp(-2 * ITX_target)) / 2 );
% initialization for RI case
% ITY_target = 0.04 : 0.04 : 0.2;
% ITX_target = -log(2 * exp(-2 * ITY_target) - 1) / 2;
lambda_target = 1 + exp(2 * ITX_target);
[~, trials] = size(lambda_target);
time_list = zeros(1, trials);
iteration_list = zeros(1, trials);
lambda_list = zeros(1, trials);
max_loops = 10;

points_cba = zeros(trials, 2);
for trial = 1 : trials
    for loop = 1 : max_loops
        tic;
        [w, r, z, d, ITY, iter, lambda_ret, ITX_ret, ITY_ret, ITX_list, ITY_list] = ...
            CBA_IB_IR(M, N, K, p, s, q, ITX_target(trial), ITY_target(trial), 10, 1e4, 1e-2, 1e-6);
%         [w, r, z, d, iter, lambda_ret, ITX_ret, ITY_ret, ITX_list, ITY_list] = ...
%             CBA_IB_RI(M, N, K, p, s, q, ITY_target(trial), 10, 1e4, 1e-2, 1e-6);
        time = toc;
        time_list(trial) = time_list(trial) + time / max_loops;
        iteration_list(trial) = iteration_list(trial) + iter / max_loops;
        lambda_list(trial) = lambda_list(trial) + lambda_ret / max_loops;
    end
%     ITX = sum(p' * (w .* log(w ./ r)));
%     ITY = -sum(p' * (w .* d)) - sum(q .* log(q));
    ITX = ITX_ret;       ITY = ITY_ret;
%     fprintf("ITX_pred = %f, ITX_ret = %f; ITY_pred = %f, ITY_ret = %f;\n  lambda_pred = %f, lambda_ret = %f; number of iteration = %f, time = %f sec\n", ...
%          ITX_target(trial), ITX, ITY_target(trial), ITY, lambda_target(trial), lambda_list(trial), iteration_list(trial), time_list(trial));
    plot(ITX, ITY, 'ro');
    points_cba(trial, 1) = ITX;     points_cba(trial, 2) = ITY;
    hold on;
end
points_abp = points_cba;
% save("parameters_gaussian.mat", 'points_abp');

% plot(x,y,'Color',[0.9,0.9,0.9],'LineWidth',5.0);
% scatter(points_ba(:,1),points_ba(:,2),'bx');
% scatter(points_cba(:,1),points_cba(:,2),'ro');
% legend('theoretical curve','BA Algorithm','ABP Algorithm');
% axis([0,2.5,0,0.35]);
% saveas(gcf,'*.png');
