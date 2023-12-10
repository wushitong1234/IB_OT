
%% Information Bottleneck: 
% Construct the Markov chain Y <--> X <--> T
% CBA for IB problem, on RI function: given I, find
%       R = min I(T, X), s.t. I(T, Y) >= I. 

%% Known variables: 
%       p: p(x), M * 1;
%       s: p(y|x), K * M;
%       q: p(y), K * 1;
%       ITY_target: I, lower bound of I(T, Y). 

%% Unknown variables:
%       w: p(t|x), M * N;
%       r: p(t), 1 * N;
%       z: p(y|t), K * N;

%% Main function CBA_IB_RI
function [w, r, z, d, outer_iter, lambda, ITX_ret, ITY_ret, ITX_list, ITY_list, lambda_list] = ...
    CBA_IB_RI(M, N, K, p, s, q, ITY_target, max_inner_it, max_outer_it, inner_eps, outer_eps)
    % 1. initialization
    w = rand(M, N);         w = w ./ sum(w, 2);         r = p' * w;
    z = ( (s .* p') * w ) ./ r;                         d = -s' * log(z);    
    Const = -sum(q .* log(q));                          
    % lambda_pred = 0.1;
    lambda = 5.0;
    I_hat = ITY_target - Const;
    ITX_pred = sum(p' * (w .* log(w ./ r)));
    ITX_list = zeros(1, max_outer_it);                  ITY_list = zeros(1, max_outer_it);
    lambda_list = zeros(1, max_outer_it);

    % 2. main iteration
    for outer_it = 1 : max_outer_it
        % 2.0. pre-process lambda by bisection method (unnecessary for RI)
        d_exp = exp(-lambda * d);

        % 2.1. update lambda by Newton's method
        for inner_it = 1 : max_inner_it
            lambda_pred = lambda;
            upper_factor = G_R(d, d_exp, p, r, I_hat);
            lower_factor = G_R_prime(d, d_exp, p, r);
            factor = upper_factor / lower_factor;
            if isnan(factor) || abs(factor) > 4
                factor = -4 * (upper_factor / abs(upper_factor));
            end
            lambda = lambda - factor;
            d_exp = exp(-lambda * d);
            if abs(lambda_pred - lambda) < inner_eps || abs(G_R(d, d_exp, p, r, I_hat)) < inner_eps
                break;
            end
        end

        % 2.2. update w, r, z, d
        w = exp(-lambda * d) .* r;
        w(w < 1e-200) = 1e-200;             w(isnan(w)) = 1e-200;           w = w ./ sum(w, 2);
        r = p' * w;                         %r(r < 1e-120) = 1e-120;
        z = ( (s .* p') * w ) ./ r;         %z(z < 1e-120) = 1e-120;
        d = -s' * log(z);                   %d(d < 1e-120) = 1e-120;

        % 2.3. stop condition
        ITX = sum(p' * (w .* log(w ./ r)));
        ITX_ret = ITX;           ITY_ret = -sum(p' * (w .* d)) + Const;
        ITX_list(outer_it) = ITX_ret;   ITY_list(outer_it) = ITY_ret;
        lambda_list(outer_it) = lambda;
        if abs(ITX_ret - ITX_pred) < outer_eps
            outer_iter = outer_it;
            return;
        end
        if isnan(ITX)
            outer_iter = outer_it - 1;
            return;
        end
        ITX_pred = ITX;
    end
    outer_iter = max_outer_it - 1;
end

function g = G_R(d, d_exp, p, r, I_hat)
%     d_exp(d_exp < 1e-120) = 1e-120;
    w_tmp = d_exp * diag(r);
    w_tmp(w_tmp < 1e-120) = 1e-120;             w_tmp(isnan(w_tmp)) = 1e-120;           
    w_tmp = w_tmp ./ sum(w_tmp, 2);
    g = sum(p' * (d .* w_tmp)) + I_hat;    
%     d_exp(d_exp < 1e-120) = 1e-120;
%     g = (p ./ (d_exp * r'))' * (d .* d_exp) * r' + I_hat;
end

function g_prime = G_R_prime(d, d_exp, p, r)
    d_exp(d_exp < 1e-200) = 1e-200;
    dr = d_exp * r';                        ddr = (d .* d_exp) * r';                  dddr = (d .* d .* d_exp) * r';
    g_prime = p' * ((ddr ./ dr) .* (ddr ./ dr) - dddr ./ dr);
end

% lambda_axis = 1 : 1 : 50;
% g_axis = zeros(1, 50);
% g_prime_axis = zeros(1, 50);
% for index = 1 : 50
%   d_exp = exp(-lambda_axis(index) * d);
%   g_axis(index) = G_R(d, d_exp, p, r, I_hat);
%   g_prime_axis(index) = G_R_prime(d, d_exp, p, r);
% end
% plot(lambda_axis, g_axis);

