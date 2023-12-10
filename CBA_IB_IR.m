
%% Information Bottleneck: 
% Construct the Markov chain Y <--> X <--> T
% CBA for IB problem, on IR function: given R, find
%       I = max I(T, Y), s.t. I(T, X) <= R. 

%% Known variables: 
%       p: p(x), M * 1;
%       s: p(y|x), K * M;
%       q: p(y), K * 1;
%       ITX_target: R, the upper bound of I(T, X). 

%% Unknown variables:
%       w: p(t|x), M * N;
%       r: p(t), 1 * N;
%       z: p(y|t), K * N;

%% Main function CBA_IB_IR
function [w, r, z, d, ITY, outer_it, lambda, ITX_ret, ITY_ret, ITX_list, ITY_list] = ...
                CBA_IB_IR(M, N, K, p, s, q, ITX_target, ~, max_inner_it, max_outer_it, inner_eps, outer_eps)
    % 1. initialization
    w = rand(M, N);         w = w ./ sum(w, 2);         r = p' * w;
    z = ( (s .* p') * w ) ./ r;                          
    d = -s' * log(z);    
    Const = -sum(q .* log(q));                          
    lambda_pred = 0.1;                                  lambda = 10.0;
    ITY_pred = -sum(p' * (w .* d)) + Const;
    ITX_list = zeros(1, max_outer_it);                  ITY_list = zeros(1, max_outer_it);

    % 2. main iteration
    for outer_it = 1 : max_outer_it
        % 2.0. pre-process lambda by bisection method
        lambda_left = min(lambda, lambda_pred);         lambda_right = max(lambda, lambda_pred);
        lambda_mid = (lambda_left + lambda_right) / 2;
        G_left = G_I(exp(-lambda_left * d), p, r, ITX_target);
        G_right = G_I(exp(-lambda_right * d), p, r, ITX_target);
        if outer_it <= 2
            for pre_iter = 1 : 10
                if G_left > 0
                    lambda_right = lambda_left;                     lambda_left = lambda_left / 2;
                    lambda_mid = (lambda_left + lambda_right) / 2;
                    G_right = G_left;                               G_left = G_I(exp(-lambda_left * d), p, r, ITX_target);
                elseif G_right < 0
                    lambda_left = lambda_right;                     lambda_right = lambda_right + 10.0;
                    lambda_mid = (lambda_left + lambda_right) / 2;
                    G_left = G_right;                               G_right = G_I(exp(-lambda_right * d), p, r, ITX_target);
                else
                    lambda_mid = (lambda_left + lambda_right) / 2;  G_mid = G_I(exp(-lambda_mid * d), p, r, ITX_target);
                    if G_mid < 0
                        lambda_left = lambda_mid;                   G_left = G_mid;
                    else
                        lambda_right = lambda_mid;                  G_right = G_mid;
                    end
                    if lambda_right - lambda_left < inner_eps
                        break;
                    end
                end
            end
            lambda = lambda_mid;
        end
        % lambda = 1 + exp(2 * ITX_target);

        % 2.1. update lambda by Newton's method
        for inner_it = 1 : max_inner_it
            d_exp = exp(-lambda * d);       d_exp(d_exp < 1e-200) = 1e-200;
            lambda_pred = lambda;
            % lambda = lambda - G_I(d_exp, p, r, ITX_target) / G_I_prime(d, d_exp, p, r, lambda);
            descent = G_I(d_exp, p, r, ITX_target) / G_I_prime(d, d_exp, p, r, lambda);
            lambda = lambda - descent / abs(descent) * min(10, abs(descent));
            if abs(lambda_pred - lambda) < inner_eps
                break;
            end
        end

        % 2.2. update w, r, z, d
        w = exp(-lambda * d) .* r;
        w(w < 1e-120) = 1e-120;             w(isnan(w)) = 1e-120;           w = w ./ sum(w, 2);
        r = p' * w;                         z = ( (s .* p') * w ) ./ r;     d = -s' * log(z);

        % 2.3. stop condition
        ITY_ret = -sum(p' * (w .* d)) + Const;      ITY = ITY_ret;
        ITX_ret = sum(p' * (w .* log(w ./ r)));
        ITX_list(outer_it) = ITX_ret;
        ITY_list(outer_it) = ITY_ret;
        % if abs(ITY - ITY_target) < outer_eps
        if abs(ITY_ret - ITY_pred) < outer_eps
            return;
        end
        ITY_pred = ITY_ret;
    end
end

function g = G_I(d_exp, p, r, R)
    w_tmp = d_exp .* r;
    % w_tmp(w_tmp == 0) = 1e-120;             w_tmp(isnan(w_tmp)) = 1e-120;           
    w_tmp = w_tmp ./ sum(w_tmp, 2);
    g = sum(p' * (w_tmp .* log(w_tmp))) - p' * w_tmp * log(r') - R;
end

function g_prime = G_I_prime(d, d_exp, p, r, lambda)
    dr = d_exp * r';                        ddr = (d .* d_exp) * r';                dddr = (d .* d .* d_exp) * r';    
    g_prime = lambda * p' * ( dddr ./ dr - (ddr ./ dr) .* (ddr ./ dr) );
end
