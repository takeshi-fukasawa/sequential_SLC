%% Verification Script: Finite Difference Test
clear; clc;

% 1. Setup Test Dimensions
n_obs = 4; 
n_unobs = 3;
n_total = n_obs * n_unobs;
n_vars = n_total - 1;

% 2. Generate Random valid psteady_1 (must sum to < 1)
% Using a Dirichlet-like approach for a stable distribution
temp_p = rand(n_total, 1);
%temp_p = temp_p / sum(temp_p);
psteady_1 = temp_p(1:end-1); 

% 3. Calculate Analytical Jacobian using your function
[Pr_analytic, Jac_analytic] = compute_Pr_s1_given_x1_func(psteady_1, n_obs, n_unobs);

% 4. Calculate Numerical Jacobian using Central Difference
epsilon = 1e-7;
Jac_numeric = zeros(n_total, n_vars);

for i = 1:n_vars
    % Perturb parameter i positively
    p_plus = psteady_1;
    p_plus(i) = p_plus(i) + epsilon;
    [Pr_plus, ~] = compute_Pr_s1_given_x1_func(p_plus, n_obs, n_unobs);
    
    % Perturb parameter i negatively
    p_minus = psteady_1;
    p_minus(i) = p_minus(i) - epsilon;
    [Pr_minus, ~] = compute_Pr_s1_given_x1_func(p_minus, n_obs, n_unobs);
    
    % Central difference: (f(x+h) - f(x-h)) / (2*h)
    % Note: we vectorize Pr to match the Jacobian rows
    diff_vec = (Pr_plus(:) - Pr_minus(:)) / (2 * epsilon);
    Jac_numeric(:, i) = diff_vec;
end

% 5. Compare Results
max_diff = max(abs(Jac_analytic(:) - Jac_numeric(:)));
mean_diff = mean(abs(Jac_analytic(:) - Jac_numeric(:)));

fprintf('--- Jacobian Verification Results ---\n');
fprintf('Max Absolute Difference: %e\n', max_diff);
fprintf('Mean Absolute Difference: %e\n', mean_diff);

if max_diff < 1e-8
    fprintf('Status: SUCCESS (The Jacobian is correct)\n');
else
    fprintf('Status: FAILURE (Difference too large)\n');
end

%% Plotting the error structure (Optional)
%figure;
%subplot(1,2,1); imagesc(Jac_analytic); title('Analytical'); colorbar;
%subplot(1,2,2); imagesc(Jac_numeric);  title('Numerical'); colorbar;