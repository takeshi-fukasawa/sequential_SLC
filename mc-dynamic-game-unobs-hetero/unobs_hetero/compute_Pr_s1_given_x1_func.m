function [Pr_s1_given_x1,Jac]=compute_Pr_s1_given_x1_func(...
    psteady_1,n_state_obs,n_state_unobs)
    %% 1. Forward Pass
    n_vars = length(psteady_1);
    
    p_vec = [psteady_1; 1 - sum(psteady_1)];
    
    is_clamped = (p_vec <= 1e-16);
    p_vec=p_vec.*(1-is_clamped)+ (1e-16).*is_clamped;
    %p_vec=max(p_vec,1e-16);%%%%%%

    % Compute row sums (marginal probabilities)
    P_mat = reshape(p_vec, n_state_obs, n_state_unobs);
    row_sums = sum(P_mat, 2); 
    
    % Conditional probabilities: Pr = P(x,s) / P(x)
    Pr_s1_given_x1 = P_mat ./ row_sums;
    
    %% 2. Derivative Setup
    % D = d(p_vec) / d(psteady_1)
    D = [eye(n_vars); -ones(1, n_vars)]; % [n_total x n_vars]

    %D(is_clamped, :) = 0;
    D=D.*(1-is_clamped);% Avoid numerical issues

    % d(row_sums) / d(p_1...p_n-1)
    % Sum derivatives of p belonging to the same row
    D_tensor = reshape(D, n_state_obs, n_state_unobs, n_vars);
    dV_dp1_reduced = squeeze(sum(D_tensor, 2)); % [n_state_obs x n_vars]
    
    %% 3. Apply Quotient Rule: (D*V - P*dV) / V^2
    V_vec = repmat(row_sums, n_state_unobs, 1);
    dV_dp1_full = repmat(dV_dp1_reduced, n_state_unobs, 1);
    
    % Element-wise quotient rule with broadcasting
    Jac = (D .* V_vec - p_vec .* dV_dp1_full) ./ (V_vec.^2);

end
