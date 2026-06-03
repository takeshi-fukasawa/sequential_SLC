function [Pr_s1_given_x1, Jac] = compute_Pr_s1_given_x1_func2(...
    psteady_1, n_state_obs, n_state_unobs)
    
    %% 1. Forward Pass
    n_vars = length(psteady_1);
    n_total = n_state_obs * n_state_unobs; % p_vecの長さ
    
    p_vec = [psteady_1; 1 - sum(psteady_1)];
    is_clamped = p_vec <= 1e-16;
    p_vec(is_clamped) = 1e-16;
    
    
    P_mat = reshape(p_vec, n_state_obs, n_state_unobs);
    row_sums = sum(P_mat, 2); % [n_state_obs x 1]
    
    Pr_s1_given_x1 = P_mat ./ row_sums;
    
    %% 2. Derivative Setup (Memory Efficient)
    % D = d(p_vec)/d(p_1...p_{n-1}) の構造を利用
    % Dは上部が単位行列、最下行が-1のベクトル。ただしis_clampedなら0。
    
    % d(row_sums)/d(p_vars) [n_state_obs x n_vars] の計算
    % D_tensorの和を効率的に求める
    dV_dp1_reduced = zeros(n_state_obs, n_vars);
    for j = 1:n_vars
        if ~is_clamped(j)
            % p_j が属する行(row)を特定
            r_idx = mod(j-1, n_state_obs) + 1;
            dV_dp1_reduced(r_idx, j) = 1;
        end
    end
    
    % 最後の要素 p_end = 1 - sum(p) の影響を全変数に反映
    last_idx = n_total;
    if ~is_clamped(last_idx)
        % 最後の要素が属する行を特定
        r_idx_last = mod(last_idx-1, n_state_obs) + 1;
        dV_dp1_reduced(r_idx_last, :) = dV_dp1_reduced(r_idx_last, :) - 1;
    end
    
    %% 3. Apply Quotient Rule: Jac = (D*V - p*dV) / V^2
    % 計算式を整理: Jac = (D ./ V) - (p_vec .* dV_full) ./ V^2
    
    % V_vec [n_total x 1]
    V_vec = repmat(row_sums, n_state_unobs, 1);
    
    % 先に第2項「- (p_vec * dV) / V^2」を計算 (巨大なブロードキャスト)
    dV_dp1_full = repmat(dV_dp1_reduced, n_state_unobs, 1);
    Jac = - (p_vec .* dV_dp1_full) ./ (V_vec.^2);
    
    % 第1項「D ./ V」を 0 でない場所（スパースなDの構造）だけに加算
    % 1. 対角成分付近 (Dの単位行列部分)
    for j = 1:n_vars
        if ~is_clamped(j)
            % Jac(j, j) += 1 / V(j)
            Jac(j, j) = Jac(j, j) + (1 / V_vec(j));
        end
    end
    
    % 2. 最終行 (Dの -1 部分)
    if ~is_clamped(last_idx)
        % Jac(last_idx, :) -= 1 / V(last_idx)
        Jac(last_idx, :) = Jac(last_idx, :) - (1 / V_vec(last_idx));
    end
    
    % 数値的安定性のための最終チェック
    %Jac(isinf(Jac) | isnan(Jac)) = 0;
end