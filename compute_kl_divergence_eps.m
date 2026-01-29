function kl_div = compute_kl_divergence_eps(P, Q, epsilon)
% 计算两个离散概率分布的KL散度，添加epsilon避免无穷大
% 输入：
%   P, Q - 两个概率分布（行向量或列向量）
%   epsilon - (可选) 平滑参数，默认值为1e-10
% 输出：
%   kl_div - KL散度值（以2为底的对数，单位bits）

% 设置默认epsilon值
if nargin < 3
    epsilon = 1e-10;
end

% 转换为列向量并确保是数值类型
P = P(:);
Q = Q(:);

% 检查输入是否有效
if length(P) ~= length(Q)
    error('两个分布的长度必须相同');
end

if any(P < 0) || any(Q < 0)
    error('概率分布不能包含负数');
end

% 归一化处理
P_sum = sum(P);
Q_sum = sum(Q);
P_norm = P / P_sum;
Q_norm = Q / Q_sum;

% 添加epsilon避免0值问题
% 方法1：直接添加epsilon（简单有效）
P_smooth = P_norm + epsilon;
Q_smooth = Q_norm + epsilon;

% 重新归一化
P_smooth = P_smooth / sum(P_smooth);
Q_smooth = Q_smooth / sum(Q_smooth);

% 方法2：仅对0元素添加epsilon（更精确但复杂）
% 也可以使用下面的代码代替上面的简单平滑
% P_smooth = P_norm;
% Q_smooth = Q_norm;
% zero_mask_Q = Q_smooth == 0;
% if any(zero_mask_Q)
%     Q_smooth(zero_mask_Q) = epsilon;
%     Q_smooth = Q_smooth / sum(Q_smooth);
% end
% % 同样处理P，确保一致性
% zero_mask_P = P_smooth == 0;
% if any(zero_mask_P)
%     P_smooth(zero_mask_P) = epsilon;
%     P_smooth = P_smooth / sum(P_smooth);
% end

% 计算KL散度（使用自然对数，最后转换为以2为底）
% 使用log2可以直接得到以2为底的结果
kl_div = sum(P_smooth .* log2(P_smooth ./ Q_smooth));

% 可选：输出一些检查信息
fprintf('输入分布检查：\n');
fprintf('  P是否归一化：%.6f (原始) -> %.6f (归一化后) -> %.6f (平滑后)\n', ...
    P_sum, sum(P_norm), sum(P_smooth));
fprintf('  Q是否归一化：%.6f (原始) -> %.6f (归一化后) -> %.6f (平滑后)\n', ...
    Q_sum, sum(Q_norm), sum(Q_smooth));
fprintf('  使用的epsilon值：%g\n', epsilon);
fprintf('KL散度计算结果：%.6f bits\n', kl_div);

% 检查原始分布中0值的情况
P_zero_count = sum(P_norm == 0);
Q_zero_count = sum(Q_norm == 0);
if P_zero_count > 0 || Q_zero_count > 0
    fprintf('注意：原始分布中有0值元素（P:%d个, Q:%d个），已应用平滑处理\n', ...
        P_zero_count, Q_zero_count);
end
end