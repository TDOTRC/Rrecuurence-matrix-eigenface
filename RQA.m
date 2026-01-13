
function [RR, DET, ENTR, L, DETV, ENTRV, LV] = RQA(RN, lmin1, lmin2)
% shared parameter
N = size(RN, 1);
total_ones = sum(RN, 'all');

%% ========== 斜向线优化 ==========
distributeVector = zeros(1, N);
for k = -(N-1):(N-1)
    d = diag(RN, k);
    d_extended = [0; d; 0];
    diffs = diff(d_extended);
    starts = find(diffs == 1);
    ends = find(diffs == -1);
    lengths = ends - starts;
    for m = 1:numel(lengths)
        l = lengths(m);
        if l >= 1 && l <= N
            distributeVector(l) = distributeVector(l) + 1;
        end
    end
end

%% ========== 纵向线优化 ==========
distributeVectorV = zeros(1, N);
for col = 1:N
    d = RN(:, col);
    d_extended = [0; d; 0];
    diffs = diff(d_extended);
    starts = find(diffs == 1);
    ends = find(diffs == -1);
    lengths = ends - starts;
    for m = 1:numel(lengths)
        l = lengths(m);
        if l >= 1 && l <= N
            distributeVectorV(l) = distributeVectorV(l) + 1;
        end
    end
end

%% ========== diag parameter ==========
% 公共参数
RR = total_ones / N^2;

% 斜向参数
valid_lengths = lmin1:N;
weighted_counts = valid_lengths .* distributeVector(valid_lengths);
DET = sum(weighted_counts) / total_ones;

prob_dist = distributeVector / sum(distributeVector);
cut_dist = prob_dist(valid_lengths);
ENTR = -sum(cut_dist .* log2(cut_dist + eps));

valid_lengths_L = lmin1:N;
weighted_L = valid_lengths_L .* distributeVector(valid_lengths_L);
L = sum(weighted_L) / sum(distributeVector(valid_lengths_L));

% 纵向参数
valid_lengthsV = lmin2:N; % 使用lmin1替代原代码中的lmin2参数
weighted_countsV = valid_lengthsV .* distributeVectorV(valid_lengthsV);
DETV = sum(weighted_countsV) / total_ones;

prob_distV = distributeVectorV/ sum(distributeVectorV);
cut_distV = prob_dist(valid_lengthsV);
ENTRV = -sum(cut_distV .* log2(cut_distV + eps));

valid_lengths_L = lmin2:N;
weighted_L = valid_lengths_L .* distributeVectorV(valid_lengths_L);
LV = sum(weighted_L) / sum(distributeVector(valid_lengths_L));
final_vector=[RR, DET, ENTR, L, DETV, ENTRV, LV];
end
