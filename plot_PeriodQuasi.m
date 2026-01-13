% time delay embedding for Koopman operator

clear; clc; close all; 

% example: Periodic system

% periodic 
t = 0:0.01:100; 
f1 = 1.0; 
y = sin(2 * pi * f1 * t) + 0.051 * (randn(1, size(t, 2))); 
% plot(y, '-o')
% set(gca, 'xlim', [0 1000]);
% grid on; 

% y = [sin(f1 * t); cos(f1 * t)] + 0.1 * randn(1, size(t, 2)); 

% quasi-periodic 
f2 = (sqrt(5.0) - 1.0) / 2.0; 
y = sin(2 * pi * f1 * t) + 0.2 * sin(2 * pi * f2 * t) + 0.051 * (randn(1, size(y, 2))); 
y = sin(f1 * t) + 0.2 * sin(f2 * t) + 0.051 * (randn(1, size(y, 2))); 
% plot(y, '-o')
% set(gca, 'xlim', [0 1000]);
% grid on; 

% R = 4.0; 
% r = 0.5; 
% f1 = 1.0; 
% f2 = (sqrt(5.0) - 1.0) / 2.0;
% x1 = (R + r * sin(f1 * t)) .* cos(f2 * t); 
% x2 = r * cos(f2 * t); 
% x3 = (R + r * sin(f1 * t)) .* sin(f2 * t); 
% y = [x1; x2; x3]; 

% example: Lorenz system
% y = load('lorenz.dat'); 
% y = y'; 

% uniform noise
% y = rand(1, 5000); %


% embedding 
yEm = y(:, 1:2000); % 5000 for quasiperiodic 3D space


% get the distance matrix for uni-variate time series
%% % distance matrix to create recurrence plots
% D = zeros(size(yEm, 2)/2, size(yEm, 2)/2); 
D = zeros(size(yEm, 2), size(yEm, 2)); 
% for i = 1:size(yEm, 2) / 2
for i = 1:size(yEm, 2)
    for j=1:size(yEm, 2)
    a = yEm(:, i); 
    b=yEm(:,j);
    % b = yEm(:, i+1:size(yEm, 2)/2+i); 
    % distAll = pdist2(a', b'); 
    distAll = abs(a-b); 
    D(i, j) = distAll; 
    end
end

% Ymean = mean(D, 2); 
% D = D - repmat(Ymean, 1, length(D)); 

Dtemp = D(:); 
sorted = sort(Dtemp, 'descend'); 
[~, index] = ismember(Dtemp, sorted); 
index = index / size(yEm, 2) / size(yEm, 2); 
Brank = reshape(index, size(D)); 

% SVD approximation or remove noise effects
[U, S, V] = svd(D, 'econ'); 

nx = size(D, 1); 
ny = size(D, 2); 

%% singular values; 
hps = figure('Position', [100 300 600 500], 'DefaultTextFontName', 'times', 'DefaultAxesFontName', 'times', 'DefaultTextFontSize', 24, 'DefaultAxesFontSize', 24);
axes('Position', [0.15 0.29 0.30 0.49]); hold on; box on; grid on; 
% subplot(1,2,1); 
plot(diag(S) / sum(diag(S)), 'ko-', 'LineWidth', 2); 
set(gca, 'yscale', 'log', 'xlim', [1 20]); 
xlabel('r'); 
ylabel('%'); 
text(-6, 1, '(a)', 'FontSize', 26); 

axes('Position', [0.57 0.29 0.30 0.49]); hold on; box on; grid on; 
% subplot(1,2,2); 
plot(cumsum(diag(S)) / sum(diag(S)), 'ko-', 'LineWidth', 2); 
set(gca, 'ylim', [0.5 1], 'yscale', 'log', 'xlim', [1 20]); 
xlabel('r'); 
ylabel('cum. sum. %'); 
text(-5, 1, '(b)', 'FontSize', 26); 

tightfig(hps); 
% print(hps, '-dpdf', 'cum_sum_r_periodic_D_RP.pdf'); 
% print(hps, '-dpdf', 'cum_sum_r_periodicEm_D_RP.pdf'); 
% print(hps, '-dpdf', 'cum_sum_r_QuasiPeriodic_D_RP.pdf'); 
% print(hps, '-dpdf', 'cum_sum_r_QuasiPeriodicEm_D_RP.pdf'); 
% print(hps, '-dpdf', 'cum_sum_r_lorenzEm_D_RP.pdf'); 
% print(hps, '-dpdf', 'cum_sum_r_noise_D_RP.pdf');


figure(2); 
subplot(2, 1, 1); hold on; box on; 
plot(50 * U(:, 1:3), 'LineWidth', 2); 
% plot(yEm, 'k-', 'LineWidth', 2); 
% plot(2 * D(1, :), 'g-', 'LineWidth', 2); 

subplot(2, 1, 2); hold on; box on; 
plot(50 * V(:, 1:3), 'LineWidth', 2); 
% plot(yEm, 'k-', 'LineWidth', 2); 

figure(3); 
plot3(V(:, 1), V(:, 2), V(:, 3), '-'); 


%% Distance plots
hps = figure('Position', [100 100 1200 420], 'DefaultTextFontName', 'times', 'DefaultAxesFontName', 'times', 'DefaultTextFontSize', 22, 'DefaultAxesFontSize', 22);
plotind = 1; 
subplot(2, 5, plotind); 
imagesc(D), axis xy; axis off; axis equal; 
colormap; 
set(gca, 'xlim', [0 length(D)], 'ylim', [0 length(D)], 'xtick', '', 'ytick', ''); 
title(['Original D']); 

plotind = 2; 
rSpace = [1 5 10 20]; 
for iR = 1:1:length(rSpace)
    r = rSpace(iR);  
    Xapprox = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)'; % approx. distance matrix
    subplot(2, 5, plotind); 
    imagesc(Xapprox), axis xy; axis off; axis equal; 
    colormap; 
    set(gca, 'xlim', [0 length(D)], 'ylim', [0 length(D)], 'xtick', '', 'ytick', ''); 
    title(['1:', num2str(r, '%d'), '^{th} mode']); 
    % text(500, -200, ['cr = ', num2str(100 * r * (nx + ny) / (nx * ny), '%2.2f'), '%']);
    plotind = plotind + 1; 
end
% tightfig(hps); 

% r = 1; 
% Xapprox = u6(:, 1:r) * s6(1:r, 1:r) * v6(:, 1:r)'; % approx. distance matrix

%% Recurrence Plots 
% hps = figure('Position', [100 100 1200 800], 'DefaultTextFontName', 'times', 'DefaultAxesFontName', 'times', 'DefaultTextFontSize', 22, 'DefaultAxesFontSize', 22);
% plotind = 1; 
subplot(2, 5, plotind); 
% Recurrence plots
ix = 0.005 * length(D(:)); 
[temp, ind] = sort(D(:), 'ascend'); 
threshold_p = temp(ix); 
RP = D; 
RP(RP <=  threshold_p) = 1; 
RP(RP ~= 1) = 0;
RP = RP - diag(diag(RP)); 
[row, col] = find(RP == 1); 
plot(row, col, 'k.'); axis equal; 
set(gca, 'xlim', [0 length(D)], 'ylim', [0 length(D)], 'xtick', '', 'ytick', ''); 
title(['Original RP']); 

% create network and statistics
Adj = RP - diag(diag(RP)); 
G = graph(Adj); 
degList = G.degree; 
return_time = []; 
for iNode = 1:length(degList)
    neibors = neighbors(G, iNode); 
    neibors = diff(sort([neibors; iNode])); 
    return_time = [return_time; neibors]; 
end
pU = unique(return_time); % for the plotting use only 
for ipU=1:1:size(pU, 1)
    B = pU(ipU, :); 
    indx = find(return_time == B); 
    rtnt(ipU) = length(indx) / length(return_time); 
end
% 
% % figure(1); 
% plot(pU, rtnt, 'o-'); 
% set(gca, 'yscale', 'log', 'xlim', [0 1500]); 
% 
% Y = fft(yEm); L = length(Y); 
% P2 = abs(Y / L); % 双边频谱的幅值
% P1 = P2(1:L/2+1); % 单边频谱的幅值（正频率部分）
% P1(2:end-1) = 2 * P1(2:end-1); % 由于对称性，其它频率的功率需要乘以2
% 
% figure(2); 
% plot(P1, '-o'); 
% set(gca, 'yscale', 'log', 'xlim', [0 80]); 
mrt = mean(return_time); % this is the same as sum(rtnt * pU); 
acc = avgClusteringCoefficient(Adj);
cpl = characteristicPathLength(Adj);
disp(['mrt = ', num2str(mrt, '%g'), ', acc = ', num2str(mrt, '%g'), ', cpl = ', num2str(mrt, '%g')]); 



plotind = plotind + 1; 

% plotind = 2; 
rSpace = [1 5 10 20]; 
for iR = 1:1:length(rSpace)
    r = rSpace(iR); 
    Xapprox = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)'; % approx. distance matrix
    [temp, ind] = sort(Xapprox(:), 'ascend'); 
    threshold_p = temp(ix); 
    RP = Xapprox; 
    RP(RP <=  threshold_p) = 1; 
    RP(RP ~= 1) = 0; 
    [row, col] = find(RP == 1); 
    subplot(2, 5, plotind); 
    plot(row, col, 'k.'); axis equal; 
    set(gca, 'xlim', [0 length(D)], 'ylim', [0 length(D)], 'xtick', '', 'ytick', ''); 
    title(['1:', num2str(r, '%d'), '^{th} mode']); 
    
    Adj = RP - diag(diag(RP)); 
    Adj_sym = tril(Adj) + triu(Adj', 1); % get a symmetric adjacency matrix 
    G = graph(Adj_sym); 
    degList = G.degree; 
    return_time = []; 
    for iNode = 1:length(degList)
        neibors = neighbors(G, iNode); 
        neibors = diff(sort([neibors; iNode])); 
        return_time = [return_time; neibors]; 
    end
    mrt = mean(return_time);
    acc = avgClusteringCoefficient(Adj_sym);
    cpl = characteristicPathLength(Adj_sym);
    disp(['iR = ', num2str(iR, '%d'), ', mrt = ', num2str(mrt, '%g'), ', acc = ', num2str(mrt, '%g'), ', cpl = ', num2str(mrt, '%g')]); 
    
    plotind = plotind + 1; 
end

tightfig(hps); 
% print(hps, '-dpdf', 'periodic_D_RP.pdf'); 
% print(hps, '-dpdf', 'periodicEm_D_RP.pdf'); 
% print(hps, '-dpdf', 'quasi_periodic_D_RP.pdf'); 
% print(hps, '-dpdf', 'quasi_periodicEm_D_RP.pdf'); 
% print(hps, '-dpdf', 'lorenz_Em_D_RP.pdf'); 
% print(hps, '-dpdf', 'noise_D_RP.pdf'); 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% distance Laplacian Matrix - SVD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the row sum for each node 
D = D - diag(diag(D)); 
Tii = mean(D, 2); 
DL = diag(Tii) - D; 
% normalized distance Laplacian
% Tr = power(diag(Tii), -1/2); 
% LsymD = Tr .* DL .* Tr; 
% DL = LsymD; 
% 

% SVD approximation or remove noise effects
[U, S, V] = svd(DL, 'econ'); 
% Ymean = mean(DL, 2); 
% DS = DL - repmat(Ymean, 1, length(D)); 
% [U, S, V] = svd(DS, 'econ'); 

nx = size(DL, 1); 
ny = size(DL, 2); 

%% singular values; 
hps = figure('Position', [100 300 600 500], 'DefaultTextFontName', 'times', 'DefaultAxesFontName', 'times', 'DefaultTextFontSize', 24, 'DefaultAxesFontSize', 24);
axes('Position', [0.15 0.29 0.30 0.49]); hold on; box on; grid on; 
% subplot(1,2,1); 
plot(diag(S) / sum(diag(S)), 'ko-', 'LineWidth', 2); 
set(gca, 'yscale', 'log', 'xlim', [1 200]); 
xlabel('r'); 
ylabel('singular values'); 
text(-6, 1, '(a)', 'FontSize', 26); 

axes('Position', [0.57 0.29 0.30 0.49]); hold on; box on; grid on; 
% subplot(1,2,2); 
plot(cumsum(diag(S)) / sum(diag(S)), 'ko-', 'LineWidth', 2); 
set(gca, 'yscale', 'log', 'xlim', [1 200]); %'ylim', [0.5 1], 
xlabel('r'); 
ylabel('cum. energy %'); 
text(-5, 1, '(b)', 'FontSize', 26); 

tightfig(hps); 
% print(hps, '-dpdf', 'cum_sum_r_periodic_D_RP.pdf'); 
% print(hps, '-dpdf', 'cum_sum_r_periodicEm_D_RP.pdf'); 
% print(hps, '-dpdf', 'cum_sum_r_QuasiPeriodic_D_RP.pdf'); 
% print(hps, '-dpdf', 'cum_sum_r_QuasiPeriodicEm_D_RP.pdf'); 
% print(hps, '-dpdf', 'cum_sum_r_lorenzEm_D_RP.pdf'); 
% print(hps, '-dpdf', 'cum_sum_r_noise_D_RP.pdf');


figure(2); 
subplot(2, 1, 1); 
plot(U(:, 1:5), 'LineWidth', 2); 

subplot(2, 1, 2); 
plot(V(:, 1:5), 'LineWidth', 2); 

figure(3); 
plot3(V(:, 1), V(:, 2), V(:, 3), '-'); 


%% Distance plots
hps = figure('Position', [100 100 1200 420], 'DefaultTextFontName', 'times', 'DefaultAxesFontName', 'times', 'DefaultTextFontSize', 22, 'DefaultAxesFontSize', 22);
plotind = 1; 
subplot(2, 5, plotind); 
imagesc(DL), axis xy; axis off; axis equal; 
colormap; 
set(gca, 'xlim', [0 length(yEm)], 'ylim', [0 length(yEm)], 'xtick', '', 'ytick', ''); 
title(['Original Lp']); 

plotind = 2; 
rSpace = [1 50 100 200]; 
for iR = 1:1:length(rSpace)
    r = rSpace(iR);  
    Xapprox = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)'; % approx. distance matrix
    subplot(2, 5, plotind); 
    imagesc(Xapprox), axis xy; axis off; axis equal; 
    colormap; 
    set(gca, 'xlim', [0 length(yEm)], 'ylim', [0 length(yEm)], 'xtick', '', 'ytick', ''); 
    title(['1:', num2str(r, '%d'), '^{th} mode']); 
    plotind = plotind + 1; 
end
% tightfig(hps); 

% r = 1; 
% Xapprox = u6(:, 1:r) * s6(1:r, 1:r) * v6(:, 1:r)'; % approx. distance matrix

%% Recurrence Plots 
% hps = figure('Position', [100 100 1200 800], 'DefaultTextFontName', 'times', 'DefaultAxesFontName', 'times', 'DefaultTextFontSize', 22, 'DefaultAxesFontSize', 22);
% plotind = 1; 
subplot(2, 5, plotind); 
% Recurrence plots
ix = 0.005 * length(DL(:)); 
[temp, ind] = sort(DL(:), 'ascend'); 
threshold_p = temp(ix); 
RP = DL; 
RP(RP <=  threshold_p) = 1; 
RP(RP ~= 1) = 0; 
[row, col] = find(RP == 1); 
plot(row, col, 'k.'); axis equal; 
set(gca, 'xlim', [0 length(yEm)], 'ylim', [0 length(yEm)], 'xtick', '', 'ytick', ''); 
title(['Original Lp']); 
% acc = avgClusteringCoefficient(RP)

plotind = plotind + 1; 

% plotind = 2; 
rSpace = [1 50 100 200]; 
for iR = 1:1:length(rSpace)
    r = rSpace(iR); 
    Xapprox = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)'; % approx. distance matrix
    [temp, ind] = sort(Xapprox(:), 'ascend'); 
    threshold_p = temp(ix); 
    RP = Xapprox; 
    RP(RP <=  threshold_p) = 1; 
    RP(RP ~= 1) = 0; 
    [row, col] = find(RP == 1); 
    subplot(2, 5, plotind); 
    plot(row, col, 'k.'); axis equal; 
    set(gca, 'xlim', [0 length(yEm)], 'ylim', [0 length(yEm)], 'xtick', '', 'ytick', ''); 
    title(['1:', num2str(r, '%d'), '^{th} mode']); 
    plotind = plotind + 1; 
end

tightfig(hps); 
