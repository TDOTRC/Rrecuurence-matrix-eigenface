function [MATCHINGINDEX,MEANDEGREE,TRANSIVITY,GLOBALEFFICIENCY,MEANCLUSTER]= CN(RN)
%% mean degree
N = size(RN, 1);
degree=sum(RN,1);
degreedist=zeros(1,N);
for i=1:N
    if degree(i)<N
        degreedist(degree(i)+1)=degreedist(degree(i)+1)+1;
    else
        degreedist(degree(i))=degreedist(degree(i))+1;
    end
end
MEANDEGREE=([0:N-1]*degreedist')/N;

% C=zeros(N,1);
% for u=1:N
%     V=find(RN(u,:));
%     k=length(V);
%     if k>=2                 
%         S=RN(V,V);
%         C(u)=sum(S(:))/(k^2-k);
%     end
% end
% MEANCLUSTER=mean(C);

%TRANSITIVITY_BU    Transitivity
%
%   T = transitivity_bu(A);
%
%   Transitivity is the ratio of 'triangles to triplets' in the network.
%   (A classical version of the clustering coefficient).
%
%   Input:      A       binary undirected connection matrix
%
%   Output:     T       transitivity scalar
%
%   Reference: e.g. Humphries et al. (2008) Plos ONE 3: e0002051
%
%
%   Alexandros Goulas, Maastricht University, 2010
%% Transivity
TRANSIVITY = trace(RN^3) / (sum(sum(RN^2)) - trace(RN^2));
%% Mean cluster
C=zeros(N,1);
for u=1:N
    V=find(RN(u,:));
    k=length(V);
    if k>=2                 
        S=RN(V,V);
        C(u)=sum(S(:))/(k^2-k);
    end
end
MEANCLUSTER=mean(C);
%% global efficiency
% 创建无向图对象
    G = graph(RN);
    
    % 计算所有节点对之间的最短路径距离矩阵
    d = distances(G);
    
    % 获取矩阵的大小
    N = size(RN, 1);
    
    % 初始化距离倒数矩阵为全0
    inv_d = zeros(N);
    
    % 找到满足以下条件的索引：i≠j 且 路径存在（即距离不为无穷大）
    valid_indices = (d ~= Inf) & (~eye(N));
    
    % 计算有效距离的倒数
    inv_d(valid_indices) = 1 ./ d(valid_indices);
    
    % 计算所有点对的倒数总和（排除自环）
    total = sum(inv_d(:));
    
    % 计算平均值，总共有N*(N-1)个有序点对
    GLOBALEFFICIENCY = total / (N * (N - 1));
%% matching index    
    % 计算节点度数向量
    degrees = sum(RN, 2);
    
    % 计算共同邻居矩阵
    common_neighbors = RN * RN;
    
    % 计算Jaccard分母矩阵
    denominator = degrees + degrees' - common_neighbors;
    
    % 初始化Jaccard矩阵
    jaccard = zeros(size(RN));
    
    % 处理正常情况（分母非零）
    valid_mask = (denominator ~= 0);
    jaccard(valid_mask) = common_neighbors(valid_mask) ./ denominator(valid_mask);
    
    % 处理特殊边界情况：两个孤立节点（分母和分子均为零）
    isolated_pairs = (denominator == 0) & (common_neighbors == 0);
    jaccard(isolated_pairs) = 1;
    
    % 提取所有无序点对（i<j）
    N = size(RN, 1);
    upper_triangle_mask = triu(true(N), 1);
    pairwise_jaccard = jaccard(upper_triangle_mask);
    
    % 计算平均值
    MATCHINGINDEX = mean(pairwise_jaccard(:));


[MATCHINGINDEX,MEANDEGREE,TRANSIVITY,GLOBALEFFICIENCY,MEANCLUSTER];
end
