clear;
load C.mat;
% [L, S, gamma] = CorrentropyBasedDecomposition(X, sigma, M)
% function [L, S, gamma] = CorrentropyBasedDecomposition(X, sigma, M)
    X=C_N;
%     sigma=sigma1;
%     [N, W] = size(X); % 假设 X 是 N x W 的矩阵
%     Z = zeros(N, 1); % 初始化累积相关性向量

%     % 计算累积相关性
%     for i = 1:N
%         for j = 1:W
%             for w = 1:W
%                 Z(i) = Z(i) + 1/W * exp(-(X(i, j) - X(i, w))^2 / (2 * sigma^2));
%             end
%         end
%     end

    % 识别低秩分量
    Lp_star = [];
    S_star = [];
    sp_star = inf;

    for j = 1:28 % 假设百分位数 p 从 1 到 100
        LX=[];
        for p = 1:99
         ind = find(X(j,:) >= prctile(X(j,:), p)); % 找到高于 p 百分位数的索引
         LX{p} = X(j,ind);
         pian(p) = abs(skewness(LX{p})); % 计算偏度
        end
        [~,p1] = min(pian); % 计算偏度
        L{j} = LX{p1}; % 对应的 X 分量
%         S{j} = LX(X(j,:) < prctile(X(j,:), p1));
        abs(L{j})
          S{j}=LX;
%           Pj = abs(S{j}) / (abs(S{j}) + abs(L{j}));
%        P{j}=[];
    end

        % 检查是否是最小偏度
%         if abs(sp) < sp_star
%             sp_star = abs(sp);
%             Lp_star = Lp;
%         end
%     end

%     L = X(Lp_star, :); % 最低秩分量
% 
%     % 提取稀疏分量
%     IS = find(Z < prctile(Z, Lp_star));
%     S = X(IS, :);

    % 识别和提取数据片段
%     P = findSnippets(S, M);
%     gamma = min(vecnorm(P, 2, 2)); % 计算非高斯度量
% 
%     % 辅助函数，找到数据片段
%     function P = findSnippets(S, M)
%         % 这里插入 findSnippets 函数的代码
%         % ...
%     end
% end