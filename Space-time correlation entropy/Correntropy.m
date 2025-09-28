% 假设X是您的多通道EEG数据矩阵[通道 x 时间点]
% T是时间片段的数量
% fs是采样频率
% sigma1是高斯核的宽度参数

% 步骤1：将X分解为连续的时间矩阵
X=EEG.data;
T = 5; % 片段长度（根据预处理步骤，以秒为单位）
J=size(X,1);
fs = 500; % 采样率
N = size(X, 2) / (T * fs); % 片段数量
sigma1 = selectSigmaForMaxCorrentropy(X, 0.1, 1, 100);
segment_length = T * fs;
Xn = reshape(X, size(X, 1), segment_length, N); % 将X重塑为[通道 x 片段长度 x N]

% 步骤2：计算每个通道的时间相关性
Cn = zeros(size(Xn, 1), N); % 初始化相关性矩阵
for j = 1:size(Xn, 1) % 对每个通道
    for n = 1:N % 对每个片段
        for k = 1:N % 与每个片段进行比较
            if k ~= n
                % 计算片段间的相关性
                Cn(j, n) = Cn(j, n) + exp(-sum((Xn(j, :, n) - Xn(j, :, k)).^2) / (2 * sigma1^2));
            end
        end
    end
    Cn(j, :) = Cn(j, :) / (T * fs * (N-1)); % 标准化，除以比较次数
end

% 步骤3：计算通道间的空间相关性
Z = zeros(size(Xn, 1), 1); % 初始化Z向量
for j = 1:size(Xn, 1) % 对每个通道
    for k = 1:size(Xn, 1) % 与每个通道进行比较
        if k ~= j
            % 计算空间相关性
            Z(j) = Z(j) + sum(exp(-1/(2*sigma1^2) * (Cn(j, :) - Cn(k, :)).^2));
        end
    end
end

% 步骤4：合并形成代表时空相关性的矩阵Z
Z = Z / (J * (J - 1)); % 标准化，除以通道比较次数

% 假设X是您的多通道脑电图数据矩阵
% Cn_j代表单个通道的相关性
% 初始化数据（这应该替换为真实的EEG数据）
X = randn(32, 1000); % 32个通道，每个通道1000个时间点
Cn_j = randn(32, 1); % 每个通道的相关性

% 计算每个通道的神经调制强度
P = zeros(size(X, 1), 1); % 向量，用于保存每个通道的神经调制强度

for j = 1:size(X, 1) % 对于每个通道
    % 找到偏度最小的信号片段
    L = X(j, :);
    [~, p_star] = min(skewness(L));
    Lj = L(p_star);
    
    % 指示神经调制活动 Sj
    Sj = X(j, Cn_j < prctile(Cn_j, p_star));
    
    % 计算 S 在所有信号片段中的比例
    Pj = abs(Sj) / (abs(Sj) + abs(Lj));
    P(j) = Pj; % 存储通道 j 的神经调制强度
end

% 计算所有通道的神经调制强度
P_global = norm(P, 1); % P的L1范数
