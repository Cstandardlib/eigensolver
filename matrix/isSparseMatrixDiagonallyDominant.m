function isDiagonallyDominant = isSparseMatrixDiagonallyDominant(a)
    % 检查输入矩阵是否为稀疏矩阵
    if ~issparse(a)
        error('Input must be a sparse matrix.');
    end
    
    % 获取矩阵的大小
    [n, m] = size(a);
    
    % 如果矩阵不是方阵，则它不能是对角占优的
    if n ~= m
        disp('矩阵不是方阵，不能判断是否对角占优');
        return;
    end
    
    % 初始化一个逻辑数组，表示每一行是否对角占优
    isDiagonallyDominant = true(1, n);
    
    % 遍历每一行，检查是否对角占优
    for i = 1:n
        % 提取第i行的对角线元素
        diagonalElement = a(i, i);
        
        % 计算第i行非对角线元素的绝对值之和
        % 使用find来获取非零元素的索引
        [~, j] = find(a(i,:));
        offDiagonalSum = sum(abs(a(i,j)), 'double');
        
        % 打印对角线元素和非对角线元素的和
        disp(['第 ' num2str(i) ' 行的对角线元素是 ' num2str(diagonalElement) ...
             '，非对角线元素的绝对值之和是 ' num2str(offDiagonalSum)]);
        
        % 比较对角线元素和非对角线元素的和
        if offDiagonalSum > abs(diagonalElement)
            isDiagonallyDominant(i) = false;
        end
    end
    
    % 如果所有行都满足对角占优条件，则整个矩阵是对角占优的
    % 否则，矩阵不是对角占优的
    isDiagonallyDominant = all(isDiagonallyDominant);
end