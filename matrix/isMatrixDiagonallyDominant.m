function isDiagonallyDominant = isMatrixDiagonallyDominant(A)
    % 检查输入是否为方阵
    [n, m] = size(A);
    if n ~= m
        disp('矩阵不是方阵，不能判断是否对角占优');
        return;
    end
    
    % 初始化一个逻辑值，表示矩阵是否对角占优
    isDiagonallyDominant = true;
    
    % 遍历每一行，检查是否对角占优
    for i = 1:n
        % 提取第i行的对角线元素
        diagonalElement = A(i, i);
        
        % 计算第i行非对角线元素的绝对值之和
        offDiagonalSum = sum(abs(A(i,:)), 'double') - abs(diagonalElement);
        
        % 打印对角线元素和非对角线元素的和
        disp(['第 ' num2str(i) ' 行的对角线元素是 ' num2str(diagonalElement) ...
             '，非对角线元素的绝对值之和是 ' num2str(offDiagonalSum)]);
        
        % 比较对角线元素和非对角线元素的和
        if offDiagonalSum >= abs(diagonalElement)
            isDiagonallyDominant = false;
            %break; % 一旦找到不满足条件的行，就终止循环
        end
        if i>10
            break;
        end
    end
end