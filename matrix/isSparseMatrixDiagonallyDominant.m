function isDiagonallyDominant = isSparseMatrixDiagonallyDominant(a)
    % �����������Ƿ�Ϊϡ�����
    if ~issparse(a)
        error('Input must be a sparse matrix.');
    end
    
    % ��ȡ����Ĵ�С
    [n, m] = size(a);
    
    % ��������Ƿ������������ǶԽ�ռ�ŵ�
    if n ~= m
        disp('�����Ƿ��󣬲����ж��Ƿ�Խ�ռ��');
        return;
    end
    
    % ��ʼ��һ���߼����飬��ʾÿһ���Ƿ�Խ�ռ��
    isDiagonallyDominant = true(1, n);
    
    % ����ÿһ�У�����Ƿ�Խ�ռ��
    for i = 1:n
        % ��ȡ��i�еĶԽ���Ԫ��
        diagonalElement = a(i, i);
        
        % �����i�зǶԽ���Ԫ�صľ���ֵ֮��
        % ʹ��find����ȡ����Ԫ�ص�����
        [~, j] = find(a(i,:));
        offDiagonalSum = sum(abs(a(i,j)), 'double');
        
        % ��ӡ�Խ���Ԫ�غͷǶԽ���Ԫ�صĺ�
        disp(['�� ' num2str(i) ' �еĶԽ���Ԫ���� ' num2str(diagonalElement) ...
             '���ǶԽ���Ԫ�صľ���ֵ֮���� ' num2str(offDiagonalSum)]);
        
        % �Ƚ϶Խ���Ԫ�غͷǶԽ���Ԫ�صĺ�
        if offDiagonalSum > abs(diagonalElement)
            isDiagonallyDominant(i) = false;
        end
    end
    
    % ��������ж�����Խ�ռ�������������������ǶԽ�ռ�ŵ�
    % ���򣬾����ǶԽ�ռ�ŵ�
    isDiagonallyDominant = all(isDiagonallyDominant);
end