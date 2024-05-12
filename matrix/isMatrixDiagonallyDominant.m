function isDiagonallyDominant = isMatrixDiagonallyDominant(A)
    % ��������Ƿ�Ϊ����
    [n, m] = size(A);
    if n ~= m
        disp('�����Ƿ��󣬲����ж��Ƿ�Խ�ռ��');
        return;
    end
    
    % ��ʼ��һ���߼�ֵ����ʾ�����Ƿ�Խ�ռ��
    isDiagonallyDominant = true;
    
    % ����ÿһ�У�����Ƿ�Խ�ռ��
    for i = 1:n
        % ��ȡ��i�еĶԽ���Ԫ��
        diagonalElement = A(i, i);
        
        % �����i�зǶԽ���Ԫ�صľ���ֵ֮��
        offDiagonalSum = sum(abs(A(i,:)), 'double') - abs(diagonalElement);
        
        % ��ӡ�Խ���Ԫ�غͷǶԽ���Ԫ�صĺ�
        disp(['�� ' num2str(i) ' �еĶԽ���Ԫ���� ' num2str(diagonalElement) ...
             '���ǶԽ���Ԫ�صľ���ֵ֮���� ' num2str(offDiagonalSum)]);
        
        % �Ƚ϶Խ���Ԫ�غͷǶԽ���Ԫ�صĺ�
        if offDiagonalSum >= abs(diagonalElement)
            isDiagonallyDominant = false;
            %break; % һ���ҵ��������������У�����ֹѭ��
        end
        if i>10
            break;
        end
    end
end