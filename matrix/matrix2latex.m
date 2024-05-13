

function latexCode = matrix2latex(matrix)
    % ��ȡ����Ĵ�С
    [rows, cols] = size(matrix);
    % ����LaTeX���󻷾�
    latexCode = '\\begin{bmatrix}\n';
    % ������Ӿ���Ԫ��
    for i = 1:rows
        for j = 1:cols
            latexCode = [latexCode, num2str(matrix(i,j), '%.4g')];
            if j ~= cols
                latexCode = [latexCode, ' & '];
            end
        end
        %if i ~= rows
            latexCode = [latexCode, ' \\\\'];
        %end
        %if i ~= rows
            latexCode = [latexCode, '\n']; % ÿ�н��������һ�����з�
        %end
    end
    latexCode = [latexCode, '\\end{bmatrix}\n'];
end

% ʾ������
%A = [1 2; 3 4];
% ���ú�������ӡLaTeX����
%latexCode = matrix2latex(A);
%disp(latexCode);