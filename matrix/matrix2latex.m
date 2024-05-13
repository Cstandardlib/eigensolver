

function latexCode = matrix2latex(matrix)
    % 获取矩阵的大小
    [rows, cols] = size(matrix);
    % 构建LaTeX矩阵环境
    latexCode = '\\begin{bmatrix}\n';
    % 逐行添加矩阵元素
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
            latexCode = [latexCode, '\n']; % 每行结束后添加一个换行符
        %end
    end
    latexCode = [latexCode, '\\end{bmatrix}\n'];
end

% 示例矩阵
%A = [1 2; 3 4];
% 调用函数并打印LaTeX代码
%latexCode = matrix2latex(A);
%disp(latexCode);