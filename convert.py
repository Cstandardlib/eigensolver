# 空白分隔的矩阵数据
matrix_data ="""
-0.045008846777  -0.012487056830   0.451774729742
-0.595457985511  -0.302271780449   0.262108946955
-0.683658739782  -0.026606629965  -0.547438361389
-0.214624268522  -0.097081891179   0.451904306374
-0.326352098345   0.603950319551   0.108625784707
-0.115209128264   0.356526801138   0.228854930509
0.075833974684   0.399936355871  -0.300634194996
0.045339444223  -0.459789809603  -0.236122012805
-0.048640002418  -0.187403125974   0.113945458005
"""
"""
-0.04501  -0.01249   0.45177
-0.59546  -0.30227   0.26211
-0.68366  -0.02661  -0.54744
-0.21462  -0.09708   0.45190
-0.32635   0.60395   0.10863
-0.11521   0.35653   0.22885
0.07583   0.39994  -0.30063
0.04534  -0.45979  -0.23612
-0.04864  -0.18740   0.11395
"""

# 按行分割数据
lines = matrix_data.strip().split('\n')

# 初始化一个空字符串用于存储转换后的Eigen格式数据
eigen_formatted_data = ""

# 遍历每一行，替换空白为逗号，并添加到新字符串
first=1
for line in lines:
    # 将空白替换为逗号，并将转换后的行添加到新字符串
    # 同时在每一行的末尾添加一个逗号
    if first:
        first=0
    else:
        eigen_formatted_data += "        "
    eigen_formatted_data += line.replace('  ', ',').replace(' -', ',') 
    eigen_formatted_data += ',' + '\n'

evec = eigen_formatted_data[:-2] + ';'
# 打印转换后的Eigen格式数据
print("evec <<", evec)