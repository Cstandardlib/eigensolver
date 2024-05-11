import re
import matplotlib.pyplot as plt

# 文件名列表
filenames = ['out/gen1000thin_1_3_nopre.txt', 'out/gen1000thin_1_3_pre.txt']
# filenames = ['out/gen1000thin_1_3_nopre.txt', 'out/gen1000thin_1_3_pre.txt']

# 初始化变量，用于存储所有文件的数据
all_first_num_lists = []
all_fourth_num_lists = []

# 循环读取每个文件
for filename in filenames:
    first_num_list = []
    fourth_num_list = []

    # 读取文件并解析每一行
    with open(filename, 'r') as file:
        for line in file:
            # 使用正则表达式匹配一行中的五个数值
            match = re.match(r'([-+]?\d+)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)', line)
            if match:
                first_num_list.append(float(match.group(1)))  # 提取第一个数
                fourth_num_list.append(float(match.group(4)))  # 提取第四个数

    # 限制数据点的数量，如果需要
    first_num_list = first_num_list[0:30]
    fourth_num_list = fourth_num_list[0:30]

    # 存储每个文件的数据
    all_first_num_lists.append(first_num_list)
    all_fourth_num_lists.append(fourth_num_list)

    # 打印提取的数值以验证
    print(f"Data from {filename}:")
    print("First number list:", first_num_list)
    print("Fourth number list:", fourth_num_list)

# 绘图
plt.figure(figsize=(10, 5))

# 循环遍历每个文件的数据并绘制
for first_nums, fourth_nums in zip(all_first_num_lists, all_fourth_num_lists):
    plt.plot(first_nums, fourth_nums, marker='o', label=f'Data from {filenames[all_first_num_lists.index(first_nums)]}')

plt.title('LOBPCG Residuals vs Iterations')
plt.xlabel('Iteration')
plt.ylabel('Residuals')
plt.grid(True)
plt.legend()
plt.show()