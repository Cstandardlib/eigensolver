# 指定编译器
CXX = g++

# 指定编译选项，包括Eigen库的头文件路径
CXXFLAGS = -I../eigen -std=c++11

# 指定链接选项，如果需要链接Eigen库或其他库，可以在这里添加
LDFLAGS = 

# 获取当前目录下的所有cpp文件
SRC = $(wildcard *.cpp)

# 从源文件生成对应的目标文件（.o文件）
OBJ = $(SRC:.cpp=.o)

# 从源文件生成对应的可执行文件
TARGET = $(SRC:.cpp=)

# 默认目标是所有的可执行文件
all: $(TARGET)

# 链接每个目标文件，生成对应的可执行文件
% : %.o
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@
# % : %.cpp
# 	$(CXX) $(CXXFLAGS) $<  -o $@

# 清理所有的对象文件和可执行文件
clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean

# # Eigen库的根目录
# EIGEN_DIR = ../eigen

# # 编译器
# CXX = g++

# # 包含Eigen头文件的路径
# CXXFLAGS = -I$(EIGEN_DIR)

# # 源文件目录
# SRC_DIR = ./ #simple

# # 源文件（根据自己的源文件列表进行修改）
# # SRC = $(wildcard $(SRC_DIR)/*.cpp)
# SRC = $(wildcard *.cpp)

# # 对象文件
# OBJ = $(SRC:.cpp=.o)

# # 目标可执行文件
# TARGET = simple_test

# # 链接标志和库路径
# LDFLAGS = 
# # LDFLAGS = -L. -leigensolver

# # 默认目标
# all: $(TARGET)

# $(TARGET): $(OBJ)
# 	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

# # 编译每个cpp文件
# %.o: %.cpp
# 	$(CXX) $(CXXFLAGS) -c $< -o $@

# .PHONY: clean

# clean:
# 	rm -f $(OBJ) $(TARGET)

# # 包含自动生成的依赖文件
# -include $(OBJ:.o=.d)


# CC = g++
# CFLAGS = -Wall -Wextra -I../eigen

# SRCS = $(wildcard *.cpp)
# OBJS = $(SRCS:.cpp=.o)
# EXEC = myprogram

# all: $(EXEC)

# $(EXEC): $(OBJS)
#     $(CC) $(CFLAGS) $^ -o $@

# %.o: %.cpp
#     $(CC) $(CFLAGS) -c $< -o $@

# clean:
#     rm -f $(OBJS) $(EXEC)