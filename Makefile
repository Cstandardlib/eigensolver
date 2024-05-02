# Eigen库的根目录
EIGEN_DIR = ../eigen

# 编译器
CXX = g++

# 包含Eigen头文件的路径
CXXFLAGS = -I$(EIGEN_DIR)

# 源文件目录
SRC_DIR = simple

# 源文件（根据自己的源文件列表进行修改）
SRC = $(wildcard $(SRC_DIR)/*.cpp)

# 对象文件
OBJ = $(SRC:.cpp=.o)

# 目标可执行文件
TARGET = simple_test

# 链接标志和库路径
LDFLAGS = 
# LDFLAGS = -L. -leigensolver

# 默认目标
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

.PHONY: clean

clean:
	rm -f $(OBJ) $(TARGET)

# 包含自动生成的依赖文件
-include $(OBJ:.o=.d)
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