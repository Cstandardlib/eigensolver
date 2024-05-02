#include <Eigen/Dense>
#include <iostream>

// 假设factYBY是YBY的Cholesky分解
// blockVectorBY, blockVectorV, blockVectorY是Eigen::MatrixXd类型的矩阵
// 这个函数将blockVectorV原地修改
void applyConstraints(
    Eigen::MatrixXd& blockVectorV,
    const Eigen::LLT<Eigen::MatrixXd>& factYBY,
    const Eigen::MatrixXd& blockVectorBY,
    const Eigen::MatrixXd& blockVectorY
) {
    // YBV = Y*BV (使用转置和共轭)
    Eigen::MatrixXd YBV = blockVectorBY.transpose() * blockVectorV;
    // Eigen::MatrixXd YBY = blockVectorY.transpose() * blockVectorBY; // 假设这是YBY的乘积
    // Eigen::LLT<Eigen::MatrixXd> factYBY(YBY);
    // tmp = x, 给定A的Cholesky分解，求解方程 YBY * x = YBV
    // std::cout << factYBY.solve(YBV) << std::endl;

    Eigen::MatrixXd tmp = factYBY.solve(YBV);

    // V = V - Y @ tmp (使用矩阵乘法更新blockVectorV)
    blockVectorV.noalias() -= blockVectorY * tmp;
}

void test_applyConstraints() {
    // 假设我们有以下矩阵
    int n = 3; // 矩阵的大小，需要根据实际情况调整
    int m = 2;
    int k = 2;
    Eigen::MatrixXd blockVectorV = Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd blockVectorBY = Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd blockVectorY = Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd YBY = blockVectorY.transpose() * blockVectorBY; // 假设这是YBY的乘积
    Eigen::MatrixXd x(n,k);
    x<<-0.846852, -0.492958,
    0.107965, -0.526819,
    -0.520755,  0.692427;
    Eigen::MatrixXd y(n,m);
    y <<-0.453376, 0.205905,
        0.599143,-0.649943,
        -0.659907,-0.731559;
    Eigen::MatrixXd yby = y.transpose() * y;
    auto by = y;

    // 对YBY进行Cholesky分解
    Eigen::LLT<Eigen::MatrixXd> factyby(yby);

    // 应用约束
    applyConstraints(x, factyby,  by, y);

    // 输出结果
    // std::cout << "Updated blockVectorV:\n" << blockVectorV << std::endl;
    std::cout << "Updated x:\n" << x << std::endl;

}



int main(){

}