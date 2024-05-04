#include "matvec.h"
#include <assert.h>
#include <iostream>
void test_avec(){
    // Test case 1: n = 4, m = 2
    int n=4;
    int m=2;
    Eigen::MatrixXd a(n, n);
    for (int i = 0; i < n; ++i) {
        a(i, i) = static_cast<double>(i + 1) + 1.0;
        for (int j = 0; j < i; ++j) {
            a(j,i) = 1.0 / static_cast<double>(i + j);
            a(i,j) = a(j,i);
        }
    }
    std::cout << "avec a = \n" << a << std::endl;
    Eigen::MatrixXd vecs1(n, m);
    vecs1 << 1, 2,
            3, 4,
            5, 6,
            7, 8;
    std::cout << "vecs1 (" << vecs1.rows() << ", " << vecs1.cols() << ")" << std::endl;
    Eigen::MatrixXd avecs1(n, m);
    Eigen::MatrixXd avecs11 = a * vecs1;
    avec(n, m, vecs1, avecs1);
    std::cout << "avec(n, m, vecs1, avecs1) =\n" << avecs1 << std::endl;
    std::cout << "a * vecs1 =\n" << avecs11 << std::endl;

    // Test case 2: n = 4, m = 3
    m = 3;
    Eigen::MatrixXd vecs2(n, m);
    vecs2 << 1, 2, 3,
            4, 5, 6,
            7, 8, 9,
            10, 11, 12;
    std::cout << "vecs2 (" << vecs2.rows() << ", " << vecs2.cols() << ")" << std::endl;
    Eigen::MatrixXd avecs2(n, m);
    Eigen::MatrixXd avecs22 = a * vecs2;
    avec(n, m, vecs2, avecs2);
    std::cout << "avec(n, m, vecs2, avecs2) =\n" << avecs2 << std::endl;
    std::cout << "a * vecs2 =\n" << avecs22 << std::endl;

    // Test case 3: n = 4, m = 1
    m=1;
    Eigen::MatrixXd vecs3(n, m); vecs3 << 1, 2, 3, 4;
    std::cout << "vecs3 (" << vecs3.rows() << ", " << vecs3.cols() << ")" << std::endl;
    Eigen::MatrixXd avecs3(n, m); avec(n, m, vecs3, avecs3);
    Eigen::MatrixXd avecs33 = a * vecs3;
    std::cout << "avec(n, m, vecs3, avecs3) =\n" << avecs3 << std::endl;
    std::cout << "a * vecs3 =\n" << avecs33 << std::endl;
}

void test_bvec(){
    std::cout << "--- test bvec where b =identity ---" << std::endl;
    int n=4;
    int m=2;

    Eigen::MatrixXd vecs1(n, m);
    vecs1 << 1, 2,
            3, 4,
            5, 6,
            7, 8;
    std::cout << "vecs1 (" << vecs1.rows() << ", " << vecs1.cols() << ")" << std::endl;
    Eigen::MatrixXd bvecs1(n, m); bvecs1.setConstant(3);
    // Eigen::MatrixXd bvecs11; bvecs11.setConstant(3);

    std::cout << "vecs1 = \n" <<vecs1 << std::endl;
    std::cout << "bvecs1 before = \n" <<bvecs1 << std::endl;

    bvec(n, m, vecs1, bvecs1);
    
    std::cout << "bvecs(n, m, vecs1, bvecs1) =\n" << bvecs1 << std::endl;
    // std::cout << "b * vecs1 =\n" << bvecs1 << std::endl;
}

int main(){
    test_avec();
    // test_bvec();
}