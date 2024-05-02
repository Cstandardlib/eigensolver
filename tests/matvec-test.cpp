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
    // Expected output:
    // avecs1 = [[  4.5,   6.0],
    //            [  9.0,  12.0],
    //            [ 13.5,  18.0]]
    // assert(avecs1.isApprox(Eigen::MatrixXd(3, 2) << 4.5, 6.0,
    //                                                 9.0, 12.0,
    //                                                 13.5, 18.0));

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
    // Expected output:
    // avecs2 = [[  8.5,  10.0,  11.5],
    //            [ 19.0,  22.0,  25.0],
    //            [ 29.5,  34.0,  38.5],
    //            [ 40.0,  46.0,  52.0]]
    // assert(avecs2.isApprox(Eigen::MatrixXd(4, 3) << 8.5, 10.0, 11.5,
    //                                                 19.0, 22.0, 25.0,
    //                                                 29.5, 34.0, 38.5,
    //                                                 40.0, 46.0, 52.0));

    // Test case 3: n = 4, m = 1
    m=1;
    Eigen::MatrixXd vecs3(n, m);
    vecs3 << 1, 2, 3, 4;
    std::cout << "vecs3 (" << vecs3.rows() << ", " << vecs3.cols() << ")" << std::endl;
    Eigen::MatrixXd avecs3(n, m);
    Eigen::MatrixXd avecs33 = a * vecs3;
    avec(n, m, vecs3, avecs3);
    std::cout << "avec(n, m, vecs3, avecs3) =\n" << avecs3 << std::endl;
    std::cout << "a * vecs3 =\n" << avecs33 << std::endl;
    // Expected output:
    // avecs3 = [[  3.5],
    //            [  6.0]]
    // assert(avecs3.isApprox(Eigen::MatrixXd(2, 1) << 3.5,
    //                                                 6.0));
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
    test_bvec();
}