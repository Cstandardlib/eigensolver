#include <cassert>
#include <iostream>
#include "ortho.h"
#include "lobpcg.h" // check_init_guess

// test for ortho and qr
void test_ortho(){
    std::cout << "\n----- Testing ortho -----" << std::endl;
    int n = 5; // number of rows
    int m = 3;//4; // number of columns
    Eigen::MatrixXd evec(n, m); // initialize evec with zeros
    // evec << 1, 2, 3, 4, 
    //         5, 6, 7, 8, 
    //         9, 10, 11, 12, 
    //         13, 14, 15, 16, 
    //         17, 18, 19, 20;
    evec << 1, 2, 3,
            4, 5, 6, 
            7, 8, 9, 
            10, 11, 12, 
            13, 14, 15;
    evec.setRandom(); evec *= 10;
    auto another_vec = evec;

    // std::cout << evec.col(2) << std::endl;
    // evec.col(2).normalize();
    // std::cout << evec.col(2) << std::endl;

    std::cout << "evec = \n" << evec << std::endl;
    ortho(n, m, evec);
    std::cout << "ortho evec = \n" << evec << std::endl;
    /* below will not work unless ortho is redefined as template of
    taking parameter of (int n, int m, Eigen::DenseBase<Derived> &u)*/
    // std::cout <<"---"<< std::endl;
    // std::cout << "another vec = \n" << another_vec << std::endl;
    // // ortho(n, 3, another_vec.leftCols(3));
    // std::cout << "ortho evec = \n" << another_vec << std::endl;
    auto overlap = evec.transpose() * evec;
    assert(overlap.isApprox(Eigen::MatrixXd::Identity(m, m), 1e-10) && "u'u should be approximately identity matrix");
}


// test for ortho and check_init_guess
void test_ortho_check_init_guess(){
    std::cout << "\n----- Testing ortho and check_init_guess -----" << std::endl;
    int n = 4; // number of rows
    int m = 3; // number of columns
    Eigen::MatrixXd evec(n, m); // initialize evec with zeros
    evec << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;

    std::cout << evec << std::endl;

    ortho(n, m, evec);
    // std::cout << evec << std::endl;
    // std::cout << "--- ortho end ---\n";

    evec << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
    // std::cout << evec << std::endl;

    // check_init_guess(n, m, evec);
    // std::cout << evec << std::endl;

    // Check if evec is orthogonal
    Eigen::MatrixXd overlap = evec.transpose() * evec;
    std::cout << "overlap = \n" << overlap << std::endl;
    // std::cout << "areColumnsOrthogonal(u)" << Eigen::areColumnsOrthogonal(u) << std::endl;
    double diag_norm = overlap.diagonal().array().square().sum();
    double out_norm = (overlap.array().square()).sum() - diag_norm;
    assert(std::abs(diag_norm - m) <= 1e-10);
    assert(std::abs(out_norm) <= 1e-10);
    std::cout << "----- end Testing ortho and check_init_guess -----\n" << std::endl;
}

// Test for b_ortho
void test_b_ortho() {
    std::cout << "\n----- Testing b_ortho -----" << std::endl;
    int n = 4; // number of rows
    int m = 3; // number of columns

    Eigen::MatrixXd u(n, m); // initialize u with zeros
    u << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
    ortho(n, m, u);
    std::cout << "u = \n" << u << std::endl;
    // assume b = diag(1,1,2,2)
    Eigen::MatrixXd b(n, n);
    b.setZero();
    b.diagonal() << 3, 1, 2, 4; 
    std::cout << "b = \n" << b << std::endl;
    Eigen::MatrixXd bu(n, m); // initialize bu with zeros
    bu = b*u;
    // std::cout << "bu = \n" << bu << std::endl;
    // bu << 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24;

    std::cout << "Before b_ortho:\n";
    // std::cout << "u:\n" << u << std::endl;
    // std::cout << "bu:\n" << bu << std::endl;
    std::cout << "u'u:\n" << u.transpose() * u << std::endl;
    std::cout << "u'bu:\n" << u.transpose() * bu << std::endl;

    b_ortho(n, m, u, bu);

    std::cout << "After b_ortho:\n";
    // std::cout << "u:\n" << u << std::endl;
    // std::cout << "bu:\n" << bu << std::endl;

    // TODO: Add assertions to check the correctness of the b_ortho function
    Eigen::MatrixXd overlap = u.transpose() * u;
    std::cout << "overlap u'u = \n" << overlap << std::endl;
    Eigen::MatrixXd overlap_ubu = bu.transpose() * u;
    std::cout << "overlap u'bu = \n" << overlap_ubu << std::endl;

    // assert(overlap.isApprox(Eigen::MatrixXd::Identity(m, m), 1e-10) && "u'u should be approximately identity matrix");
    assert(overlap_ubu.isApprox(Eigen::MatrixXd::Identity(m, m), 1e-10) && "u'bu should be approximately identity matrix");
    std::cout << "----- end Testing b_ortho -----\n" << std::endl;
}

// Test for ortho_against_x
void test_ortho_against_y() {
    std::cout << "\n----- Testing ortho_against_y, un-orthoganalized input -----" << std::endl;
    int n = 10; // number of rows
    int m = 3; // number of columns
    int k = 2; // number of vectors

    // Initialize x and y matrices
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(n, k);
    x << 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25;
    // ortho(n, k, x);

    Eigen::MatrixXd y = Eigen::MatrixXd::Random(n, m);
    // y << 1, 2, 3, 4, 5, 6;
    // y << 1, 2, 3, 12, 5, 7, 6, 4, 8, 9, 10, 11;
    y << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30;
    
    // ortho(n, m, y);

    std::cout << "Before ortho_against_x:\n";
    // std::cout << "x:\n" << x << std::endl;
    // std::cout << "y:\n" << y << std::endl;
    std::cout << "x'x = \n" << x.transpose()*x << std::endl;
    std::cout << "y'y = \n" << y.transpose()*y << std::endl;

    ortho_against_y(n, m, k, x, y);

    std::cout << "After ortho_against_y:\n";
    // std::cout << "x:\n" << x << std::endl;
    // std::cout << "y:\n" << y << std::endl;

    // TODO: Add assertions to check the correctness of the ortho_against_x function
    Eigen::MatrixXd overlap = x.transpose() * x;
    std::cout << "overlap x'x = \n" << overlap << std::endl;
    Eigen::MatrixXd overlap_yx = y.transpose() * x;
    std::cout << "overlap y'x = \n" << overlap_yx << std::endl;

    assert(overlap.isApprox(Eigen::MatrixXd::Identity(k, k), 1e-10) && "x'x should be approximately Identity matrix");
    // assert(overlap_yx.array().isMuchSmallerThan(Eigen::MatrixXd::Constant(m,k), ASSERT_EPSILON));
    // assert(overlap_yx.array().isMuchSmallerThan(ASSERT_EPSILON) && "y'x should be approximately Zero matrix");
    assert((overlap_yx.norm() < ASSERT_EPSILON) && "y'x should be approximately Zero matrix");
    std::cout << "----- end  Testing ortho_against_y, un-orthoganalized input -----\n" << std::endl;
}

void test_b_ortho_against_y(){
    std::cout << "\n----- Testing b_ortho_against_y, b-orthogonalize x against given y and by -----" << std::endl;
    int n = 10; // number of rows
    int m = 3; // number of columns
    int k = 2; // number of vectors

    // Initialize x and y matrices
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(n, k);
    // x << 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25;
    // ortho(n, k, x);

    Eigen::MatrixXd y = Eigen::MatrixXd::Random(n, m);
    // y << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30;
    
    ortho(n, m, y);

    Eigen::MatrixXd b(n, n);
    b.setZero();
    b.diagonal() << 3, 1, 2, 4, 6, 3, 5, 9, 3, 10; 
    // std::cout << "b = \n" << b << std::endl;
    Eigen::MatrixXd by(n, m);
    by = b*y;
    // std::cout << "by = \n" << by << std::endl;

    std::cout << "Before b_ortho_against_y:\n";
    // std::cout << "x:\n" << x << std::endl << "y:\n" << y << std::endl;
    Eigen::MatrixXd overlap = y.transpose() * y;
    std::cout << "x'x = \n" << x.transpose()*x << std::endl;
    std::cout << "overlap y'y = \n" << overlap << std::endl;
    Eigen::MatrixXd overlap_by = by.transpose() * y;
    std::cout << "overlap y'by = \n" << overlap_by << std::endl;

    b_ortho_against_y(n, m, k, x, y, by);

    std::cout << "After b_ortho_against_y:\n";
    // std::cout << "x:\n" << x << std::endl << "y:\n" << y << std::endl;

    // check the correctness after the b-ortho
    overlap = x.transpose() * x;
    std::cout << "overlap x'x = \n" << overlap << std::endl;
    Eigen::MatrixXd overlap_ybx = y.transpose() * b * x;
    std::cout << "overlap y^Tbx = \n" << overlap_ybx << std::endl;
    assert(overlap.isApprox(Eigen::MatrixXd::Identity(k, k), 1e-10) && "x'x should be approximately Identity matrix");
    assert(overlap_ybx.isApprox(Eigen::MatrixXd::Zero(m, k), 1e-10) && "y'bx should be approximately Zero matrix");
    std::cout << "----- end Testing b_ortho_against_y, b-orthogonalize x against given y and by -----\n" << std::endl;
}

int main() {
    test_ortho();
    // test_ortho_check_init_guess();
    // test_b_ortho();
    // test_ortho_against_y();
    // test_b_ortho_against_y();
}