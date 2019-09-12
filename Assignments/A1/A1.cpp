#include <eigen3/Eigen/Dense>
#include <stdlib.h>
#include <iostream>

int main (int argc, char **argv)
{
    Eigen::MatrixXd m(2,2);
    m(0,0) = 1;
    m(0,1) = 2;
    m(1,0) = 3;
    m(1,1) = 4;
    std::cout << "My matrix is:\n" << m << std::endl;
    return(0);
}
