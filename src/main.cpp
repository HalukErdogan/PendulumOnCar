#include "DynamicModel/PendulumOnCar/PendulumOnCar.h"
#include "IntegrationMethod/EulerIntegration/EulerIntegration.h"
#include "IntegrationMethod/RungeKuttaIntegration/RungeKuttaIntegration.h"
#include "IntegrationMethod/RKDPIntegration/RKDPIntegration.h"
#include "ODESolver/SimpleODESolver/SimpleODESolver.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

int main(){
    // system model
    DynamicModel::PendulumOnACar model;
    
    // integration method
    IntegrationMethod::RungeKuttaIntegration method;

    // ode solver
    ODESolver::SimpleODESolver solver;

    // time span
    const double t0 = 0;    // initial time
    const double tf = 5;  // final time
    const double dt = 0.01;  // delta time
    const int n = (int)((tf - t0) / dt);    // number of time steps
    Eigen::VectorXd times = Eigen::VectorXd::LinSpaced(n+1, t0, t0 + n * dt);

    // initial state
    Eigen::VectorXd state(4);
    state(0) = 0;       // linear positon of the car
    state(1) = M_PI_2;  // angular position of the pendulum
    state(2) = 0;       // linear velocity of the car
    state(3) = 0;       // angular velocity of the pendulum

    // control inputs
    Eigen::MatrixXd controls = Eigen::MatrixXd::Zero(n,1);
    
    // solve the ode
    auto result = solver.solve(model, method, state, controls, times);

    // Create a file output stream object
    std::ofstream file;
    file.open("output.csv");
    file << "t,q1,q2,dq1,dq2" << std::endl;

    // std::cout << result.rows() << "  " << result.cols() << std::endl;
    // print the results
    for(int i=0; i<result.rows(); ++i){
	    // Write the current state to the file
        const auto & x = result.row(i);
        file << times(i) << "," << x(0) << "," << x(1) << "," << x(2) << "," <<  x(3) << std::endl;
    }
    file.close();
    return 0;
}
