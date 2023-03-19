#include "DynamicModel/PendulumOnCar/PendulumOnCar.h"
#include "IntegrationMethod/EulerIntegration/EulerIntegration.h"
#include "IntegrationMethod/RungeKuttaIntegration/RungeKuttaIntegration.h"
#include "IntegrationMethod/RKDPIntegration/RKDPIntegration.h"
#include "ODESolver/SimpleODESolver/SimpleODESolver.h"
#include "TrajectoryOptimizer/DDP/DDP.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

int main(){
    // system model
    DynamicModel::PendulumOnCar model;
    
    // integration method
    IntegrationMethod::RungeKuttaIntegration integrator;

    // ode solver
    ODESolver::SimpleODESolver ode;

    // trajectory optimization problem
    TrajectoryOptimizer::Problem problem;

    // trajectory optimization solver
    TrajectoryOptimizer::DDP optimizer;

    // time span
    const double t0 = 0;        // initial time
    const double tf = 5;        // final time
    const double dt = 0.01;     // delta time
    const int n = (int)((tf - t0) / dt);    // number of time steps
    std::vector<double> times(n);
    times.at(0) = t0;
    for(int i=0; i<n-1; ++i){
        times.at(i+1) = times.at(i) + dt;
    }
    
    // Set the problem
    Eigen::VectorXd u0(1);
    u0 << 0;
    problem.x0 = Eigen::VectorXd(4,1);
    problem.xt = Eigen::VectorXd(4,1);
    problem.R = Eigen::MatrixXd::Zero(1,1);
    problem.Q = Eigen::MatrixXd::Zero(4,4);
    problem.Qf = Eigen::MatrixXd::Zero(4,4);
    
    problem.x0 << 0.0, 0.0, 0.0, 0.0;
    problem.xt << 0.0, M_PI, 0.0, 0.0;
    problem.R << 0.01;
    // problem.Q.diagonal() << 1, 1, 1, 1;
    problem.Qf.diagonal() << 1, 1, 1, 1;
    problem.T = times;
    problem.X0 = std::vector<Eigen::VectorXd>(n, problem.x0);
    problem.U0 = std::vector<Eigen::VectorXd>(n, u0);

    auto res = optimizer.solve(problem, model, integrator, ode);

    // solve the ode
    auto result = ode.solve(model, integrator, problem.x0, res.U, problem.T);

    // Create a file output stream object
    std::ofstream file;
    file.open("output.csv");
    file << "t,q1,q2,dq1,dq2" << std::endl;

    // print the results
    for(int i=0; i<result.size(); ++i){
	    // Write the current state to the file
        const auto & x = result.at(i);
        file << problem.T.at(i) << "," << x(0) << "," << x(1) << "," << x(2) << "," <<  x(3) << std::endl;
    }
    file.close();
    return 0;
}
