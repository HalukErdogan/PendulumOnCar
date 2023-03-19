#ifndef TRAJECTORY_OPTIMIZER_BASE_H
#define TRAJECTORY_OPTIMIZER_BASE_H

#include "DynamicModel/DynamicModelBase.h"
#include "IntegrationMethod/IntegrationMethodBase.h"
#include "ODESolver/ODESolverBase.h"

#include <vector>

namespace TrajectoryOptimizer
{
    using namespace DynamicModel;
    using namespace IntegrationMethod;
    using namespace ODESolver;

    struct Problem
    {
        Eigen::VectorXd x0; // initial state
        Eigen::VectorXd xt; // target state
        Eigen::MatrixXd R;  // weight of controls in cost function
        Eigen::MatrixXd Q;  // weight of states except the final state in cost function
        Eigen::MatrixXd Qf; // weight of final state in cost function
        std::vector<double> T; // times
        std::vector<Eigen::VectorXd> X0;  // initial state vector
        std::vector<Eigen::VectorXd> U0;  // initial input vector
    };

    struct Solution
    {
        std::vector<Eigen::VectorXd> X;
        std::vector<Eigen::VectorXd> U;
    };

    class TrajectoryOptimizerBase
    {
    public:
        virtual ~TrajectoryOptimizerBase() {}
        
        virtual Solution solve(const Problem &problem, const DynamicModelBase &model, const IntegrationMethodBase &integrator,  const ODESolverBase &ode) = 0;
    };

} // namespace ODESolver

#endif