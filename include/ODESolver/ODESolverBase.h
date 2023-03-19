#ifndef ODE_SOLVER_BASE_H
#define ODE_SOLVER_BASE_H

#include "DynamicModel/DynamicModelBase.h"
#include "IntegrationMethod/IntegrationMethodBase.h"

#include <vector>

namespace ODESolver
{
    using namespace DynamicModel;
    using namespace IntegrationMethod;

    class ODESolverBase
    {
    public:
        virtual ~ODESolverBase() {}
        
        virtual std::vector<Eigen::VectorXd> solve(const DynamicModelBase &model, const IntegrationMethodBase &method, const Eigen::VectorXd &state, const std::vector<Eigen::VectorXd> &controls, const std::vector<double> &times) const = 0;
    };

} // namespace ODESolver

#endif