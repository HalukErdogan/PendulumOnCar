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
        
        virtual Eigen::MatrixXd solve(const DynamicModelBase &model, const IntegrationMethodBase &method, const Eigen::VectorXd &state, const Eigen::MatrixXd &controls, const Eigen::VectorXd &dts) const = 0;
    };

} // namespace ODESolver

#endif