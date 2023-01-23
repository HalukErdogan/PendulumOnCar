#ifndef SIMPLE_ODE_SOLVER_H
#define SIMPLE_ODE_SOLVER_H

#include "ODESolver/ODESolverBase.h"

namespace ODESolver
{
    class SimpleODESolver : public ODESolverBase
    {
    public:
        Eigen::MatrixXd solve(const DynamicModelBase &model, const IntegrationMethodBase &method, const Eigen::VectorXd &state, const Eigen::MatrixXd &controls, const Eigen::VectorXd &times) const override
        {
            Eigen::MatrixXd result(controls.size(), state.size());
            result.row(0) = state;
            double dt = 0;
            for (int i = 0; i < controls.size() - 1; i++)
            {   
                dt = times(i+1) - times(i);
                result.row(i + 1) = method.integrate(model, result.row(i), controls.row(i), times(i), dt);
            }

            return result;
        }
    };
} // namespace ODESolver

#endif