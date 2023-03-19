#ifndef SIMPLE_ODE_SOLVER_H
#define SIMPLE_ODE_SOLVER_H

#include "ODESolver/ODESolverBase.h"

namespace ODESolver
{
    class SimpleODESolver : public ODESolverBase
    {
    public:
        std::vector<Eigen::VectorXd> solve(const DynamicModelBase &model, const IntegrationMethodBase &method, const Eigen::VectorXd &state, const std::vector<Eigen::VectorXd> &controls, const std::vector<double> &times) const override
        {
            int n_sample = times.size();
            int n_state = state.size();
            int n_control = controls.at(0).size();

            std::vector<Eigen::VectorXd> result(n_sample, Eigen::VectorXd(n_state));
            result.at(0) = state;
            double dt = 0;
            for (int i = 0; i < controls.size() - 1; i++)
            {   
                dt = times.at(i+1) - times.at(i);
                result.at(i + 1) = method.integrate(model, result.at(i), controls.at(i), times.at(i), dt);
            }

            return result;
        }
    };
} // namespace ODESolver

#endif