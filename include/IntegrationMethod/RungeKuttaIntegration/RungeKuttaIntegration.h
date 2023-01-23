#ifndef RUNGE_KUTTA_INTEGRATION_H
#define RUNGE_KUTTA_INTEGRATION_H

#include "IntegrationMethod/IntegrationMethodBase.h"

namespace IntegrationMethod
{
    class RungeKuttaIntegration : public IntegrationMethodBase
    {
    public:
        Eigen::VectorXd integrate(const DynamicModelBase &model, const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time, const double &dt) const override
        {
            Eigen::VectorXd k1 = model.f(state, control, time);
            Eigen::VectorXd k2 = model.f(state + dt / 2.0 * k1, control, time + dt / 2.0);
            Eigen::VectorXd k3 = model.f(state + dt / 2.0 * k2, control, time + dt / 2.0);
            Eigen::VectorXd k4 = model.f(state + dt * k3, control, time + dt);

            return state + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        }
    };
} // namespace IntegrationMethod

#endif