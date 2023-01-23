#ifndef EULER_INTEGRATION_H
#define EULER_INTEGRATION_H

#include "IntegrationMethod/IntegrationMethodBase.h"

namespace IntegrationMethod
{
    class EulerIntegration : public IntegrationMethodBase
    {
    public:
        Eigen::VectorXd integrate(const DynamicModelBase &model, const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time, const double &dt) const override
        {
            Eigen::VectorXd state_dot = model.f(state, control, time);
            return state + dt * state_dot;
        }
    };
} // namespace IntegrationMethod
#endif