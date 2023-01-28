#ifndef RKDP_INTEGRATION_H
#define RKDP_INTEGRATION_H

#include "IntegrationMethod/IntegrationMethodBase.h"

namespace IntegrationMethod
{
    class RKDPIntegration : public IntegrationMethodBase
    {
    private:
        // Coefficients for the RKDP method
        const double a2 = 1.0 / 5.0, a3 = 3.0 / 10.0, a4 = 4.0 / 5.0, a5 = 8.0 / 9.0, a6 = 1.0;
        const double b21 = 1.0 / 5.0, b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 44.0 / 45.0, b42 = -56.0 / 15.0, b43 = 32.0 / 9.0;
        const double b51 = 19372.0 / 6561.0, b52 = -25360.0 / 2187.0, b53 = 64448.0 / 6561.0, b54 = -212.0 / 729.0;
        const double b61 = 9017.0 / 3168.0, b62 = -355.0 / 33.0, b63 = 46732.0 / 5247.0, b64 = 49.0 / 176.0, b65 = -5103.0 / 18656.0;
        const double c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 4.0 / 5.0, c5 = 8.0 / 9.0, c6 = 1.0;
        const double d1 = 35.0 / 384.0, d3 = 500.0 / 1113.0, d4 = 125.0 / 192.0, d5 = -2187.0 / 6784.0, d6 = 11.0 / 84.0;

    public:
        Eigen::VectorXd integrate(const DynamicModelBase &model, const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time, const double &dt) const override
        {
            Eigen::VectorXd k1 = model.f(state, control, time);
            Eigen::VectorXd k2 = model.f(state + dt * a2 * k1, control, time + a2 * dt);
            Eigen::VectorXd k3 = model.f(state + dt * (b31 * k1 + b32 * k2), control, time + a3 * dt);
            Eigen::VectorXd k4 = model.f(state + dt * (b41 * k1 + b42 * k2 + b43 * k3), control, time + a4 * dt);
            Eigen::VectorXd k5 = model.f(state + dt * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4), control, time + a5 * dt);
            Eigen::VectorXd k6 = model.f(state + dt * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5), control, time + a6 * dt);
            
            return state + dt * (d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6);
        }
    };
} // namespace IntegrationMethod

#endif