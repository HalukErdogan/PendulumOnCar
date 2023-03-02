#ifndef PENDULUM_ON_CAR_H
#define PENDULUM_ON_CAR_H

#include "DynamicModel/DynamicModelBase.h"

#include <cmath>

namespace DynamicModel
{
    class PendulumOnACar : public DynamicModelBase
    {
    private:
        double m1;  // mass of cart
        double m2;  // mass of pendulum
        double l;   // length of the pendulum
        double g;   // gravitational constant
        
    public:
        PendulumOnACar(double mass_cart = 1.0, double mass_pendulum = 1.0, double length = 1.0, double g = 9.81) : m1(mass_cart), m2(mass_pendulum), l(length), g(g) {}

        Eigen::VectorXd f(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
        {
            const double &q1 = state(0);  // x
            const double &q2 = state(1);  // theta
            const double &dq1 = state(2); // dx/dt
            const double &dq2 = state(3); // dtheta/dt

            const double &u = control(0); // force

            Eigen::VectorXd result(4);
            result << 
                dq1,
                dq2,
                (l * m2 * std::sin(q2) * std::pow(dq2, 2) + u + m2 * g * std::cos(q2) * std::sin(q2)) / (m1 + m2 * (1 - std::pow(std::cos(q2), 2))),
                -(l * m2 * std::cos(q2) * std::sin(q2) * std::pow(dq2, 2) + u * std::cos(q2) + (m1 + m2) * g * std::sin(q2)) / (l * m1 + l * m2 * (1 - std::pow(std::cos(q2), 2)));
            return result;
        }
    };

} // namespace DynamicModel

#endif