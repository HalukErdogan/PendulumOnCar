#ifndef PENDULUM_ON_CAR_H
#define PENDULUM_ON_CAR_H

#include "DynamicModel/DynamicModelBase.h"

#include <cmath>

namespace DynamicModel
{
    class PendulumOnACar : public DynamicModelBase
    {
    private:
        double m1; // mass of cart
        double m2; // mass of pendulum
        double l;  // length of the pendulum
        double g;  // gravitational constant

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
            result(0) = dq1;
            result(1) = dq2;
            result(2) = (l * m2 * std::sin(q2) * std::pow(dq2, 2) + u + m2 * g * std::cos(q2) * std::sin(q2)) / (m1 + m2 * (1 - std::pow(std::cos(q2), 2)));
            result(3) = -(l * m2 * std::cos(q2) * std::sin(q2) * std::pow(dq2, 2) + u * std::cos(q2) + (m1 + m2) * g * std::sin(q2)) / (l * m1 + l * m2 * (1 - std::pow(std::cos(q2), 2)));
            return result;
        }

        Eigen::MatrixXd df_dx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
        {
            const double &theta = state(2);
            const double &theta_dot = state(3);

            Eigen::MatrixXd result(4, 4);
            result << 0, 1, 0, 0,
                (m2 * l * cos(theta) * theta_dot * theta_dot + m2 * g * sin(theta)) / (m1 + m2 * sin(theta) * sin(theta)), 0, 0, 0,
                0, 0, 0, 1,
                (m2 * l * cos(theta) * theta_dot * theta_dot - (m1 + m2) * g) / (l * (m1 + m2 * sin(theta) * sin(theta))), 0, 0, 0;
            return result;
        }

        Eigen::MatrixXd df_du(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
        {
            const double &theta = state(2);

            Eigen::MatrixXd result(4, 1);
            result << 1 / (m1 + m2 * sin(theta) * sin(theta)),
                0,
                -cos(theta) / (l * (m1 + m2 * sin(theta) * sin(theta))),
                -sin(theta) / (l * (m1 + m2 * sin(theta) * sin(theta)));
            return result;
        }
    };

} // namespace DynamicModel

#endif