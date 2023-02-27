#ifndef PENDULUM_ON_CAR_H
#define PENDULUM_ON_CAR_H

#include "DynamicModel/DynamicModelBase.h"

#include <cmath>

namespace DynamicModel
{
    class PendulumOnACar : public DynamicModelBase
    {
    private:
        int n = 4;  // number of system states
        int m = 1;  // number of control inputs
        double m1;  // mass of cart
        double m2;  // mass of pendulum
        double l;   // length of the pendulum
        double g;   // gravitational constant
        double eps = 1.0e-3; // epsilon for numerical differentiation

    public:
        PendulumOnACar(double mass_cart = 1.0, double mass_pendulum = 1.0, double length = 1.0, double g = 9.81) : m1(mass_cart), m2(mass_pendulum), l(length), g(g) {}

        Eigen::VectorXd f(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
        {
            const double &q1 = state(0);  // x
            const double &q2 = state(1);  // theta
            const double &dq1 = state(2); // dx/dt
            const double &dq2 = state(3); // dtheta/dt

            const double &u = control(0); // force

            Eigen::VectorXd result = Eigen::VectorXd::Zero(4, 1);
            result(0) = dq1;
            result(1) = dq2;
            result(2) = (l * m2 * std::sin(q2) * std::pow(dq2, 2) + u + m2 * g * std::cos(q2) * std::sin(q2)) / (m1 + m2 * (1 - std::pow(std::cos(q2), 2)));
            result(3) = -(l * m2 * std::cos(q2) * std::sin(q2) * std::pow(dq2, 2) + u * std::cos(q2) + (m1 + m2) * g * std::sin(q2)) / (l * m1 + l * m2 * (1 - std::pow(std::cos(q2), 2)));
            return std::move(result);
        }

        Eigen::MatrixXd fx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
        {   
            Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, n);
            Eigen::MatrixXd deltas = eps * Eigen::MatrixXd::Identity(n, n);

            // central numerical differentiation
            for(int i=0; i<n; ++i){
                result.col(i) = (f(state + deltas.col(i), control, time) - f(state - deltas.col(i), control, time)) / 2 / eps;
            }
            return std::move(result);
        }

        Eigen::MatrixXd fu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
        {
            Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, m);
            Eigen::MatrixXd deltas = eps * Eigen::MatrixXd::Identity(m, m);

            // central numerical differentiation
            for(int i=0; i<m; ++i){
                result.col(i) = (f(state + deltas.col(i), control, time) - f(state - deltas.col(i), control, time)) / 2 / eps;
            }
            return std::move(result);
        }

        std::vector<Eigen::MatrixXd> fxx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
        {
            std::vector<Eigen::MatrixXd> result(n, Eigen::MatrixXd::Zero(n, n));
            Eigen::MatrixXd deltas = eps * Eigen::MatrixXd::Identity(n, n);

            // central numerical differentiation
            for(int i = 0; i<n; ++i){
                result.at(i).col(i) = (fx(state + deltas.col(i), control, time) - fx(state - deltas.col(i), control, time)) / 2 / eps;
            }
            return std::move(result);
        }

    };

} // namespace DynamicModel

#endif