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

            Eigen::VectorXd result(4);
            result << 
                dq1,
                dq2,
                (l * m2 * std::sin(q2) * std::pow(dq2, 2) + u + m2 * g * std::cos(q2) * std::sin(q2)) / (m1 + m2 * (1 - std::pow(std::cos(q2), 2))),
                -(l * m2 * std::cos(q2) * std::sin(q2) * std::pow(dq2, 2) + u * std::cos(q2) + (m1 + m2) * g * std::sin(q2)) / (l * m1 + l * m2 * (1 - std::pow(std::cos(q2), 2)));
            return result;
        }

        Eigen::MatrixXd fx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
        {   
            Eigen::MatrixXd result(n, n);
            Eigen::VectorXd state_0(n), state_1(n);

            // central numerical differentiation
            for(int i=0; i<n; ++i){
                state_0 = state;
                state_1 = state;
                state_0(i) -= eps;
                state_1(i) += eps;
                result.col(i) = (f(state_1, control, time) - f(state_0, control, time)) / 2 / eps;
            }
            return result;
        }

        Eigen::MatrixXd fu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
        {
            Eigen::MatrixXd result(n, m);
            Eigen::VectorXd control_0(m), control_1(m);

            // first order central numerical differentiation
            for(int i=0; i<m; ++i){
                control_0 = control;
                control_1 = control;
                control_0(i) -= eps;
                control_1(i) += eps;
                result.col(i) = (f(state, control_1, time) - f(state, control_0, time)) / 2 / eps;
            }
            return result;
        }

        std::vector<Eigen::MatrixXd> fxx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
        {
            std::vector<Eigen::MatrixXd> result(n, Eigen::MatrixXd::Zero(n, n));
            Eigen::VectorXd state_0(n), state_1(n);

            // central numerical differentiation
            for(int i=0; i<n; ++i){
                state_0 = state;
                state_1 = state;
                state_0(i) -= eps;
                state_1(i) += eps;

                // dfx/dxi = [dfx1/dxi, dfx2/dxi, ... , dfxn/dxi]
                auto temp = (fx(state_1, control, time) - fx(state_0, control, time)) / 2 / eps;

                // dfx/dxi -> d2fi/dx2
                for(int j=0; j<n; ++j){
                    result.at(j).col(i) = temp.row(j);
                }
            }
            return result;
        }

        std::vector<Eigen::MatrixXd> fux(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
        {
            std::vector<Eigen::MatrixXd> result(n, Eigen::MatrixXd::Zero(m, n));
            Eigen::VectorXd state_0(n), state_1(n);

            // central numerical differentiation
            for(int i=0; i<n; ++i){
                state_0 = state;
                state_1 = state;
                state_0(i) -= eps;
                state_1(i) += eps;
                
                // dfu/dxi = [dfu1/dxi, dfu2/dxi, ... , dfum/dxi]
                auto temp = (fu(state_1, control, time) - fu(state_0, control, time)) / 2 / eps;
                
                // dfu/dxi -> d2fi/dudx
                for(int j=0; j<n; ++j){
                    result.at(j).col(i) = temp.row(j);
                }
            }
            return result;
        }

        std::vector<Eigen::MatrixXd> fuu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
        {
            std::vector<Eigen::MatrixXd> result(n, Eigen::MatrixXd::Zero(m, m));
            Eigen::VectorXd control_0(m), control_1(m);

            // central numerical differentiation
            for(int i=0; i<m; ++i){
                control_0 = control;
                control_1 = control;
                control_0(i) -= eps;
                control_1(i) += eps;

                // dfu/dui = [dfu1/dui, dfu2/dui, ... , dfum/dui]
                auto temp = (fu(state, control_1, time) - fu(state, control_0, time)) / 2 / eps;
                
                // dfu/dxi -> d2fi/du2
                for(int j=0; j<n; ++j){
                    result.at(j).col(i) = temp.row(j);
                }
            }
            return result;
        }
    };

} // namespace DynamicModel

#endif