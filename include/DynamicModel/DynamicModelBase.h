#ifndef DYNAMIC_MODEL_BASE_H
#define DYNAMIC_MODEL_BASE_H

#include <vector>

#include <Eigen/Dense>

namespace DynamicModel
{
    class DynamicModelBase
    {
    public:
        virtual ~DynamicModelBase() {}
        virtual Eigen::VectorXd f(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const = 0;
        virtual Eigen::MatrixXd fx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const;
        virtual Eigen::MatrixXd fu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const;
        virtual std::vector<Eigen::MatrixXd> fxx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const;
        virtual std::vector<Eigen::MatrixXd> fux(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const;
        virtual std::vector<Eigen::MatrixXd> fuu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const;

    private:
        double eps = 1.0e-4;
    };

    Eigen::MatrixXd DynamicModelBase::fx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const
    {   
        int n = state.size();
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

    Eigen::MatrixXd DynamicModelBase::fu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const
    {
        int n = state.size();
        int m = control.size();
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

    std::vector<Eigen::MatrixXd> DynamicModelBase::fxx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const
    {
        int n = state.size();
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

    std::vector<Eigen::MatrixXd> DynamicModelBase::fux(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const
    {
        int n = state.size();
        int m = control.size();
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

    std::vector<Eigen::MatrixXd> DynamicModelBase::fuu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const
    {
        int n = state.size();
        int m = control.size();
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

} // namespace DynamicModel

#endif