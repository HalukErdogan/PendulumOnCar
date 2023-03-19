#ifndef COST_FUNCTION_BASE_H
#define COST_FUNCTION_BASE_H

#include <vector>

#include <Eigen/Dense>

namespace CostFunction
{
    class CostFunctionBase
    {
    public:
        virtual ~CostFunctionBase() {}
        virtual double cost(const Eigen::VectorXd &target_state, const std::vector<Eigen::VectorXd> states, const std::vector<Eigen::VectorXd> controls, const std::vector<double> times) const;
        virtual double l(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const = 0;
        virtual Eigen::VectorXd lx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const;
        virtual Eigen::VectorXd lu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const;
        virtual Eigen::MatrixXd lxx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const;
        virtual Eigen::MatrixXd lux(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const;
        virtual Eigen::MatrixXd luu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const;
        virtual double lf(const Eigen::VectorXd &target_state, const Eigen::VectorXd &final_state, const double &time) const = 0;
        virtual Eigen::MatrixXd lfx(const Eigen::VectorXd &target_state, const Eigen::VectorXd &final_state, const double &time) const;
        virtual Eigen::MatrixXd lfxx(const Eigen::VectorXd &target_state, const Eigen::VectorXd &final_state, const double &time) const;
        
    private:
        double eps = 1.0e-4;
    };

    double CostFunctionBase::cost(const Eigen::VectorXd &target_state, const std::vector<Eigen::VectorXd> states, const std::vector<Eigen::VectorXd> controls, const std::vector<double> times) const 
    {   
        double cost = 0;
        int n_sample = states.size();
        
        // running cost
        for(int i=0; i<(n_sample-1); ++i)
        {   
            const Eigen::VectorXd &state = states.at(i);
            const Eigen::VectorXd &control = controls.at(i);
            const double &time = times.at(i);
            const double &next_time = times.at(i+1);
            const double dt = next_time - time;
            cost += l(state, control, time)*dt;
        }

        // terminal cost
        const Eigen::VectorXd &final_state = states.at(n_sample-1);
        const double &final_time = times.at(n_sample-1); 
        cost += lf(target_state, final_state, final_time);

        return cost;
    }

    Eigen::VectorXd CostFunctionBase::lx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const
    {
        int n = state.size();
        Eigen::VectorXd result(n);
        Eigen::VectorXd state_0(n), state_1(n);

        // first order central numerical differentiation
        for(int i=0; i<n; ++i){
            state_0 = state;
            state_1 = state;
            state_0(i) -= eps;
            state_1(i) += eps;
            result(i) = (l(state_1, control, time) - l(state_0, control, time)) / 2.0 / eps;
        }
        return result;
    }

    Eigen::VectorXd CostFunctionBase::lu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const
    {
        int m = control.size();
        Eigen::VectorXd result(m);
        Eigen::VectorXd control_0(m), control_1(m);

        // first order central numerical differentiation
        for(int i=0; i<m; ++i){
            control_0 = control;
            control_1 = control;
            control_0(i) -= eps;
            control_1(i) += eps;
            result(i) = (l(state, control_1, time) - l(state, control_0, time)) / 2.0 / eps;
        }
        return result;
    }

    Eigen::MatrixXd CostFunctionBase::lxx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const
    {
        int n = state.size();
        Eigen::MatrixXd result(n, n);
        Eigen::VectorXd state_0(n), state_1(n);

        // first order central numerical differentiation
        for(int i=0; i<n; ++i){
            state_0 = state;
            state_1 = state;
            state_0(i) -= eps;
            state_1(i) += eps;
            result.col(i) = (lx(state_1, control, time) - lx(state_0, control, time)) / 2.0 / eps;
        }
        return result;
    }

    Eigen::MatrixXd CostFunctionBase::lux(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const
    {
        int n = state.size();
        int m = control.size();
        Eigen::MatrixXd result(m, n);
        Eigen::VectorXd state_0(n), state_1(n);

        // first order central numerical differentiation
        for(int i=0; i<n; ++i){
            state_0 = state;
            state_1 = state;
            state_0(i) -= eps;
            state_1(i) += eps;
            result.col(i) = (lu(state_1, control, time) - lu(state_0, control, time)) / 2.0 / eps;
        }
        return result;
    }

    Eigen::MatrixXd CostFunctionBase::luu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const
    {
        int m = control.size();
        Eigen::MatrixXd result(m, m);
        Eigen::VectorXd control_0(m), control_1(m);

        // first order central numerical differentiation
        for(int i=0; i<m; ++i){
            control_0 = control;
            control_1 = control;
            control_0(i) -= eps;
            control_1(i) += eps;
            result.col(i) = (lu(state, control_1, time) - lu(state, control_0, time)) / 2.0 / eps;
        }
        return result;
    }

    Eigen::MatrixXd CostFunctionBase::lfx(const Eigen::VectorXd &target_state, const Eigen::VectorXd &final_state, const double &time) const
    {
        int n = final_state.size();
        Eigen::VectorXd result(n);
        Eigen::VectorXd final_state_0(n), final_state_1(n);

        // first order central numerical differentiation
        for(int i=0; i<n; ++i){
            final_state_0 = final_state;
            final_state_1 = final_state;
            final_state_0(i) -= eps;
            final_state_1(i) += eps;
            result(i) = (lf(target_state, final_state_1, time) - lf(target_state, final_state_0, time)) / 2.0 / eps;
        }
        return result;
    }

    Eigen::MatrixXd CostFunctionBase::lfxx(const Eigen::VectorXd &target_state, const Eigen::VectorXd &final_state, const double &time) const
    {
        int n = final_state.size();
        Eigen::MatrixXd result(n, n);
        Eigen::VectorXd final_state_0(n), final_state_1(n);

        // first order central numerical differentiation
        for(int i=0; i<n; ++i){
            final_state_0 = final_state;
            final_state_1 = final_state;
            final_state_0(i) -= eps;
            final_state_1(i) += eps;
            result.col(i) = (lfx(target_state, final_state_1, time) - lfx(target_state, final_state_0, time)) / 2.0 / eps;
        }
        return result;
    }

} // namespace ConstFunction

#endif