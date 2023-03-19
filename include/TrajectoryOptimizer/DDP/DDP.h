#ifndef DDP_H
#define DDP_H

#include "TrajectoryOptimizer/TrajectoryOptimizerBase.h"

#include <iostream>

namespace TrajectoryOptimizer
{
    class DDP : public TrajectoryOptimizerBase
    {
    public:
        Solution solve(const Problem &problem, const DynamicModelBase &model, const IntegrationMethodBase &integrator,  const ODESolverBase &ode) override
        {
            // config params
            int n_iter = 200;
            double alpha = 0.01;
            
            // number of node
            int n_sample = problem.T.size();

            // number of states
            int n_state = problem.X0.at(0).rows();
            int n_control = problem.U0.at(0).rows();

            // set up the data
            std::vector<Node> data(n_sample);

            // fill time, state, and control inputs
            for(int i=0; i<n_sample; ++i){
                data.at(i).t = problem.T.at(i);
                data.at(i).x = problem.X0.at(i);
                data.at(i).u = problem.U0.at(i);
            }

            // set initial state
            data.at(0).x = problem.x0;

            // cost of the trajectory
            double prev_cost = __DBL_MAX__;

            // iterate until converge
            for(int iter=0; iter<n_iter; ++iter){
                // simulate dynamics and calculate cost functions
                for(int i=0; i<(n_sample-1); ++i){
                    auto &node = data.at(i);
                    auto &next_node = data.at(i+1);
                    const auto &state = node.x;
                    const auto &control = node.u;
                    const auto &time = node.t;
                    const auto &next_time = next_node.t;
                    const auto &dt = next_time - time;

                    // next state
                    next_node.x = integrator.integrate(model, state, control, time, dt);

                    // continuous dynamics to discete dynamics
                    node.fx = Eigen::MatrixXd::Identity(n_state, n_state) + model.fx(state, control, time) * dt;
                    node.fu = model.fu(state, control, time) * dt ;
                    node.fxx = model.fxx(state, control, time);
                    std::transform(node.fxx.begin(), node.fxx.end(), node.fxx.begin(), [&dt](const auto &fixx){return fixx*dt;});
                    node.fux = model.fux(state, control, time);
                    std::transform(node.fux.begin(), node.fux.end(), node.fux.begin(), [&dt](const auto &fiux){return fiux*dt;});
                    node.fuu = model.fuu(state, control, time);
                    std::transform(node.fuu.begin(), node.fuu.end(), node.fuu.begin(), [&dt](const auto &fiuu){return fiuu*dt;});

                    // cost function
                    node.lx = problem.Q*state*dt;
                    node.lu = problem.R*control*dt;
                    node.lxx = problem.Q*dt;
                    node.lux = Eigen::MatrixXd::Zero(n_control, n_state)*dt;
                    node.luu = problem.R*dt;
                }

                // Set final values of value function
                const Eigen::VectorXd &final_state = data.at(n_sample-1).x;
                const Eigen::VectorXd &target_state = problem.xt;
                data.at(n_sample-1).v = 0.5 * (final_state-target_state).transpose()*problem.Qf*(final_state-target_state);
                data.at(n_sample-1).vx = problem.Qf*(final_state-target_state);
                data.at(n_sample-1).vxx = problem.Qf;

                // backward pass
                for(int i=n_sample-2; i>=0; --i){
                    Eigen::VectorXd qx;
                    Eigen::VectorXd qu;
                    Eigen::MatrixXd qxx;
                    Eigen::MatrixXd qux;
                    Eigen::MatrixXd quu;
                    Eigen::MatrixXd quu_inv;

                    auto &node = data.at(i);
                    auto &next_node = data.at(i+1);
                    
                    // calculate the intermediate variables
                    qx = node.lx + node.fx.transpose()*next_node.vx;
                    qu = node.lu + node.fu.transpose()*next_node.vx;
                    
                    qxx = node.lxx 
                        + node.fx.transpose()*next_node.vxx*node.fx;
                    // for(int j=0; j<node.x.size(); ++j){
                    //     qxx += next_node.vx(j) * node.fxx.at(j); 
                    // }
                    
                    qux = node.lux 
                        + node.fu.transpose()*next_node.vxx*node.fx;
                    // for(int j=0; j<node.x.size(); ++j){
                    //     qux += next_node.vx(j) * node.fux.at(j); 
                    // }
                    
                    quu = node.luu 
                        + node.fu.transpose()*next_node.vxx*node.fu;
                    // for(int j=0; j<node.x.size(); ++j){
                    //     quu += next_node.vx(j) * node.fuu.at(j); 
                    // }
                    
                    // calculate the control gains
                    quu_inv = quu.inverse();
                    node.k = - quu_inv*qu;
                    node.K = - quu_inv*qux;

                    // calculate the value function
                    node.v = next_node.v + node.k.transpose()*qu + 0.5*node.k.transpose()*quu*node.k;
                    node.vx = qx + node.K.transpose()*qu + qux.transpose()*node.k + node.K.transpose()*quu*node.k;
                    node.vxx = qxx + node.K.transpose()*qux + qux.transpose()*node.K + node.K.transpose()*quu*node.K;
                }

                // forward pass
                Eigen::VectorXd x_hat = data.at(0).x;
                Eigen::VectorXd u_hat = Eigen::VectorXd::Zero(n_control);
                for(int i=0; i<n_sample-1; ++i){
                    auto &node = data.at(i);
                    auto &next_node = data.at(i+1);
                    const double &time = node.t;
                    const double &next_time = next_node.t;
                    const double delta_time = next_time - time;

                    // update the delta control
                    node.du = alpha * (node.k + node.K*(x_hat - node.x));

                    // calculate the new control
                    u_hat = node.u + node.du;

                    // calculate the new state
                    x_hat = integrator.integrate(model, x_hat, u_hat, time, delta_time);
                }
                
                // Line search
                for(int i=0; i<n_sample-1; ++i){
                    auto &node = data.at(i);
                    node.u += node.du;
                    // std::cout << node.u << std::endl;
                }

                // Calculate the cost
                double cost = 0.0;
                for(int i=0; i<n_sample-1; ++i){
                    auto &node = data.at(i);
                    auto &next_node = data.at(i+1);
                    double dt = next_node.t - node.t;
                    cost += (node.x.transpose()*problem.Q*node.x + node.u.transpose()*problem.R*node.u).value() * dt;
                }
                cost += ((final_state - target_state).transpose()*problem.Qf*(final_state - target_state)).value();
                std::cout << "Iteration: " << iter << ", " << "Cost: " << cost << std::endl;
                
                if(std::abs(prev_cost - cost) < 1e-5){
                    std::cout << "DDP is converged!" << std::endl;
                    break;
                }

                prev_cost = cost;
            }

            Solution sol;
            for(int i=0; i<n_sample; ++i){
                auto &node = data.at(i);
                sol.X.push_back(node.x);
                sol.U.push_back(node.u);
            }
            return sol;
        }

    private:
        struct Node
        {
            // time
            double t;

            // state
            Eigen::VectorXd x;

            // control
            Eigen::VectorXd u;
            Eigen::VectorXd du;

            // descete dynamics
            Eigen::MatrixXd fx;
            Eigen::MatrixXd fu;
            std::vector<Eigen::MatrixXd> fxx;
            std::vector<Eigen::MatrixXd> fux;
            std::vector<Eigen::MatrixXd> fuu;

            // cost function
            Eigen::VectorXd lx;
            Eigen::VectorXd lu;
            Eigen::MatrixXd lxx;
            Eigen::MatrixXd lux;
            Eigen::MatrixXd luu;

            // value function
            double v;
            Eigen::VectorXd vx;
            Eigen::MatrixXd vxx;

            // intermediate calculations
            // Eigen::VectorXd qx;
            // Eigen::VectorXd qu;
            // Eigen::MatrixXd qxx;
            // Eigen::MatrixXd qux;
            // Eigen::MatrixXd quu;

            // control gains
            Eigen::VectorXd k;  // feedforward gain
            Eigen::MatrixXd K;  // feedback gain
        };
    };

} // namespace TrajectoryOptimizer

#endif