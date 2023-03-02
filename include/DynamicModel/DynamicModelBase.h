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
        virtual Eigen::MatrixXd fx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const = 0;
        virtual Eigen::MatrixXd fu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const = 0;
        virtual std::vector<Eigen::MatrixXd> fxx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const = 0;
        virtual std::vector<Eigen::MatrixXd> fux(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const = 0;
        virtual std::vector<Eigen::MatrixXd> fuu(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const = 0;
    };
} // namespace DynamicModel

#endif