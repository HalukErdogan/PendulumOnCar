#ifndef DYNAMIC_MODEL_BASE_H
#define DYNAMIC_MODEL_BASE_H

#include <eigen3/Eigen/Dense>

namespace DynamicModel
{
    class DynamicModelBase
    {
    public:
        virtual ~DynamicModelBase() {}

        virtual Eigen::VectorXd f(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const = 0;
        virtual Eigen::MatrixXd df_dx(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const = 0;
        virtual Eigen::MatrixXd df_du(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const = 0;
    };
} // namespace DynamicModel

#endif