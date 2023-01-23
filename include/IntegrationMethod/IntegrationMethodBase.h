#ifndef INTEGRATION_METHOD_BASE_H
#define INTEGRATION_METHOD_BASE_H

#include "DynamicModel/DynamicModelBase.h"

namespace IntegrationMethod
{
    using DynamicModel::DynamicModelBase;

    class IntegrationMethodBase
    {
    public:
        virtual ~IntegrationMethodBase() {}

        virtual Eigen::VectorXd integrate(const DynamicModelBase &model, const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time, const double &dt) const = 0;
    };
} // namespace IntegrationMethod

#endif