#include <gtest/gtest.h>
#include <DynamicModel/DynamicModelBase.h>

// Define a concrete class that derives from DynamicModelBase
class TestDynamicModel : public DynamicModel::DynamicModelBase
{
public:
    // Implement the pure virtual function f()
    Eigen::VectorXd f(const Eigen::VectorXd &state, const Eigen::VectorXd &control, const double &time) const override
    {
        // Return a constant vector for testing purposes
        auto const &x0 = state(0);
        auto const &x1 = state(1);
        auto const &u0 = control(0);
        auto const &u1 = control(1);

        Eigen::VectorXd result(2);
        result << x0*x0 + u0*u0*u0, x0*x1 + u0*u1;
        return result;
    }
};

TEST(DynamicModelBaseTest, fxTest)
{
    TestDynamicModel test_model;
    Eigen::VectorXd state(2);
    state << 1.0, 2.0;
    Eigen::VectorXd control(2);
    control << 3.0, 4.0;
    double time = 0.0;
    Eigen::MatrixXd fx_expected(state.size(), state.size());
    fx_expected << 
        2*state(0),     0.0, 
        state(1),       state(0);

    Eigen::MatrixXd fx_result = test_model.fx(state, control, time);

    ASSERT_TRUE(fx_result.isApprox(fx_expected, 1e-5));
}

TEST(DynamicModelBaseTest, fuTest)
{
    TestDynamicModel test_model;
    Eigen::VectorXd state(2);
    state << 1.0, 2.0;
    Eigen::VectorXd control(2);
    control << 3.0, 4.0;
    double time = 0.0;
    Eigen::MatrixXd fu_expected(state.size(), control.size());
    fu_expected << 
        3*control(0)*control(0),    0.0, 
        control(1),                 control(0);

    Eigen::MatrixXd fu_result = test_model.fu(state, control, time);

    ASSERT_TRUE(fu_result.isApprox(fu_expected, 1e-5));
}

TEST(DynamicModelBaseTest, fxxTest)
{
    TestDynamicModel test_model;
    Eigen::VectorXd state(2);
    state << 1.0, 2.0;
    Eigen::VectorXd control(2);
    control << 3.0, 4.0;
    double time = 0.0;
    std::vector<Eigen::MatrixXd> fxx_expected(state.size(), Eigen::MatrixXd(state.size(), state.size()));
    fxx_expected.at(0) << 
        2.0,    0.0,
        0.0,    0.0;

    fxx_expected.at(1) << 
        0.0,    1.0,
        1.0,    0.0;

    std::vector<Eigen::MatrixXd> fxx_result = test_model.fxx(state, control, time);
    for(size_t i = 1; i<fxx_result.size(); ++i){
        ASSERT_TRUE(fxx_result.at(i).isApprox(fxx_expected.at(i), 1e-5));
    }
}

TEST(DynamicModelBaseTest, fuxTest)
{
    TestDynamicModel test_model;
    Eigen::VectorXd state(2);
    state << 1.0, 2.0;
    Eigen::VectorXd control(2);
    control << 3.0, 4.0;
    double time = 0.0;
    std::vector<Eigen::MatrixXd> fux_expected(state.size(), Eigen::MatrixXd(control.size(), state.size()));
    fux_expected.at(0) << 
        0.0,    0.0,
        0.0,    0.0;

    fux_expected.at(1) << 
        0.0,    0.0,
        0.0,    0.0;

    std::vector<Eigen::MatrixXd> fux_result = test_model.fux(state, control, time);
    for(size_t i = 1; i<fux_result.size(); ++i){
        ASSERT_TRUE(fux_result.at(i).isApprox(fux_expected.at(i), 1e-5));
    }
}

TEST(DynamicModelBaseTest, fuuTest)
{
    TestDynamicModel test_model;
    Eigen::VectorXd state(2);
    state << 1.0, 2.0;
    Eigen::VectorXd control(2);
    control << 3.0, 4.0;
    double time = 0.0;
    std::vector<Eigen::MatrixXd> fuu_expected(state.size(), Eigen::MatrixXd(control.size(), control.size()));
    fuu_expected.at(0) << 
        6*control(0),   0.0,
        0.0,            0.0;

    fuu_expected.at(1) << 
        0.0,    1.0,
        1.0,    0.0;

    std::vector<Eigen::MatrixXd> fuu_result = test_model.fuu(state, control, time);
    for(size_t i = 1; i<fuu_result.size(); ++i){
        ASSERT_TRUE(fuu_result.at(i).isApprox(fuu_expected.at(i), 1e-5));
    }
}

// Add more test cases for the other virtual methods of DynamicModelBase

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
