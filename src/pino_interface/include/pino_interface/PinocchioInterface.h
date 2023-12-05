
#pragma once

// #include <Eigen/Eigen>
#include <Eigen/Dense>
#include <iosfwd>
#include <memory>
#include <pinocchio/fwd.hpp> // always include it before any other header
#include <string>
// #include <type_traits>
// #include <vector>

// #include "common/models/EndEffectorKinematics.h"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/multibody/fwd.hpp"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/spatial/explog.hpp"

class PinocchioInterface
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  public:
    PinocchioInterface() = default;
    explicit PinocchioInterface(std::string urdf_filename, bool verbose = false);

    ~PinocchioInterface() = default;

    void ReadRobotEndEffectName(std::vector<std::string> eename, bool verbose = false);
    void ReadRobotTrunkName(std::string bodyname, bool verbose = false);
    void ReadOffsetName(std::vector<std::string> offsetname, bool verbose = false);

    void forward_kin(const Eigen::VectorXd &q_, const Eigen::VectorXd &dq_,
                     const Eigen::VectorXd &ddq_, bool verbose = false);
    void forward_kin_local(const Eigen::VectorXd &q_, const Eigen::VectorXd &dq_,
                           const Eigen::VectorXd &ddq_, bool verbose = false);
    void inv_dyn(const Eigen::VectorXd &q_, const Eigen::VectorXd &dq_, const Eigen::VectorXd &ddq_,
                 const Eigen::VectorXd &F_ext, bool verbose = false);

    Eigen::Vector3d inv_kin_local(const Eigen::Vector3d &p_foot_local, Eigen::Vector3d q_fdb,
                                  int leg_num, bool verbose = false);
    pinocchio::Model robot_model;
    pinocchio::Data robot_data;

    // q
    Eigen::VectorXd q_state, dq_state, ddq_state;
    // tau
    Eigen::VectorXd tau;
    // x
    Eigen::VectorXd x_state, dx_state, ddx_state;
    Eigen::VectorXd x_state_local, dx_state_local, ddx_state_local;
    Eigen::VectorXd hip_offset_local;

    // jocobian
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
        allj2EE_jcb; // alljoint to endeffect jocobian
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> allj2EE_jcb_local;

    Eigen::Matrix<double, 6, Eigen::Dynamic> allj2fbase_jcb;       // alljoint to floatbase jocobian
    Eigen::Matrix<double, 6, Eigen::Dynamic> allj2fbase_jcb_local; // alljoint to floatbase jocobian

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
        allj2EE_J; // 6-D vector ,but we only need the linear sector
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
        allj2EE_J_local;  // 6-D vector ,but we only need the linear sector
    pinocchio::SE3 T_W2B; // the T in inertia frame

  private:
    int EndEffect_num;
    std::vector<std::string> endEffectName;
    std::vector<std::string> OffsetName;

    std::string TrunkName;
};
