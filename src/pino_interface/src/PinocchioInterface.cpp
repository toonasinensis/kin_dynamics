
#include "pino_interface/PinocchioInterface.h"
#include "pino_interface/ik_fast_comm.h"
#include <Eigen/Dense>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/spatial/fwd.hpp>

#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Core/util/Constants.h"
#include "pinocchio/multibody/fwd.hpp"

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/

PinocchioInterface::PinocchioInterface(std::string urdf_filename,
                                       bool verbose) {
  // robotModelPtr_ = std::make_shared<pinocchio::Model>(urdf_filename);
  pinocchio::urdf::buildModel(urdf_filename, robot_model);
  robot_data = pinocchio::Data(robot_model);
  // joint_num = robot_model.njoints;

  q_state.resize(robot_model.nq);
  dq_state.resize(robot_model.nv);
  ddq_state.resize(robot_model.nv);
  tau.resize(robot_model.nv);

  if (verbose) {
    // std::cout << "robot_model.supports[1] : " << (int)robot_model.supports[1]
    // << std::endl;

    std::cout << "model name: " << robot_model.name << std::endl;
    std::cout << "nq: " << robot_model.nq << std::endl;
    std::cout << "nv: " << robot_model.nv << std::endl;

    std::cout << "number of joints = " << robot_model.njoints;
    std::cout << "number of frames = " << robot_model.nframes;
    std::cout << "number of bodies = " << robot_model.nbodies;
    std::cout << "lower joint configuration limit = "
              << robot_model.lowerPositionLimit;
    std::cout << "upper joint configuration limit = "
              << robot_model.upperPositionLimit;
    std::cout << "maximal joint velocities = " << robot_model.velocityLimit;
    std::cout << "maximal joint torques = " << robot_model.effortLimit;

    std::cout << "joint name & joint position:" << std::endl;

    for (size_t i = 0; i < robot_model.njoints; i++) {
      std::cout << robot_model.names.at(i) << std::endl;
      std::cout << robot_data.oMi.at(i).translation().transpose() << std::endl;
    }
    std::cout << ("frame name & frame position:") << std::endl;

    // pinocchio::FrameVector frames1 = robot_model.frames;
    std::vector<pinocchio::Frame, Eigen::aligned_allocator<pinocchio::Frame>>
        frames2 = robot_model.frames;
    for (size_t i = 0; i < frames2.size(); ++i) {
      std::cout << frames2.at(i).name << std::endl;
      std::cout << robot_data.oMf.at(i).translation().transpose() << std::endl;
    }
  }
}

void PinocchioInterface::ReadRobotEndEffectName(std::vector<std::string> eename,
                                                bool verbose) {
  endEffectName = eename;
  EndEffect_num = endEffectName.size();

  allj2fbase_jcb.resize(6, robot_model.nv);
  allj2fbase_jcb_local.resize(6, robot_model.nv);

  allj2EE_J.resize(6 * EndEffect_num, robot_model.nv);
  allj2EE_J_local.resize(6 * EndEffect_num, robot_model.nv);

  allj2EE_jcb.resize(EndEffect_num * 3, robot_model.nv);
  allj2EE_jcb_local.resize(EndEffect_num * 3, robot_model.nv);

  x_state.resize(EndEffect_num * 3, 1);
  dx_state.resize(EndEffect_num * 3, 1);
  ddx_state.resize(EndEffect_num * 3, 1);

  hip_offset_local.resize(EndEffect_num * 3, 1);

  x_state_local.resize(EndEffect_num * 3, 1);
  dx_state_local.resize(EndEffect_num * 3, 1);
  ddx_state_local.resize(EndEffect_num * 3, 1);

  if (verbose) {
    for (size_t i = 0; i < EndEffect_num; i++) {
      std::cout << "load endeffect name ,NO." << i << "   "
                << endEffectName.at(i) << std::endl;
    }
    std::cout << "size of allj2fbase_jcb" << allj2fbase_jcb.size() << std::endl;
    std::cout << "size of allj2ee_J" << allj2EE_J.size() << std::endl;
    std::cout << "size of allj2EE_jcb" << allj2EE_jcb.size() << std::endl;
    std::cout << "size of x_state" << x_state.size() << std::endl;
    std::cout << "size of dx_state" << dx_state.size() << std::endl;
    std::cout << "size of ddx_state" << ddx_state.size() << std::endl;
  }
}

void PinocchioInterface::ReadRobotTrunkName(std::string bodyname,
                                            bool verbose) {
  TrunkName = bodyname;
  if (verbose) {
    std::cout << "load body name: " << bodyname << std::endl;
  }
}

void PinocchioInterface::ReadOffsetName(std::vector<std::string> offsetname,
                                        bool verbose) {
  OffsetName = offsetname;
  if (OffsetName.size() != endEffectName.size()) {
    std::cout
        << "The size of offsetname is not equal to the size of endEffectName"
        << std::endl;
    exit(0);
  }
  q_state.setZero();
  q_state.block(3, 0, 4, 1) = Eigen::Vector4d(0, 0, 0, 1);
  dq_state.setZero();
  ddq_state.setZero();

  pinocchio::ReferenceFrame rf = pinocchio::LOCAL_WORLD_ALIGNED;
  pinocchio::forwardKinematics(robot_model, robot_data, q_state, dq_state,
                               ddq_state);
  pinocchio::updateFramePlacements(robot_model, robot_data);

  for (size_t i = 0; i < EndEffect_num; i++) {
    // get ee pos in inertia frame
    auto Frame_id = robot_model.getFrameId(OffsetName.at(i));
    // std::cout << robot_data.oMf[Frame_id].translation() << std::endl;
    hip_offset_local.block(3 * i, 0, 3, 1) =
        robot_data.oMf[Frame_id].translation();
  }
  for (size_t i = 0; i < EndEffect_num; i++) {
    std::cout << "load endeffect name ,NO." << i << "   " << OffsetName.at(i)
              << std::endl;
    std::cout << "hip_offset_local "
              << hip_offset_local.block(3 * i, 0, 3, 1).transpose()
              << std::endl;
  }
}

void PinocchioInterface::forward_kin(const Eigen::VectorXd &q_,
                                     const Eigen::VectorXd &dq_,
                                     const Eigen::VectorXd &ddq_,
                                     bool verbose) {
  if (q_.size() != robot_model.nq && dq_.size() != robot_model.nv &&
      ddq_.size() != robot_model.nv) {
    std::cout << "q_.size(): " << q_.size()
              << "robot_model.nq: " << robot_model.nq << "dq_.size()"
              << dq_.size() << " robot_model.nv" << robot_model.nv << std::endl;
    throw std::runtime_error(" the input is not equal the robotconfig!!!!!!!");
    std::cerr << " the input is not equal the robotconfig!!!!!!!" << std::endl;
    exit(0);
  }

  q_state = q_;
  dq_state = dq_;
  ddq_state = ddq_;
  pinocchio::ReferenceFrame rf = pinocchio::LOCAL_WORLD_ALIGNED;

  pinocchio::forwardKinematics(robot_model, robot_data, q_state, dq_state,
                               ddq_state);
  pinocchio::computeJointJacobians(robot_model, robot_data, q_state);
  pinocchio::computeJointJacobiansTimeVariation(robot_model, robot_data,
                                                q_state, dq_state);
  pinocchio::updateFramePlacements(robot_model, robot_data);
  // pinocchio::framesForwardKinematics(model, data, q);

  pinocchio::crba(robot_model, robot_data, q_state);
  pinocchio::nonLinearEffects(robot_model, robot_data, q_state, dq_state);

  for (size_t i = 0; i < EndEffect_num; i++) {
    // get ee pos in inertia frame
    auto Frame_id = robot_model.getFrameId(endEffectName.at(i));
    // std::cout << robot_data.oMf[Frame_id].translation() << std::endl;
    x_state.block(3 * i, 0, 3, 1) = robot_data.oMf[Frame_id].translation();

    // get ee vel in inertia frame
    dx_state.block(3 * i, 0, 3, 1) =
        pinocchio::getFrameVelocity(robot_model, robot_data, Frame_id, rf)
            .linear();

    // get ee acc in inertia frame
    ddx_state.block(3 * i, 0, 3, 1) =
        pinocchio::getFrameAcceleration(robot_model, robot_data, Frame_id, rf)
            .linear();

    Eigen::Matrix<double, 6, Eigen::Dynamic> J_base;
    J_base.resize(6, robot_model.nv);
    J_base.setZero();

    pinocchio::computeFrameJacobian(robot_model, robot_data, q_state, Frame_id,
                                    rf, J_base);
    // std::cout << J_base << "J_base" << std::endl;

    allj2EE_J.block(6 * i, 0, 6, robot_model.nv) = J_base;
    allj2EE_jcb.block(3 * i, 0, 3, robot_model.nv) =
        J_base.block(0, 0, 3, robot_model.nv);

    if (verbose) {
      std::cout << "the end effector's state :" << endEffectName.at(i) << "\n"
                << "x_state"
                << "\n"
                << x_state.block(3 * i, 0, 3, 1) << "\n"
                << "dx_state"
                << "\n"
                << dx_state.block(3 * i, 0, 3, 1) << "\n"
                << "ddx_state"
                << "\n"
                << ddx_state.block(3 * i, 0, 3, 1) << std::endl;
    }
  }
  if (verbose) {
    for (int i = 0; i < EndEffect_num; i++) {
      if (dx_state.block(3 * i, 0, 3, 1)
              .isApprox(
                  (allj2EE_jcb.block(3 * i, 0, 3, robot_model.nv) * dq_))) {
        std::cout << "dx = J * dq" << std::endl;

        std::cout << endEffectName.at(i) << " dx == j*dq" << std::endl;
      } else {

        std::cout << "dx_state.block(3 * i, 0, 3, 1)"
                  << dx_state.block(3 * i, 0, 3, 1) << std::endl;
        std::cout << "allj2EE_jcb.block(3 * i, 0, 3, robot_model.nv) "
                  << allj2EE_jcb.block(3 * i, 0, 3, robot_model.nv)
                  << std::endl;
        std::cout << "dq_" << dq_ << std::endl;
        std::cout << "j*dq "
                  << allj2EE_jcb.block(3 * i, 0, 3, robot_model.nv) * dq_
                  << std::endl;
        std::cerr << "dx !=j*dq calucae error !" << std::endl;
        exit(0);
      }
    }
  }

  auto body_id = robot_model.getFrameId(TrunkName);
  allj2fbase_jcb.setZero();
  if (verbose) {
    std::cout << TrunkName << std::endl;
    std::cout << "q_state" << q_state << std::endl;
    std::cout << "allj2fbase_jcb" << allj2fbase_jcb << std::endl;
  }
  pinocchio::computeFrameJacobian(robot_model, robot_data, q_state, body_id, rf,
                                  allj2fbase_jcb);
  // std::cout << "allj2fbase_jcb  \n" << allj2fbase_jcb << std::endl;
  assert(allj2fbase_jcb.size() == 6 * robot_model.nv);

  // test the body!!!!
  if (verbose) {
    // std::cout << "allj2fbase_jcb*dq" << allj2fbase_jcb * dq_state <<
    // std::endl;
    if (dq_.block(0, 0, 6, 1).isApprox(allj2fbase_jcb * dq_state)) {
      std::cout << "dx_body == j * dq " << std::endl;
    } else {

      std::cout << "body" << dq_state.block(6, 0, 6, 1) << std::endl;
      std::cout << "j" << allj2fbase_jcb << std::endl;
      std::cout << "dq" << dq_state << std::endl;

      std::cout << "j*dq" << allj2fbase_jcb * dq_state << std::endl;
      std::cout << "dx_body != j * dq !!!!calucae error !" << std::endl;
      exit(0);
    }
  }
  //   std::cout << "allj2EE_jcb\n" << allj2EE_jcb << std::endl;
  Eigen::Vector4d quad = q_state.block(3, 0, 4, 1);
  Eigen::Vector3d pos = q_state.block(0, 0, 3, 1);
  T_W2B = pinocchio::SE3(Eigen::Quaterniond(quad), pos);

  if (verbose) {
    std::cout << "the body's state : quad"
              << "\n"
              << quad << "\n"
              << "pos" << pos << std::endl;
  }
}

void PinocchioInterface::forward_kin_local(const Eigen::VectorXd &q_,
                                           const Eigen::VectorXd &dq_,
                                           const Eigen::VectorXd &ddq_,
                                           bool verbose) {
  if (q_.size() != robot_model.nq && dq_.size() != robot_model.nv &&
      ddq_.size() != robot_model.nv) {
    std::cout << "q_.size(): " << q_.size()
              << "robot_model.nq: " << robot_model.nq << "dq_.size()"
              << dq_.size() << " robot_model.nv" << robot_model.nv << std::endl;
    throw std::runtime_error(" the input is not equal the robotconfig!!!!!!!");
    std::cerr << " the input is not equal the robotconfig!!!!!!!" << std::endl;
    exit(0);
  }

  q_state = q_;
  q_state.block(3, 0, 3, 1) = Eigen::Vector3d(0, 0, 0);
  q_state.block(4, 0, 4, 1) = Eigen::Vector4d(0, 0, 0, 1);
  dq_state = dq_;
  ddq_state = ddq_;
  pinocchio::ReferenceFrame rf = pinocchio::LOCAL_WORLD_ALIGNED;

  pinocchio::forwardKinematics(robot_model, robot_data, q_state, dq_state,
                               ddq_state);
  pinocchio::computeJointJacobians(robot_model, robot_data, q_state);
  pinocchio::computeJointJacobiansTimeVariation(robot_model, robot_data,
                                                q_state, dq_state);
  pinocchio::updateFramePlacements(robot_model, robot_data);
  // pinocchio::framesForwardKinematics(model, data, q);

  pinocchio::crba(robot_model, robot_data, q_state);
  pinocchio::nonLinearEffects(robot_model, robot_data, q_state, dq_state);

  for (size_t i = 0; i < EndEffect_num; i++) {
    // get ee pos in inertia frame
    auto Frame_id = robot_model.getFrameId(endEffectName.at(i));
    // std::cout << robot_data.oMf[Frame_id].translation() << std::endl;
    x_state_local.block(3 * i, 0, 3, 1) =
        robot_data.oMf[Frame_id].translation();

    // get ee vel in inertia frame
    dx_state_local.block(3 * i, 0, 3, 1) =
        pinocchio::getFrameVelocity(robot_model, robot_data, Frame_id, rf)
            .linear();

    // get ee acc in inertia frame
    ddx_state_local.block(3 * i, 0, 3, 1) =
        pinocchio::getFrameAcceleration(robot_model, robot_data, Frame_id, rf)
            .linear();

    Eigen::Matrix<double, 6, Eigen::Dynamic> J_base_local;
    J_base_local.resize(6, robot_model.nv);
    J_base_local.setZero();

    pinocchio::computeFrameJacobian(robot_model, robot_data, q_state, Frame_id,
                                    rf, J_base_local);
    // std::cout << J_base << "J_base" << std::endl;

    allj2EE_J_local.block(6 * i, 0, 6, robot_model.nv) = J_base_local;
    allj2EE_jcb_local.block(3 * i, 0, 3, robot_model.nv) =
        J_base_local.block(0, 0, 3, robot_model.nv);

    if (verbose) {
      std::cout << "the end effector's state :" << endEffectName.at(i) << "\n"
                << "x_state_local"
                << "\n"
                << x_state_local.block(3 * i, 0, 3, 1) << "\n"
                << "dx_state_local"
                << "\n"
                << dx_state_local.block(3 * i, 0, 3, 1) << "\n"
                << "ddx_state_local"
                << "\n"
                << ddx_state_local.block(3 * i, 0, 3, 1) << std::endl;
    }
  }
  if (verbose) {
    for (int i = 0; i < EndEffect_num; i++) {
      if (dx_state_local.block(3 * i, 0, 3, 1)
              .isApprox((allj2EE_jcb_local.block(3 * i, 0, 3, robot_model.nv) *
                         dq_))) {
        std::cout << "dx_local = J_local * dq" << std::endl;

        std::cout << endEffectName.at(i) << " dx_local == j_local * dq"
                  << std::endl;
      } else {
        std::cout << "dx_state.block(3 * i, 0, 3, 1)"
                  << dx_state_local.block(3 * i, 0, 3, 1) << std::endl;
        std::cout << "allj2EE_jcb.block(3 * i, 0, 3, robot_model.nv) "
                  << allj2EE_jcb_local.block(3 * i, 0, 3, robot_model.nv)
                  << std::endl;
        std::cout << "dq_" << dq_ << std::endl;
        std::cout << "j*dq "
                  << allj2EE_jcb_local.block(3 * i, 0, 3, robot_model.nv) * dq_
                  << std::endl;
        std::cerr << "dx !=j*dq calucae error !" << std::endl;
        exit(0);
      }
    }
  }

  auto body_id = robot_model.getFrameId(TrunkName);
  allj2fbase_jcb_local.setZero();
  std::cout << TrunkName << std::endl;
  std::cout << "q_state" << q_state << std::endl;
  std::cout << "allj2fbase_jcb" << allj2fbase_jcb_local << std::endl;
  pinocchio::computeFrameJacobian(robot_model, robot_data, q_state, body_id, rf,
                                  allj2fbase_jcb_local);
  // std::cout << "allj2fbase_jcb  \n" << allj2fbase_jcb << std::endl;
  assert(allj2fbase_jcb_local.size() == 6 * robot_model.nv);
}

void PinocchioInterface::inv_dyn(const Eigen::VectorXd &q_,
                                 const Eigen::VectorXd &dq_,
                                 const Eigen::VectorXd &ddq_,
                                 const Eigen::VectorXd &F_ext, bool verbose) {
  if (F_ext.size() != EndEffect_num * 3) {
    throw std::runtime_error(
        " the ext force  is not equal the number of  end effector!!!!!!!");
    std::cerr
        << " the ext force  is not equal the number of  end effector!!!!!!!"
        << std::endl;
    exit(0);
  }
  forward_kin(q_, dq_, ddq_);
  pinocchio::rnea(robot_model, robot_data, q_, dq_, ddq_);
  tau = robot_data.tau;

  tau += allj2EE_jcb.transpose() * F_ext;
  if (verbose) {
    std::cout << "tau" << tau << std::endl;
  }
}

// void inv_kin(const Eigen::VectorXd &p_leg, Eigen::Vector3d q_fdb, int
// leg_num);

Eigen::Vector3d
PinocchioInterface::inv_kin_local(const Eigen::Vector3d &p_foot_local,
                                  Eigen::Vector3d q_fdb, int leg_num,
                                  bool verbose) {
  Eigen::Vector3d q, foot_local_hip;
  IkSolutionList<IkReal> solutions;
  Eigen::Vector3d hipoffset = hip_offset_local.block(3 * leg_num, 0, 3, 1);
  double trans[3], rot[9];
  // std::vector<IkReal> vfree(GetNumFreeParameters());

  foot_local_hip = p_foot_local - hipoffset;
  trans[0] = foot_local_hip[0];
  trans[1] = foot_local_hip[1];
  trans[2] = foot_local_hip[2];

  std::cout << "leg_num" << leg_num << std::endl;
  if (leg_num == 0) {
    auto bSuccess = el_mini_ComputeIk_lb(trans, rot, NULL, solutions);
  } else if (leg_num == 1 || leg_num == 2) {
    auto bSuccess = el_mini_ComputeIK_left(trans, rot, NULL, solutions);
  } else if (leg_num == 3) {
    auto bSuccess = el_mini_ComputeIK_rb(trans, rot, NULL, solutions);

  } else {
    auto bSuccess = el_mini_ComputeIK_right(trans, rot, NULL, solutions);
  }
  if (verbose) {
    // std::cout << "solutions.nr" << solutions.nrSolutions() << std::endl;
    std::cout << "eeName:" << leg_num << endEffectName.at(leg_num) << std::endl;
    std::cout << "p_foot_local" << p_foot_local.transpose() << std::endl;
    std::cout << "hipoffset" << hipoffset << std::endl;
    std::cout << "foot_local_hip " << trans[0] << " " << trans[1] << " "
              << trans[2] << std::endl;
    printf("Found %d ik solutions:\n", (int)solutions.GetNumSolutions());
  }

  double min = 1000;
  int solution_num = 0;
  for (std::size_t i = 0; i < solutions.GetNumSolutions(); ++i) {
    std::vector<IkReal> solvalues(GetNumJoints());
    const IkSolutionBase<IkReal> &sol = solutions.GetSolution(i);

    std::vector<IkReal> vsolfree(sol.GetFree().size());
    sol.GetSolution(&solvalues[0], NULL);
    for (std::size_t j = 0; j < solvalues.size(); ++j) {
      if (verbose) {
        std::cout << solvalues[j] << " ";
      }
    }
    std::cout << std::endl;
    Eigen::Vector3d temp_solution(solvalues[0], solvalues[1], solvalues[2]);
    if (min > (q_fdb - temp_solution).norm()) {
      min = (q_fdb - temp_solution).norm();
      solution_num = i;
    }
  }
  if (solutions.GetNumSolutions() != 0) {
    const IkSolutionBase<IkReal> &sol = solutions.GetSolution(solution_num);
    std::vector<IkReal> real_solution(GetNumJoints());
    std::vector<IkReal> vsolfree(sol.GetFree().size());
    sol.GetSolution(&real_solution[0],
                    vsolfree.size() > 0 ? &vsolfree[0] : NULL);
    q[0] = real_solution[0];
    q[1] = real_solution[1];
    q[2] = real_solution[2];
    return q;
  } else {
    return q_fdb;
  }
}
