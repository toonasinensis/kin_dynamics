// #include <Eigen/src/Core/VectorBlock.h>
#include <cmath>
#include <memory>
#include <ostream>
#include <string>

// #include "Eigen/src/Core/util/Constants.h"
#include "iostream"
#include "pino_interface/PinocchioInterface.h"
#include "string"

#include "iostream"
#include <geometry_msgs/TransformStamped.h>
#include <geometry_msgs/Twist.h>
#include <ros/ros.h>
#include <sensor_msgs/JointState.h>
#include <string>
#include <tf2/LinearMath/Quaternion.h>
#include <tf2_ros/transform_broadcaster.h>
#include <turtlesim/Pose.h>

#include <tf2_ros/transform_listener.h>
#include <thread>
#include <vector>
#define BACKWARD_HAS_DW 1
#include "backward.hpp"
#include <unistd.h>

namespace backward {
backward::SignalHandling sh;
}

std::string quadruped_path =
    "/home/hexapod/Software/HexapodSoftware/xcy/a1_pinocchio.urdf";

std::string hex_path =
    "/home/hexapod/xie/src/robot/urdf/el_mini/urdf/el_mini.urdf";

using namespace std;
using namespace pinocchio;
std::vector<std::string> Hex_EEname = {"LB_FOOT", "LF_FOOT", "LM_FOOT",
                                       "RB_FOOT", "RF_FOOT", "RM_FOOT"};

std::vector<std::string> Hex_HIPname = {"LB_HIP", "LF_HIP", "LM_HIP",
                                        "RB_HIP", "RF_HIP", "RM_HIP"};

std::vector<std::string> quadruped_EEname = {"FL_FOOT", "FR_FOOT", "RL_FOOT",
                                             "RR_FOOT"};

std::vector<std::string> Hex_prex = {"LB", "LF", "LM", "RB", "RF", "RM"};
std::vector<std::string> JointName = {"HAA", "HFE", "KFE"};
std::vector<std::string> Hex_JointName;
bool if_frist_in = false;
Eigen::Matrix<double, 3, 6> foot_x;
std::vector<double> q_limit_lower = {
    -9999, -9999, -9999, -1.01, -1.01, -1.01, -1.01, -1.5, 0, 0,    -1.5, 0, 0,
    -1.5,  0,     0,     -1.5,  0,     0,     -1.5,  0,    0, -1.5, 0,    0};
std::vector<double> q_limit_upper = {
    999, 999, 999, 1.01, 1.01, 1.01, 1.01, 1.5, 1.2, 1.2, 1.5, 1.2, 1.2,
    1.5, 1.2, 1.2, 1.5,  1.2,  1.2,  1.5,  1.2, 1.2, 1.5, 1.2, 1.2};
double max_error;
int main(int argc, char **argv) {
  ros::init(argc, argv, "test_pino");
  ros::NodeHandle n;
  ros::Publisher joint_pub =
      n.advertise<sensor_msgs::JointState>("joint_states", 1);
  // tf::TransformBroadcaster broadcaster;
  ros::Rate loop_rate(10);
  tf2_ros::TransformBroadcaster Broadcaster_odom;

  // message declarations
  geometry_msgs::TransformStamped odom_trans;
  sensor_msgs::JointState joint_state;
  odom_trans.header.frame_id = "odom";
  odom_trans.child_frame_id = "BASE";
  tf2_ros::Buffer tfBuffer;
  tf2_ros::TransformListener tfListener(tfBuffer);
  geometry_msgs::TransformStamped transformStamped;

  double yaw = 0, pitch = 0, roll = 0;
  int num = 6;
  for (int i = 0; i < num; i++) {
    for (int j = 0; j < 3; j++) {
      Hex_JointName.push_back(Hex_prex.at(i) + "_" + JointName.at(j));
    }
  }

  std::string Hex_Bodyname = "BASE";
  std::string quadruped_Bodyname = "trunk";
  std::shared_ptr<PinocchioInterface> p_ptr = nullptr;
  static long double time = 0;
  if (num == 6) {

    joint_state.name.resize(18);
    joint_state.position.resize(18);

    p_ptr = std::make_shared<PinocchioInterface>(hex_path, true);
    PinocchioInterface pin(hex_path, true);
    pin.ReadRobotEndEffectName(Hex_EEname, true);
    pin.ReadRobotTrunkName(Hex_Bodyname, true);
    pin.ReadOffsetName(Hex_HIPname, true);
    // Eigen::Vector3d leg_inv_test_right(0.10679, -0.16418, -0.15286);
    // exit(0);
    Eigen::Matrix<double, 25, 1> q;
    Eigen::Matrix<double, 24, 1> dq;
    Eigen::Matrix<double, 24, 1> ddq;
    Eigen::Matrix<double, 18, 1> F;
    double theta = 0;

    while (ros::ok()) {
      q.setZero();
      dq.setZero();
      ddq.setZero();
      F.setZero();
      time += 0.0002;
      if (time > 1) {
        time = 0;
      }
      theta = 1.2 * sin(time * M_PI / 2);
      //   theta = 0.3;
      Eigen::Vector3d ea0(theta, theta, theta);
      Eigen::Matrix3d R;
      R = Eigen::AngleAxisd(ea0[2], ::Eigen::Vector3d::UnitZ()) *
          Eigen::AngleAxisd(ea0[1], ::Eigen::Vector3d::UnitY()) *
          Eigen::AngleAxisd(ea0[0], ::Eigen::Vector3d::UnitX());
      Eigen::Quaterniond quad;
      quad = R;

      Eigen::Vector4d vquad;

      // vquad(0) = quad.x();
      // vquad(1) = quad.y();
      // vquad(2) = quad.z();
      // vquad(3) = quad.w();

      vquad(0) = 0;
      vquad(1) = 0;
      vquad(2) = 0;
      vquad(3) = 1;
      for (int i = 0; i < 18; i++) {

        q(i + 7) = theta;
        if (q(i + 7) > q_limit_upper.at(i + 7)) {
          q(i + 7) = q_limit_upper.at(i + 7);
        }
        if (q(i + 7) < q_limit_lower.at(i + 7)) {
          q(i + 7) = q_limit_lower.at(i + 7);
        }
        dq(i + 6) = 0.0 * i;
      }
      joint_state.header.stamp = ros::Time::now();
      // q(0) = theta;
      // q(1) = theta;
      // q(2) = theta;
      // set visualize struct
      for (int i = 0; i < 18; i++) {
        joint_state.name[i] = Hex_JointName.at(i);
        joint_state.position[i] = q(i + 7);
      }
      // q(6) = 1;
      q.block(3, 0, 4, 1) = vquad;

      odom_trans.header.stamp = ros::Time::now();
      odom_trans.transform.translation.x = q(0);
      odom_trans.transform.translation.y = q(1);
      odom_trans.transform.translation.z = q(2);
      odom_trans.transform.rotation.x = q(3);
      odom_trans.transform.rotation.y = q(4);
      odom_trans.transform.rotation.z = q(5);
      odom_trans.transform.rotation.w = q(6);

      Broadcaster_odom.sendTransform(odom_trans);

      joint_pub.publish(joint_state);
      // broadcaster.sendTransform(odom_trans);

      std::cout << std::setw(24) << std::left << std::fixed
                << std::setprecision(5);
      pin.forward_kin(q, dq, ddq, false);
      pin.inv_dyn(q, dq, ddq, F, false);

      std::thread thread([&]() {
        if (n.ok()) {
          for (int i = 0; i < 6; i++) {
            try {
              transformStamped = tfBuffer.lookupTransform(
                  "odom", Hex_EEname.at(i), ros::Time(0));
              // std::cout
              //     << Hex_EEname.at(i)
              //     << "Translation: X:" <<
              //     transformStamped.transform.translation.x
              //     << "Y:" << transformStamped.transform.translation.y
              //     << "Z:" << transformStamped.transform.translation.z <<
              //     std::endl;
              Eigen::Vector3d FOOT_POS;
              FOOT_POS(0) = transformStamped.transform.translation.x;
              FOOT_POS(1) = transformStamped.transform.translation.y;
              FOOT_POS(2) = transformStamped.transform.translation.z;
              // foot_x(0, i) = transformStamped.transform.translation.x;
              // foot_x(1, i) = transformStamped.transform.translation.y;
              // foot_x(2, i) = transformStamped.transform.translation.z;
              // pin.forward_kin(q, dq, ddq, false);

              if (pin.x_state.block(3 * i, 0, 3, 1).isApprox(FOOT_POS, 0.02)) {
                std::cout << "pino_fk_right" << std::endl;

                Eigen::Vector3d q_fdb = q.block(3 * i + 7, 0, 3, 1);
                // std::cout << "q_fdb:" << q_fdb.transpose() << std::endl;
                Eigen::Vector3d inv_q =
                    pin.inv_kin_local(FOOT_POS, q_fdb, i, false);
                if (inv_q.isApprox(q.block(3 * i + 7, 0, 3, 1), 0.01)) {
                  std::cout << "ik_fast_right" << std::endl;
                } else {
                  std::cout << "ik_fast_error,pin:" << std::endl;
                }

              } else {
                std::cout << Hex_EEname.at(i) << "pino_error,pin:"
                          << pin.x_state.block(3 * i, 0, 3, 1)
                          << "TF:" << FOOT_POS << std::endl;
                double err =
                    (pin.x_state.block(3 * i, 0, 3, 1) - FOOT_POS).norm();
                std::cout << "error" << err << "\n" << std::endl;
                if (max_error < err) {
                  max_error = err;
                }
                std::cout << "max_error" << max_error << std::endl;
                std::cout << pin.q_state.transpose() << std::endl;
                if (max_error > 2) {
                  exit(0);
                }
              }
            } catch (tf2::TransformException &ex) {
              ROS_WARN("%s", ex.what());
              // ros::Duration(1.0).sleep();
              // continue;
            }
          }
          usleep(1000000);
        }
      });
      thread.detach();

      loop_rate.sleep();
    }
  }
}
