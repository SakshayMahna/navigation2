// Copyright (c) 2024 Sakshay Mahna
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <cmath>
#include "nav2_util/interpolation_utils.hpp"

namespace nav2_util
{
geometry_msgs::msg::Point circleSegmentIntersection(
  const geometry_msgs::msg::Point & p1,
  const geometry_msgs::msg::Point & p2,
  double r)
{
  // Formula for intersection of a line with a circle centered at the origin,
  // modified to always return the point that is on the segment between the two points.
  // https://mathworld.wolfram.com/Circle-LineIntersection.html
  // This works because the poses are transformed into the robot frame.
  // This can be derived from solving the system of equations of a line and a circle
  // which results in something that is just a reformulation of the quadratic formula.
  // Interactive illustration in doc/circle-segment-intersection.ipynb as well as at
  // https://www.desmos.com/calculator/td5cwbuocd
  double x1 = p1.x;
  double x2 = p2.x;
  double y1 = p1.y;
  double y2 = p2.y;

  double dx = x2 - x1;
  double dy = y2 - y1;
  double dr2 = dx * dx + dy * dy;
  double D = x1 * y2 - x2 * y1;

  // Augmentation to only return point within segment
  double d1 = x1 * x1 + y1 * y1;
  double d2 = x2 * x2 + y2 * y2;
  double dd = d2 - d1;

  geometry_msgs::msg::Point p;
  double sqrt_term = std::sqrt(r * r * dr2 - D * D);
  p.x = (D * dy + std::copysign(1.0, dd) * dx * sqrt_term) / dr2;
  p.y = (-D * dx + std::copysign(1.0, dd) * dy * sqrt_term) / dr2;
  return p;
}

geometry_msgs::msg::PoseStamped getLookAheadPoint(
  const double & lookahead_dist,
  const nav_msgs::msg::Path & transformed_plan,
  bool interpolate_after_goal)
{
  // Find the first pose which is at a distance greater than the lookahead distance
  auto goal_pose_it = std::find_if(
    transformed_plan.poses.begin(), transformed_plan.poses.end(), [&](const auto & ps) {
      return hypot(ps.pose.position.x, ps.pose.position.y) >= lookahead_dist;
    });

  // If the no pose is not far enough, take the last pose
  if (goal_pose_it == transformed_plan.poses.end()) {
    if (interpolate_after_goal) {
      auto last_pose_it = std::prev(transformed_plan.poses.end());
      auto prev_last_pose_it = std::prev(last_pose_it);

      double end_path_orientation = atan2(
        last_pose_it->pose.position.y - prev_last_pose_it->pose.position.y,
        last_pose_it->pose.position.x - prev_last_pose_it->pose.position.x);

      // Project the last segment out to guarantee it is beyond the look ahead
      // distance
      auto projected_position = last_pose_it->pose.position;
      projected_position.x += cos(end_path_orientation) * lookahead_dist;
      projected_position.y += sin(end_path_orientation) * lookahead_dist;

      // Use the circle intersection to find the position at the correct look
      // ahead distance
      const auto interpolated_position = circleSegmentIntersection(
        last_pose_it->pose.position, projected_position, lookahead_dist);

      // Calculate orientation towards interpolated position
      // Convert yaw to quaternion
      double interpolated_yaw = atan2(
        interpolated_position.y - last_pose_it->pose.position.y,
        interpolated_position.x - last_pose_it->pose.position.x);

      geometry_msgs::msg::Quaternion interpolated_orientation;
      interpolated_orientation.x = 0.0;
      interpolated_orientation.y = 0.0;
      interpolated_orientation.z = sin(interpolated_yaw / 2.0);
      interpolated_orientation.w = cos(interpolated_yaw / 2.0);

      geometry_msgs::msg::PoseStamped interpolated_pose;
      interpolated_pose.header = last_pose_it->header;
      interpolated_pose.pose.position = interpolated_position;
      interpolated_pose.pose.orientation = interpolated_orientation;

      return interpolated_pose;
    } else {
      goal_pose_it = std::prev(transformed_plan.poses.end());
    }
  } else if (goal_pose_it != transformed_plan.poses.begin()) {
    // Find the point on the line segment between the two poses
    // that is exactly the lookahead distance away from the robot pose (the origin)
    // This can be found with a closed form for the intersection of a segment and a circle
    // Because of the way we did the std::find_if, prev_pose is guaranteed to be inside the circle,
    // and goal_pose is guaranteed to be outside the circle.
    auto prev_pose_it = std::prev(goal_pose_it);
    auto point = circleSegmentIntersection(
      prev_pose_it->pose.position,
      goal_pose_it->pose.position, lookahead_dist);

    // Calculate orientation towards interpolated position
    // Convert yaw to quaternion
    double yaw = atan2(
      point.y - prev_pose_it->pose.position.y,
      point.x - prev_pose_it->pose.position.x);

    geometry_msgs::msg::Quaternion orientation;
    orientation.x = 0.0;
    orientation.y = 0.0;
    orientation.z = sin(yaw / 2.0);
    orientation.w = cos(yaw / 2.0);

    geometry_msgs::msg::PoseStamped pose;
    pose.header.frame_id = prev_pose_it->header.frame_id;
    pose.header.stamp = goal_pose_it->header.stamp;
    pose.pose.position = point;
    pose.pose.orientation = orientation;
    return pose;
  }

  return *goal_pose_it;
}
}  // namespace nav2_util
