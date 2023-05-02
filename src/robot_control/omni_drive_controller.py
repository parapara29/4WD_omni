#!/usr/bin/env python3

import math
import rospy
import numpy as np
from geometry_msgs.msg import Twist
from std_msgs.msg import Float64

class OmniDriveController:
    def __init__(self):
        rospy.init_node('omni_drive_controller')

        # Wheel controller publishers
        self.wheel1_pub = rospy.Publisher('/robot/Revolute_34_controller/command', Float64, queue_size=1)
        self.wheel2_pub = rospy.Publisher('/robot/Revolute_42_controller/command', Float64, queue_size=1)
        self.wheel3_pub = rospy.Publisher('/robot/Revolute_43_controller/command', Float64, queue_size=1)
        self.wheel4_pub = rospy.Publisher('/robot/Revolute_44_controller/command', Float64, queue_size=1)

        # Subscribe to the cmd_vel topic
        rospy.Subscriber('/robot/cmd_vel', Twist, self.cmd_vel_callback)

    def cmd_vel_callback(self, msg):
        vx = msg.linear.x
        vy = msg.linear.y
        omega = msg.angular.z

        r = 0.21  # Wheel radius in meters
        R = 0.2   # Distance from the center of the robot to the wheel in meters

        A = math.sin(omega + math.pi / 4)
        B = math.cos(omega + math.pi / 4)
        C = math.sin(omega + 3 * math.pi / 4)
        D = math.cos(omega + 3 * math.pi / 4)
        E = math.sin(omega + 5 * math.pi / 4)
        F = math.cos(omega + 5 * math.pi / 4)
        G = math.sin(omega + 7 * math.pi / 4)
        H = math.cos(omega + 7 * math.pi / 4)

        Vw = (1 / r) * np.array([
            [-A, B, R],
            [-C, D, R],
            [-E, F, R],
            [-G, H, R]
        ])

        wheel_velocities = Vw.dot(np.array([vx, vy, omega]))

        # Publish the calculated wheel velocities
        self.wheel1_pub.publish(Float64(wheel_velocities[0]))
        self.wheel2_pub.publish(Float64(wheel_velocities[1]))
        self.wheel3_pub.publish(Float64(wheel_velocities[2]))
        self.wheel4_pub.publish(Float64(wheel_velocities[3]))

if __name__ == '__main__':
    try:
        omni_drive_controller = OmniDriveController()
        rospy.spin()
    except rospy.ROSInterruptException:
        pass
