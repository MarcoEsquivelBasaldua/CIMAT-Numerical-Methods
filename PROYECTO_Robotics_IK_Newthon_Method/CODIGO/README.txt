Code name: Numerical_IK.cpp
Author: Marco Antonio Esquivel Basaldua
Date: 18/12/19

Code description:
    Given an 3 DoF articulated robotic arm in three dimensions which lenght of each link is specified in the
    file "lengths.txt", this code determines the set of nedded angles to each link so the end
    effector of the robot will pass throught the points specified in the file "Path.txt".
    The angles generated can be used to simulate the robotic arm in the software MATLAB.
    
Inputs:
- "lengths.txt" lengths of eacho of three links
- "Path.txt" each of desired point (in coordinates x y z, separated by a space) to be visited by the end-efector. 
    The first value in this file must state the number of desired points to be visited.
    
Outputs:
- "angles.txt": set of calculated angles to be used on every articulation to pass throught the desired points

Compiler instructions:
- By positioning in the folder where "Numerical_IK.cpp", "Path.txt" and "lengths.txt" are stored, use the Makefile placed in this folder, just write "make" and press enter in the terminal       

Visualisation instructions:
- Once "Numerical_IK.cpp" excecuted and the file "angles.txt" generated, run on MATLAB the file "Plot_IK.m". Results will be ploted.
