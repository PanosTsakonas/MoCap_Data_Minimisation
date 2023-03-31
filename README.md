# Parameter estimation from least square minimisation

This code is used to determine the torsional spring and damper parameters of a torsional mass-spring-damper (IBK) system. The parameters are determined by minimising the square differences between the numerical solution of the IBK system and the filtered Motion capture data recorded during a cylindrical grasp. The moment of Inertia parameter is determined uniquely from a cylindrical approximation of the digits which has been verified in a paper that is current under review. Since the moment of inertia parameter is known the spring and damper parameters can be determined uniquely from the data.


#Digit Constants Minimisation

This is the main script that you must run. You need to specify the sampling frequency of your Motion Capture system (fs), and whether the task performed is within the plane of the gravitational field. The variable grav=[] is chosen when the gravity is perpendicular to the plane of motion. You will need to fill the excel file Moments_of_interia.xlsx based on the landmarks provided in the "In-vivo hand mass determination and an anthropometric investigation on segment length and radius for prosthetic segment design" paper. You need to select the directory folder where you have saved this file and paste in line 84. 
