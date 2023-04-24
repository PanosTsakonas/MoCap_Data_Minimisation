# Parameter estimation from least square minimisation

This code is used to determine the torsional spring and damper parameters of a torsional mass-spring-damper (IBK) system. The parameters are determined by minimising the square differences between the numerical solution of the IBK system and the filtered Motion capture data recorded during a cylindrical grasp. The moment of Inertia parameter is determined uniquely from a cylindrical approximation of the digits which has been verified in a paper that is currently under review. Since the moment of inertia parameter is known the spring and damper parameters can be determined uniquely from the data.


# Digit Constants Minimisation

This is the main script that you must run. You need to specify the sampling frequency of your Motion Capture system (fs), and it must be an proper divisor of the Electromyograph (EMG) sampling frequency, and whether the task performed is within the plane of the gravitational field. The variable grav=[] is chosen when the gravity is perpendicular to the plane of motion. You will need to fill the excel file Moments_of_interia.xlsx based on the landmarks provided in the "In-vivo hand mass determination and an anthropometric investigation on segment length and radius for prosthetic segment design" paper. You need to select the directory folder where you have saved this file and paste in line 105. Setting the equilibrium angles is important and must be done manually. Initial conditions are the starting angle from the motion capture data and assumed 0 initial velocity before the experiment starts. Within this script all the functions are called, where the muscle moments for each joint are determined and are used as inputs for the IBK system. The parameters determined from the minimisation process are the spring and damper constants alongside scaling parameters for the muscle moment functions. These function are determined by fitting polynomial approximations to the data obtained from the OpenSim model found in the 
"A Musculoskeletal Model of the Hand and Wrist Capable of Simulating
Functional Tasks" paper from their Sign Language model. The muscles that are used in this model are the Extensor Digitorum Communis (EDC), Flexor Digitorum Profundus (FDP) and Flexor DIgitorum Superficialis (FDS). It is quite important to set the maximum voltage readings from the Maximum Voluntary Contraction (MVC) dataset in order to normalise the EMG data from the cylindrical trials. Within this script the raw angular data for each finger are filtered using the technique found in Section 3.4.4.3 in Biomechanics and motor control of human movement. Once the angular data for each joint are filtered ($\theta_{f,i}$) then the muscle moments are calculated based on the polynomial approximations determined from the OpenSim data mentioned above. Then the EMG data are band-passed filtered between 15 and 450 Hz as suggested in "Feasibility of using 
combined EMG and kinematic signals for prosthesis control: A simulation study 
using a virtual reality environment". Then the filtered EMG data are rectified, normalised and the envelope is obtained using MATLAB's envelope function. Once the envelope is obtained a cubic spline interpolation is used to fit the latter in order to obtain the muscle level activation ($\alpha$) using the formula in "Real-time simulation of hand motion for prosthesis control". The initial condition for the 
muscle activation, since there is no motion before the experiment, is
assumed to be equal to zero. $T_{act}$ and $T_{deact}$ are taken from the same
paper. Following this calculation, the activations are then downsampled so that the activations have the same length as the angular data and the muscle moments are calculated using the formula found in "Muscle and Tendon: Propoerties, models, scaling, and application to biomechanics and motor control" by Zajac, where the muscle moment ($\tau_m$) is equal to the product of the muscle level activation and the active muscle force plus the passive force and all that multiplied by the muscle moment arm ($r_m$). All the functions are in general functions of time.

$$ \tau_m=(\alpha_m*F_{active,m}+F_{passive,m})*r_m $$


Once the muscle moments are calculated for each joint then a cubic spline interpolation is to generate the piecewise polynomial expression of the muscle moments. Then these functions are used as the actuators in the IBK system. The optimisation procedure is as follows: Let $p=[B,K,\rho_1,\rho_2,\rho_3]$ be the parameter vector that we need to determine. The 2nd ODE is:

$$\rho_1*\tau_{EDC}+\rho_2*\tau_{FDP}+\rho_3*\tau_{FDS}=I_i \ddot \theta_i (t) +B_i \dot \theta_i (t)+ K_i (\theta_i (t)-\theta_{i,eq})$$

The parameter vector $p$ is determined by minimising the square differences between the filtered data and the solution to the IBK model for each joint for each finger

Minimise $\sum (\theta_{f,i}-\theta_i)^2$ using the following constraints:

For the MCP joint


$$0.1<=\rho_1<=2$$

$$0.1<=\rho_2<=2$$

$$0.1<=\rho_3<=2$$

$$ 0< B <=2 $$

$$ 0< K <=10 $$

For the PIP joint

$$0.1<=\rho_1<=1.8$$

$$0.1<=\rho_2<=1.8$$

$$0.1<=\rho_3<=1.8$$

$$0< B <=2$$

$$0< K <=10$$

For the DIP joint

$$0.1<=\rho_1<=1.8$$

$$0.1<=\rho_2<=1.8$$

$$0< B <=2$$

$$0< K <=5$$

