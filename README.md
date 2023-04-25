# Parameter estimation from least square minimisation

This code is used to determine the torsional spring and damper parameters of a torsional mass-spring-damper (IBK) system. The parameters are determined by minimising the square differences between the numerical solution of the IBK system and the filtered Motion capture data recorded during a cylindrical grasp. The moment of Inertia parameter is determined uniquely from a cylindrical approximation of the digits which has been verified in a paper that is currently under review. Since the moment of inertia parameter is known the spring and damper parameters can be determined uniquely from the data.


# Digit Constants Minimisation

To run this script, you'll need to specify the sampling frequency of your Motion Capture system (fs) and ensure it's a proper divisor of the Electromyograph (EMG) sampling frequency (fEMG). You'll also need to indicate whether the task performed is within the plane of the gravitational field. If the gravity is perpendicular to the plane of motion, choose grav=[].

To fill the Excel file Moments_of_interia.xlsx, based on the landmarks provided in the "In-vivo hand mass determination and an anthropometric investigation on segment length and radius for prosthetic segment design" paper, select the directory folder where you saved the file and paste it into line 105.

Setting the equilibrium angles manually is crucial. The initial conditions are the starting angle from the motion capture data and the assumed 0 initial velocity before the experiment begins.

All functions are called within this script. Muscle moments for each joint are determined, and these are used as inputs for the IBK system. These functions are determined by fitting polynomial approximations to the data obtained from the OpenSim model found in the "A Musculoskeletal Model of the Hand and Wrist Capable of Simulating Functional Tasks" paper from their Sign Language model. The muscles used in this model are the Extensor Digitorum Communis (EDC), Flexor Digitorum Profundus (FDP), and Flexor Digitorum Superficialis (FDS). The minimization process determines the spring (K) and damper (B) constants alongside scaling parameters for the muscle moment functions. 

To normalize the EMG data from the cylindrical trials, set the maximum voltage readings from the Maximum Voluntary Contraction (MVC) dataset using the MVC_Analysis script. The raw angular data for each finger is filtered using the technique found in Section 3.4.4.3 in "Biomechanics and motor control of human movement." Once the angular data for each joint is filtered, muscle moments are calculated based on the polynomial approximations determined from the OpenSim data mentioned above. The raw EMG data are band-passed filtered between 15 and 450 Hz, as suggested in "Feasibility of using combined EMG and kinematic signals for prosthesis control: A simulation study using a virtual reality environment." Then, the filtered EMG data are rectified, normalized, and the envelopes are obtained using MATLAB's envelope function. A cubic spline interpolation is then used to fit the envelope to obtain the muscle level activation using the formula in "Real-time simulation of hand motion for prosthesis control." The initial condition for the muscle activation is assumed to be zero, and the formula for the activation is taken from the same paper. $T_{act}$ and $T_{deact}$ are taken from the same
paper. The activations are then downsampled so that the activations have the same length as the angular data, and the muscle moments are calculated using the formula found in "Muscle and Tendon: Properties, models, scaling, and application to biomechanics and motor control" by Zajac, where the muscle moment ($\tau_m$) is equal to the product between the muscle level activation and the active muscle force plus the passive force and all that multiplied by the muscle moment arm ($r_m$). All the functions are in general functions of time.

$$ \tau_m=(\alpha_m*F_{active,m}+F_{passive,m})*r_m $$


Once the muscle moments are calculated for each joint, a cubic spline interpolation is used to generate the piecewise polynomial expression of the muscle moments. These functions are used as the actuators in the IBK system. The 2nd ODE is:

$$\rho_1*\tau_{EDC}+\rho_2*\tau_{FDP}+\rho_3*\tau_{FDS}=I_i \ddot \theta_i (t) +B_i \dot \theta_i (t)+ K_i (\theta_i (t)-\theta_{i,eq})$$

Let $p=[K,B,\rho_1,\rho_2,\rho_3]$ be the parameter vector. The parameter vector $p$ is determined by minimising the square differences between the filtered data and the solutions to the IBK model for each joint for each finger. For the DIP joint only the EDC and FDP muscles are used.

Minimise $\sum (\theta_{f,i}-\theta_i)^2$ using the following constraints:

For the MCP joint


$$0.1<=\rho_1<=2$$

$$0.1<=\rho_2<=2$$

$$0.1<=\rho_3<=2$$

$$ 0< B <=2 $$

$$ 0< K <=10 $$

For the PIP joint

$$0.1<=\rho_1<=2$$

$$0.1<=\rho_2<=2$$

$$0.1<=\rho_3<=2$$

$$0< B <=2$$

$$0< K <=10$$

For the DIP joint

$$0.1<=\rho_1<=2$$

$$0.1<=\rho_2<=2$$

$$0< B <=2$$

$$0< K <=5$$


The $\rho_1 ,\rho_2, \rho_3$ are scaling constants for the muscle moments. According to the "A Musculoskeletal Model of the Hand and Wrist Capable of Simulating
Functional Tasks" these muscle moments are representative of the average adult male. Since the experimental procedure was open to everyone we wanted to scale this model based on the collected data up to twice the maximum muscle moment.

The $K, B$ constants are constrained to be always positive and less than twice the values found in "Real-Time simulation of three-dimensional shoulder girdle and arm dynamics" 
since these constants represent the passive moment parameters for the shoulder and arm joints.

# Structural Identifiability

Let $p=[K,B,\rho_1,\rho_2,\rho_3]$ be the parameter vector and the 2nd ODE solved is:

$$\rho_1*\tau_{EDC}+\rho_2*\tau_{FDP}+\rho_3*\tau_{FDS}=I_i \ddot \theta_i (t) +B_i \dot \theta_i (t)+ K_i (\theta_i (t)-\theta_{i,eq})$$

Taking the Laplace Transform of the above and solving for \theta yields.

$$\theta (s) = \rho_1 / I_{i} *\tau_{EDC} (s)+\rho_2 / I_{i} *\tau_{FDP} (s)+ \rho_3 / I_{i} *\tau_{FDS} (s) -\theta (0)*(s+B_i/I_i)/(s^2+B_i s/I_i+K_i/I_i)+K_i*\theta_{i,eq}/(I_i*(s^2+B_i s/I_i+K_i/I_i)$$

From the above the muscle moments are known functions and the parameters that are determined uniquely are:

$\rho_1/I_i,\rho_2/I_i,\rho_3/I_i,\theta (0),\theta_{i,eq}, K_i/I_i, B_i/I_i$

This model is unidentifiable and an a-priori knowledge of at least one parameter is needed to determine the rest. Because the moment of inertia $I_i$ is known from the cylindrical approximation of the segments, hence it is unique then the parameter vector $p$ is determined uniquely from the input-output relationship and our model is structurally identifiable. 

# User defined parameters

In this code you will be prompted to input different numbers based on the plots provided.

Set the sampling frequency of the MoCap at fs.

Set the sampling frequency of the EMG at the EMG_2_Fit function at fEMG.

Set the aEDCMax, aFDPMax, and aFDSMax from a MVC trial after using the MVC_Analysis script to determine these values.

Set the aEDC, aFDP and aFDS matrices from the raw EMG data during the cylindrical grasp.

Set th1, th2 and th3 the raw angular data for the MCP, PIP and DIP joint movements of the finger you will be analysing.

Set the theq vector that will hold the equilibrium angles of each joint determined from a static capture.

When you run the main script:

Initially you will be asked to determine the linear part of the cut-off frequency for each joint in order to fit a linear model so that the residual at 0 Hz can be determined,  as shown in section 3.4.4.3 in "Biomechanics and motor control of human movement" in Figure 3.20. The cut-off frequency is then determined automatically and used in the creation of a 4th order low pass Butterworth filter. You will be asked if your MCP joint data show a sharp edge at 0 degrees. If this is the case press Y or y. This sharp edge corresponds to the extension angles of the MCP joint and the OpenSim convention is that these angles are negative. Because the way angles are calculated in ProCalc from the developed model both flexion and extension angles are positive. However, a visual inspection on the C3D file in Vicon nexus can give you an idea on which part of the signal is the extension. In this case there were participants that extended their fingers well above 0 degrees, before they performed the cylindrical grasp which is a typical finger flexion movement. That's why the sharp edge at zero degrees is visible. If not press any other button. 

Next the participant number and finger will be asked. Fingers are labelled as numbers and correspond to 2 (index), 3 (middle), 4 (ring), 5 (little). These numbers will be used to obtain the moments of inertia of the respective finger segments from the Moments of Inertia excel file. Furthermore, the respective muscle forces and moment arm functions are going to be used for the specific finger motion. One finger motion (i.e. Index MCP, Index PIP and Index DIP) data are used at a time. 

Then the peak number for determining the envelope of the filtered and rectified EMG data for each muscle is going to be asked. This is a number typically taken between 5-12. This is a visual representation of the actual envelope. Choose the peak that encloses the EMG signals the best.

You will then be asked to press 1 to continue. A set of plots will be shown that will give you information about the muscle moments, activations with different methods and the envelopes of the EMG signals.

Lastly, a prompt will be given in order to save the muscle moment data, time, filtered angles, and IBK solutions with optimised values in order to be used with the Mathematica code to visualise the response of the Lagrangian model and compare the results.





