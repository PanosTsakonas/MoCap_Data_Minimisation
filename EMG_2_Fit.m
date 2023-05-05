function [SF_MCP,SF_PIP,SF_DIP,tEMG,Res,aEDCf,aFDPf,aFDSf]=EMG_2_Fit(th1f,th2f,th3f,aEDC,aFDP,aFDS,aEDCMax,aFDPMax,aFDSMax,tim,fs,in)
global mn1 mn2 mn3;
syms x
%Set the sampling frequency of the EMG
fsEMG=1000;
%Nyquist Frequency
fnEMG=fsEMG/2;

%Setting the time vector for the EMG signals
tEMG=zeros(length(aEDC),1);
for i=2:length(aEDC)
    tEMG(i)=tEMG(i-1)+1/fsEMG;
end

%For the finger of interest calculate the moment arm active and passive
%forces per degree of freedom. These functions will be used in the
%calculation of the muscle moments. These function are determined by
%fitting polynomial degrees to the data obtained from the OpenSim model of
%"A Musculoskeletal Model of the Hand and Wrist Capable of Simulating
%Functional Tasks" and the data came from the Sign Language model.
 
if in==2
    disp("You are fitting the Index finger muscle data");
    pause(0.5);
%Index Finger
%MCP
%Moment arm
rEDC_MCP=feval(matlabFunction(0.001608*x^4-0.005081*x^3+0.004694*x^2+0.001651-0.01035),th1f);
rFDP_MCP=feval(matlabFunction(-7.611*10^-5*x^4-0.0002441*x^3+0.0004334*x^2+0.0007466*x+0.009723),th1f);
rFDS_MCP=feval(matlabFunction(-0.0003308*x^3+0.0003252*x^2+0.00115*x+0.01081),th1f);
%Active Force
FA_EDC_MCP=feval(matlabFunction(0.5357*x^4-1.25*x^3-0.8672*x^2+3.179*x+47.27),th1f);
FA_FDP_MCP=feval(matlabFunction(-2.958*x^4+9.664*x^3-13.5*x^2-1.132*x+198.2),th1f);
FA_FDS_MCP=feval(matlabFunction(-1.98*x^4+10.4*x^3-22.18*x^2+12.26*x+161.2),th1f);
%Passive Force
PA_EDC_MCP=feval(matlabFunction(0.047*x^4+.001064*x^3-0.03169*x^2-0.003156*x+0.002138),th1f);
PA_FDP_MCP=feval(matlabFunction(0.1201*x^4-0.6952*x^3+1.182*x^2-0.6844*x+0.09092),th1f);
PA_FDS_MCP=feval(matlabFunction(-0.129*x^4-0.2022*x^3+1.446*x^2-1.604*x+0.4801),th1f);

%PIP
%Moment arm
rEDC_PIP=feval(matlabFunction(0.0002685*x^4-0.0003149*x^3-0.001513*x^2+0.003234*x-0.005506),th2f);
rFDP_PIP=feval(matlabFunction(-0.001802*x^4+0.004753*x^3-0.00347*x^2+0.00154*x+0.00712),th2f);
rFDS_PIP=feval(matlabFunction(-0.001694*x^2+0.001553*x+0.007548),th2f);
%Active Force
FA_EDC_PIP=feval(matlabFunction(-0.4909*x^2+1.796*x+47.19),th2f);
FA_FDP_PIP=feval(matlabFunction(2.828*x^3-11*x^2+4.239*x+197),th2f);
FA_FDS_PIP=feval(matlabFunction(1.229*x^4-4.755*x^3+3.314*x^2+0.2573*x+162.1),th2f);
%Passive Force
PA_EDC_PIP=zeros(length(th2f),1);
PA_FDP_PIP=zeros(length(th2f),1);
PA_FDS_PIP=feval(matlabFunction(0.2426*x^4-1.3*x^3+2.495*x^2-2.01*x+0.5659),th2f);

%DIP
%Moment arm
rEDC_DIP=feval(matlabFunction(0.0001762*x^2+0.001466*x-0.003962),th3f);
rFDP_DIP=feval(matlabFunction(-0.0006361*x^2+0.002305*x+0.003045),th3f);
%Active Force
FA_EDC_DIP=feval(matlabFunction(-0.2784*x^2+1.223*x+47.22),th3f);
FA_FDP_DIP=feval(matlabFunction(-2.014*x^3+1.788*x^2-0.4118*x+197.3),th3f);
%Passive Force
PA_EDC_DIP=zeros(length(th3f),1);
PA_FDP_DIP=zeros(length(th3f),1);

elseif in==3
    disp("You are fitting the Middle finger muscle data");
    pause(0.5);
%Middle
%MCP
%Moment arm
rEDC_MCP=feval(matlabFunction(-0.002325*x^3+0.007136*x^2-0.003553*x-0.007848),th1f);
rFDP_MCP=feval(matlabFunction(-0.0002092*x^4-0.002243*x^3+0.002023*x^2+0.004282*x+0.008124),th1f);
rFDS_MCP=feval(matlabFunction(0.0008655*x^4-0.002502*x^3+0.0002376*x^2+0.005472*x+0.009507),th1f);
%Active Force
FA_EDC_MCP=feval(matlabFunction(3.864*x^4-7.771*x^3-6.391*x^2+1.513*x+94.91),th1f);
FA_FDP_MCP=feval(matlabFunction(1.72*x^4-2.342*x^3-6.651*x^2+0.1405*x+212.8),th1f);
FA_FDS_MCP=feval(matlabFunction(0.704*x^4+1.456*x^3-13.49*x^2-1.198*x+259.5),th1f);
%Passive Force
PA_EDC_MCP=feval(matlabFunction(-0.09657*x^4+0.1466*x^3+0.5417*x^2+0.3078*x+0.03137),th1f);
PA_FDP_MCP=feval(matlabFunction(0.03674*x^4-0.3602*x^3+0.731*x^2-0.4595*x+0.06231),th1f);
PA_FDS_MCP=feval(matlabFunction(-0.002824*x^4-0.4485*x^3+1.154*x^2-0.8103*x+0.1237),th1f);


%PIP
%Moment arm
rEDC_PIP=feval(matlabFunction(0.001683*x^4-0.006261*x^3+0.006067*x^2+0.0001409*x-0.00531),th2f);
rFDP_PIP=feval(matlabFunction(-0.002828*x^4+0.007248*x^3-0.007375*x^2+0.006607*x+0.006253),th2f);
rFDS_PIP=feval(matlabFunction(-0.01151*x^6+0.05281*x^5-0.09135*x^4+0.07332*x^3-0.0282*x^2+0.00846*x+0.006936),th2f);
%Active Force
FA_EDC_PIP=feval(matlabFunction(-2.425*x^3+2.22*x^2-0.7452*x+94.45),th2f);
FA_FDP_PIP=feval(matlabFunction(1.6*x^4-2.026*x^3-7.438*x^2+3.635*x+212.1),th2f);
FA_FDS_PIP=feval(matlabFunction(-0.3559*x^4+6.384*x^3-20.57*x^2+7.344*x+258.4),th2f);
%Passive Force
PA_EDC_PIP=feval(matlabFunction(0.08758*x^3-0.1471*x^2+0.4855*x-0.03225),th2f);
PA_FDP_PIP=zeros(length(th2f),1);
PA_FDS_PIP=zeros(length(th2f),1);


%DIP 
%Moment arm
rEDC_DIP=feval(matlabFunction(0.00131*x^4-0.003982*x^3+0.003192*x^2+0.0002258*x-0.003475),th3f);
rFDP_DIP=feval(matlabFunction(-0.0007404*x^3+0.0001684*x^2+0.002696*x+0.004166),th3f);
%Active Force
FA_EDC_DIP=feval(matlabFunction(-0.4358*x^3+0.3924*x^2-0.1479*x+94.41),th3f);
FA_FDP_DIP=feval(matlabFunction(3.89*x^4-12.32*x^3+8.914*x^2-2.035*x+212.5),th3f);
%Passive Force
PA_EDC_DIP=feval(matlabFunction(0.1872*x^4-0.5826*x^3+0.6213*x^2+0.02121*x-0.003635),th3f);
PA_FDP_DIP=zeros(length(th3f),1);

elseif in==4
    disp("You are fitting the Ring finger muscle data");
    pause(0.5);
%Ring
%MCP
%Moment arm
rEDC_MCP=feval(matlabFunction(-0.0006839*x^4-0.001755*x^3+0.006419*x^2-0.003621*x-0.005524),th1f);
rFDP_MCP=feval(matlabFunction(0.0001133*x^5-0.001472*x^4-0.0009706*x^3+0.002839*x^2+0.002742*x+0.008508),th1f);
rFDS_MCP=feval(matlabFunction(-5.466*10^-5*x^5-0.0007669*x^4-0.001297*x^3+0.002332*x^2+0.003638*x+0.008996),th1f);
%Active Force
FA_EDC_MCP=feval(matlabFunction(-3.868*x^3-2.812*x^2+0.4682*x+109.5),th1f);
FA_FDP_MCP=feval(matlabFunction(0.6106*x^4+0.4321*x^3-7.098*x^2-0.1922*x+173.4),th1f);
FA_FDS_MCP=feval(matlabFunction(-0.1076*x^4+2.636*x^3-9.668*x^2-0.8327*x+171.8),th1f);
%Passive Force
PA_EDC_MCP=feval(matlabFunction(0.06846*x^3+0.4289*x^2+0.3043*x+0.03531),th1f);
PA_FDP_MCP=feval(matlabFunction(0.05169*x^4-0.3896*x^3+0.7339*x^2-0.4454*x+0.0592),th1f);
PA_FDS_MCP=feval(matlabFunction(0.06219*x^4-0.4495*x^3+0.8366*x^2-0.5069*x+0.06835),th1f);

%PIP 
%Moment arm
rEDC_PIP=feval(matlabFunction(-0.003131*x^6+0.0178*x^5-0.03703*x^4+0.03306*x^3-0.01186*x^2+0.002918*x-0.004137),th2f);
rFDP_PIP=feval(matlabFunction(-0.001625*x^4+0.004082*x^3-0.004308*x^2+0.004603*x+0.004498),th2f);
rFDS_PIP=feval(matlabFunction(-0.0003201*x^4+0.0002555*x^3-0.001011*x^2+0.003327*x+0.00302),th2f);
%Active Force
FA_EDC_PIP=feval(matlabFunction(-1.895*x^3+2.647*x^2-1.304*x+109.3),th2f);
FA_FDP_PIP=feval(matlabFunction(2.478*x^4-7.338*x^3+3.05*x^2-0.03967*x+172.9),th2f);
FA_FDS_PIP=feval(matlabFunction(2.353*x^4-7.979*x^3+6.121*x^2-1.455*x+171.3),th2f);
%Passive Force
PA_EDC_PIP=feval(matlabFunction(0.1665*x^4-0.5094*x^3+0.4837*x^2+0.2482*x-0.01589),th2f);
PA_FDP_PIP=zeros(length(th2f),1);
PA_FDS_PIP=zeros(length(th2f),1);

%DIP
%Moment arm
rEDC_DIP=feval(matlabFunction(0.000406*x^2+0.001107*x-0.003328),th3f);
rFDP_DIP=feval(matlabFunction(-0.0006442*x^2+0.002374*x+0.00279),th3f);
%Active Force
FA_EDC_DIP=feval(matlabFunction(-0.3268*x^2+0.1253*x+109.2),th3f);
FA_FDP_DIP=feval(matlabFunction(-0.4454*x^4-0.3845*x^3+0.7531*x^2-0.2591*x+172.9),th3f);
%Passive Force
PA_EDC_DIP=feval(matlabFunction(0.2123*x^4-0.7213*x^3+0.7467*x^2+0.05183*x-0.005412),th3f);
PA_FDP_DIP=zeros(length(th3f),1);

else
    disp("You are fitting the Little finger muscle data");
    pause(0.5);
%Little
%MCP
%Moment arm
rEDC_MCP=feval(matlabFunction(0.001326*x^4-0.003376*x^3+0.001346*x^2+0.001553*x-0.006571),th1f);
rFDP_MCP=feval(matlabFunction(-0.002*x^4+0.003392*x^3-0.0001561*x^2+0.0004944*x+0.006389),th1f);
rFDS_MCP=feval(matlabFunction(0.0007292*x^4-0.003182*x^3+0.002078*x^2+0.004899*x+0.006576),th1f);
%Active Force
FA_EDC_MCP=feval(matlabFunction(0.1615*x^4-1.588*x^3-1.332*x^2+0.4985*x+39.52),th1f);
FA_FDP_MCP=feval(matlabFunction(0.996*x^4-2.225*x^3-4.049*x^2+0.7381*x+237.1),th1f);
FA_FDS_MCP=feval(matlabFunction(-0.1581*x^4+0.9791*x^3-4.369*x^2-0.6144*x+75.54),th1f);
%Passive Force
PA_EDC_MCP=feval(matlabFunction(0.01449*x^4+0.01199*x^3+0.1422*x^2+0.1139*x+0.01443),th1f);
PA_FDP_MCP=feval(matlabFunction(-0.01485*x^4-0.259*x^3+0.7214*x^2-0.509*x+0.0734),th1f);
PA_FDS_MCP=feval(matlabFunction(0.02651*x^4-0.1897*x^3+0.3542*x^2-0.2174*x+0.03071),th1f);

%PIP 
%Moment arm
rEDC_PIP=feval(matlabFunction(0.00277*x^4-0.01016*x^3+0.009658*x^2+7.733*10^-5*x-0.005168),th2f);
rFDP_PIP=feval(matlabFunction(-0.008171*x^4+0.02129*x^3-0.01926*x^2+0.008913*x+0.005823),th2f);
rFDS_PIP=feval(matlabFunction(-0.001302*x^4+0.003318*x^3-0.003714*x^2+0.004776*x+0.004874),th2f);
%Active Force
FA_EDC_PIP=feval(matlabFunction(-0.3888*x^4+0.2797*x^3-0.2101*x^2+0.01148*x+39.4),th2f);
FA_FDP_PIP=feval(matlabFunction(2.84*x^4-4.608*x^3-6.366*x^2+3.374*x+236.6),th2f);
FA_FDS_PIP=feval(matlabFunction(1.529*x^3-5.794*x^2+1.927*x+75.2),th2f);
%Passive Force
PA_EDC_PIP=feval(matlabFunction(0.0055*x^3-0.1218*x^2+0.2407*x-0.01559),th2f);
PA_FDP_PIP=zeros(length(th2f),1);
PA_FDS_PIP=zeros(length(th2f),1);

%DIP Flexion
%Moment arm
rEDC_DIP=feval(matlabFunction(0.0002091*x^2+0.002192*x-0.002416),th3f);
rFDP_DIP=feval(matlabFunction(-0.0007132*x^3-0.000198*x^2+0.002966*x+0.003281),th3f);
%Active Force
FA_EDC_DIP=feval(matlabFunction(-0.003527*x^4+0.0209*x^3-0.02745*x^2+0.005582*x+39.4),th3f);
FA_FDP_DIP=feval(matlabFunction(-2.033*x^3+0.5217*x^2+0.3778*x+236.7),th3f);
%Passive Force
PA_EDC_DIP=feval(matlabFunction(0.04959*x^4-0.1804*x^3+0.1769*x^2-0.01252*x+7.734*10^-5),th3f);
PA_FDP_DIP=zeros(length(th3f),1);
end

%The next thing to be imported will be the EMG data. According to "Feasibility of using 
%combined EMG and kinematic signals for prosthesis control: A simulation study 
%using a virtual reality environment" the EMG signals were ampli-
%fied, band-pass filtered between 15 and 450 Hz. This is also supported by
%plotting the Hilbert-Huang transform of the EMG signals


%The band pass filter

[bemg,aemg]=butter(4,[15 450]./fnEMG,'bandpass');
aEDCf=filtfilt(bemg,aemg,aEDC);
aFDPf=filtfilt(bemg,aemg,aFDP);
aFDSf=filtfilt(bemg,aemg,aFDS);
clear bemg aemg
%Removing the 50Hz power line frequency
[bemg,aemg]=butter(4,[49 51]./fnEMG,'stop');
aEDCf=filtfilt(bemg,aemg,aEDCf);
aFDPf=filtfilt(bemg,aemg,aFDPf);
aFDSf=filtfilt(bemg,aemg,aFDSf);
clear bemg aemg

%Rectify by taking the absolute of the EMG data as suggested
%in "Finger joint coordination during tapping"
aEDCf=abs(aEDCf);
aFDPf=abs(aFDPf);
aFDSf=abs(aFDSf);

%Normalize
if isempty(aEDCMax)==1
aEDCf=aEDCf./max(aEDCf);
else
    aEDCf=aEDCf./aEDCMax;
end

if isempty(aFDPMax)==1
aFDPf=aFDPf./max(aFDPf);
else
    aFDPf=aFDPf./aFDPMax;
end

if isempty(aFDSMax)==1
aFDSf=aFDSf./max(aFDSf);
else
    aFDSf=aFDSf./aFDSMax;
end

%obtain the envelope by low pass filtering the rectified data
[bemg,aemg]=butter(4,10/fnEMG,'low');
yedc1=filtfilt(bemg,aemg,aEDCf);
yfdp1=filtfilt(bemg,aemg,aFDPf);
yfds1=filtfilt(bemg,aemg,aFDSf);

%The next two steps are taken from "A fast implementation for EMG signal 
%linear envelope computation".
%{
%Obtain the smoothed rectified signal from moving average filter for 5
%data points
yedc1=smooth(aEDCf,5);
yfdp1=smooth(aFDPf,5);
yfds1=smooth(aFDSf,5);

%Low pass at 30 Hz
[bemg,aemg]=butter(4,30/fnEMG,'low');
yedc1=filtfilt(bemg,aemg,yedc1);
yfdp1=filtfilt(bemg,aemg,yfdp1);
yfds1=filtfilt(bemg,aemg,yfds1);
%}


if (in==2)
ll=1;
mn1=10;
%Obtain the envelope of the normalised original EMG signal
[yedc,~]=envelope(aEDCf,mn1,'peak');
figure
while(ll==1)
plot(tEMG,aEDCf,tim,downsample(yedc,ceil(fsEMG/fs)));
legend("EMG Raw Data","EMG envelope");
xlabel("Time (s)");
ylabel("Normalised EMG");
title("Raw EMG data and envelope for EDC muscle"); 
ll=input("Is the plot fine? Select 1 to repeat any other number to continue ");

if (ll==1)
    mn1=input("Give a different peak number: ");
    [yedc,~]=envelope(aEDCf,mn1,'peak');
    close all
end
end
mn2=5;
[yfdp,~]=envelope(aFDPf,mn2,'peak');
figure
ll=1;
while(ll==1)
plot(tEMG,aFDPf,tim,downsample(yfdp,ceil(fsEMG/fs)));
legend("EMG Raw Data","EMG envelope");
xlabel("Time (s)");
ylabel("Normalised EMG");
title("Raw EMG data and envelope for FDP muscle");
ll=input("Is the plot fine? Select 1 to repeat any other number to continue ");

if (ll==1)
    mn2=input("Give a different peak number: ");
    [yfdp,~]=envelope(aFDPf,mn2,'peak');
    close all
end
end
mn3=6;
[yfds,~]=envelope(aFDSf,mn3,'peak');
figure
ll=1;
while(ll==1)
    
plot(tEMG,aFDSf,tim,downsample(yfds,ceil(fsEMG/fs)));
legend("EMG Raw Data","EMG envelope");
xlabel("Time (s)");
ylabel("Normalised EMG");
title("Raw EMG data and envelope for FDS muscle");
ll=input("Is the plot fine? Select 1 to repeat any other number to continue ");

if (ll==1)
    mn3=input("Give a different peak number: ");
    [yfds,~]=envelope(aFDSf,mn3,'peak');
    close all
end
end

else
    [yedc,~]=envelope(aEDCf,mn1,'peak');
    [yfdp,~]=envelope(aFDPf,mn2,'peak');
    [yfds,~]=envelope(aFDSf,mn3,'peak');
end

%Use the activation dynamics differential equation to obtain the
%activations from integration using ode45 and spline interpolation of the
%normalised muscle excitations
figure
plot(tEMG,aEDCf,tEMG,yedc);
legend("EMG Raw Data","EMG envelope");
xlabel("Time (s)");
ylabel("Normalised EMG");
title("Raw EMG data and envelope for EDC muscle");
figure
plot(tEMG,aFDPf,tEMG,yfdp);
legend("EMG Raw Data","EMG envelope");
xlabel("Time (s)");
ylabel("Normalised EMG");
title("Raw EMG data and envelope for FDP muscle");
figure
plot(tEMG,aFDSf,tEMG,yfds);
legend("EMG Raw Data","EMG envelope");
xlabel("Time (s)");
ylabel("Normalised EMG");
title("Raw EMG data and envelope for FDS muscle");

[b3,a3]=butter(4,10/fnEMG,'low');
[b4,a4]=butter(4,80/fnEMG,'low');
figure
plot(tEMG,filtfilt(b3,a3,aEDCf),tEMG,yedc,tEMG,filtfilt(b4,a4,aEDCf));
legend("EMG Linear envelope at 10 Hz","EMG envelope","EMG Linear envelope at 80 Hz");
xlabel("Time (s)");
ylabel("Normalised EMG");
title("EMG envelopes for EDC muscle");
figure
plot(tEMG,filtfilt(b3,a3,aFDPf),tEMG,yfdp,tEMG,filtfilt(b4,a4,aFDPf));
legend("EMG Linear envelope at 10 Hz","EMG envelope","EMG Linear envelope at 80 Hz");
xlabel("Time (s)");
ylabel("Normalised EMG");
title("EMG envelopes for FDP muscle");
figure
plot(tEMG,filtfilt(b3,a3,aFDSf),tEMG,yfds,tEMG,filtfilt(b4,a4,aFDSf));
legend("EMG Linear envelope at 10 Hz","EMG envelope","EMG Linear envelope at 80 Hz");
xlabel("Time (s)");
ylabel("Normalised EMG");
title("EMG envelopes for FDS muscle");


pedc=spline(tEMG,yedc);
pfdp=spline(tEMG,yfdp);
pfds=spline(tEMG,yfds);

ppedc=mkpp(pedc.breaks,pedc.coefs);
ppfdp=mkpp(pfdp.breaks,pfdp.coefs);
ppfds=mkpp(pfds.breaks,pfds.coefs);

%The following are taken from "Real-time simulation of hand motion for prosthesis control". The initial
%muscle activation since there is no motion before the experiment is
%assumed to be equal to zero. Tact and Tdeact are taken from the same
%paper.

Tact=15*10^-3;
Tdeact=50*10^-3;


[~,a_edc_ode]=ode45(@(t,a) activation(t,a,ppedc),tEMG,0);
[~,a_fdp_ode]=ode45(@(t,a) activation(t,a,ppfdp),tEMG,0);
[~,a_fds_ode]=ode45(@(t,a) activation(t,a,ppfds),tEMG,0);

%Downsample the raw EMG data.
yedc=downsample(aEDCf,ceil(fsEMG/fs));
yfdp=downsample(aFDPf,ceil(fsEMG/fs));
yfds=downsample(aFDSf,ceil(fsEMG/fs));

%Itterative method for the muscle level activation following the work of
%"Real-time simulation of hand motion forprosthesis control". The initial
%muscle activation since there is no motion before the experiment is
%assumed to be equal to zero.

ac_edc=zeros(length(yedc),1);
ac_fdp=zeros(length(yedc),1);
ac_fds=zeros(length(yedc),1);

for i=1:length(yedc)-1
    ac_edc(i+1)=ac_edc(i)+(tim(i+1)-tim(i))*(yedc(i)/Tact + (1-yedc(i))/Tdeact)*(yedc(i)-ac_edc(i));
    ac_fdp(i+1)=ac_fdp(i)+(tim(i+1)-tim(i))*(yfdp(i)/Tact + (1-yfdp(i))/Tdeact)*(yfdp(i)-ac_fdp(i));
    ac_fds(i+1)=ac_fds(i)+(tim(i+1)-tim(i))*(yfds(i)/Tact + (1-yfds(i))/Tdeact)*(yfds(i)-ac_fds(i));
end

%The following figures show what it is suggested in Zajac's paper that the
%envelope of the rectified EMG signal can be considered as the u(t)
%function and the low pass filtered rectified EMG can be related to the net
%activation of the muscle a(t).
figure
plot(tim,ac_edc,tEMG,yedc1,tEMG,a_edc_ode);
xlabel("Time (s)");
ylabel("Muscle activation");
legend("activation formula","activation from filtering","activation ode");
title("Muscle activation of EDC Muscle");


figure
plot(tim,ac_fdp,tEMG,yfdp1,tEMG,a_fdp_ode);
xlabel("Time (s)");
ylabel("Muscle activation");
legend("activation formula","activation from filtering","activation ode");
title("Muscle activation of FDP Muscle");

figure
plot(tim,ac_fds,tEMG,yfds1,tEMG,a_fds_ode);
xlabel("Time (s)");
ylabel("Muscle activation");
legend("activation formula","activation from filtering","activation ode");
title("Muscle activation of FDS Muscle");

    
    clear ac_edc ac_fdp ac_fds
%Use the activations from the ODE solver

ac_edc=downsample(a_edc_ode,ceil(fsEMG/fs));
ac_fdp=downsample(a_fdp_ode,ceil(fsEMG/fs));
ac_fds=downsample(a_fds_ode,ceil(fsEMG/fs));


figure
plot(tim,ac_edc,tim,downsample(yedc1,ceil(fsEMG/fs)));
xlabel("Time (s)");
ylabel("Muscle activation");
legend("activation formula","activation from filtering");
title("Muscle activation of EDC Muscle");


figure
plot(tim,ac_fdp,tim,downsample(yfdp1,ceil(fsEMG/fs)));
xlabel("Time (s)");
ylabel("Muscle activation");
legend("activation formula","activation from filtering");
title("Muscle activation of FDP Muscle");

figure
plot(tim,ac_fds,tim,downsample(yfds1,ceil(fsEMG/fs)));
xlabel("Time (s)");
ylabel("Muscle activation");
legend("activation formula","activation from filtering");
title("Muscle activation of FDS Muscle");




%From Zajac the muscle force production is a*Fact+Fpass=Fm and the muscle
%moment is M=Fm*r

M_EDC_MCP=rEDC_MCP.*(ac_edc.*FA_EDC_MCP+PA_EDC_MCP);
M_FDP_MCP=rFDP_MCP.*(ac_fdp.*FA_FDP_MCP+PA_FDP_MCP);
M_FDS_MCP=rFDS_MCP.*(ac_fds.*FA_FDS_MCP+PA_FDS_MCP);

M_EDC_PIP=rEDC_PIP.*(ac_edc.*FA_EDC_PIP+PA_EDC_PIP);
M_FDP_PIP=rFDP_PIP.*(ac_fdp.*FA_FDP_PIP+PA_FDP_PIP);
M_FDS_PIP=rFDS_PIP.*(ac_fds.*FA_FDS_PIP+PA_FDS_PIP);

M_EDC_DIP=rEDC_DIP.*(ac_edc.*FA_EDC_DIP+PA_EDC_DIP);
M_FDP_DIP=rFDP_DIP.*(ac_fdp.*FA_FDP_DIP+PA_FDP_DIP);
Res=[M_EDC_MCP M_FDP_MCP M_FDS_MCP M_EDC_PIP M_FDP_PIP M_FDS_PIP M_EDC_DIP M_FDP_DIP];

%The sum of all the moments. The extensor moments are already negative!
SF_MCP=M_EDC_MCP+M_FDP_MCP+M_FDS_MCP;
SF_PIP=M_EDC_PIP+M_FDP_PIP+M_FDS_PIP;
SF_DIP=M_EDC_DIP+M_FDP_DIP;

clear FA_EDC_MCP FA_EDC_PIP FA_EDC_DIP FA_FDS_MCP FA_FDS_PIP FA_FDP_MCP FA_FDP_PIP FA_FDP_DIP 

clear PA_EDC_MCP PA_EDC_PIP PA_EDC_DIP PA_FDS_MCP PA_FDS_PIP PA_FDP_MCP PA_FDP_PIP PA_FDP_DIP

clear rEDC_MCP rEDC_PIP rEDC_DIP rFDS_MCP rFDS_PIP rFDP_MCP rFDP_PIP rFDP_DIP
    
end

