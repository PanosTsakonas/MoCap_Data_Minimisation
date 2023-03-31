function [tEMG,Res,aEPLf,aFPLf]=Thumb_EMG_2_Fit(th1f,th2f,th3f,aEPL,aFPL,aEPLMax,aFPLMax,fs)
disp("I am at the Thumb_EMG_2_Fit function");
syms x

%EMG sampling frequency
fsEMG=1000;
fnEMG=fsEMG/2;

%Setting the time vector for the EMG signals
tEMG=zeros(length(aEPL),1);
for i=2:length(aEPL)
    tEMG(i)=tEMG(i-1)+1/fsEMG;
end

%Thumb finger parameters taken from "A Musculoskeletal Model of the Hand and Wrist Capable of Simulating
%Functional Tasks" and the data came from the Sign Language model.

 %CMC
 
 %Moment arm
 rEPL_CMC=feval(matlabFunction(0.0155*x^2-0.004558*x-0.008711),th1f);
 rFPL_CMC=feval(matlabFunction(-0.09666*x^4-0.1898*x^3-0.1383*x^2-0.04315*x+0.01526),th1f);
 
 %Active force
 FA_EPL_CMC=feval(matlabFunction(-5.092*x^3+0.09232*x^2+5.79*x+86.8),th1f);
 FA_FPL_CMC=feval(matlabFunction(-151.8*x^3-167*x^2+117.1*x+172.9),th1f);
 
 %Passive force
 PA_EPL_CMC=zeros(length(th1f),1);
 PA_FPL_CMC=feval(matlabFunction(-30.45*x^3+47.18*x^2-11.3*x+2.358),th1f);
 
 %MCP
 
 %Moment arm
 rEPL_MCP=feval(matlabFunction(0.002074*x^3+0.000608*x^2-0.004506*x-0.009047),th2f);
 rFPL_MCP=feval(matlabFunction(-0.00593*x^3+0.002096*x^2+0.006769*x+0.01065),th2f);
 %Active force
 FA_EPL_MCP=feval(matlabFunction(0.3081*x^4-4.36*x^3-2.788*x^2+5.565*x+86.87),th2f);
 FA_FPL_MCP=feval(matlabFunction(9.44*x^4-55.54*x^3-32.96*x^2+83.96*x+175.2),th2f);
 
 %Passive Force
 PA_EPL_MCP=feval(matlabFunction(0.8953*x^4+0.4616*x^3-0.189*x^2-0.0714*x+0.0042),th2f);
 PA_FPL_MCP=feval(matlabFunction(-2.499*x^3+7.217*x^2-7.578*x+2.574),th2f);
 
 %IP
 
 %Moment arm
 rEPL_IP=feval(matlabFunction(-0.00102*x^3+0.003302*x^2-0.001067*x-0.00515),th3f);
 rFPL_IP=feval(matlabFunction(-0.004318*x^5+0.005549*x^4+0.002303*x^3-0.003717*x^2-0.0008377*x+0.008827),th3f);
 
 %Active force
 FA_EPL_IP=feval(matlabFunction(1.854*x^4-3.256*x^3-0.8478*x^2+3.7*x+86.83),th3f);
 FA_FPL_IP=feval(matlabFunction(32.2*x^4-48.89*x^3-32.32*x^2+74.82*x+174.6),th3f);
 
 %Passive Force
 PA_EPL_IP=feval(matlabFunction(0.1575*x^5-0.2144*x^4+0.09181*x^3+0.04446*x^2-0.004481*x-0.00125),th3f);
 PA_FPL_IP=feval(matlabFunction(-1.662*x^3+5.531*x^2-6.501*x+2.679),th3f);
 
 
%The next thing to be imported will be the EMG data. According to "Feasibility of using 
%combined EMG and kinematic signals for prosthesiscontrol: A simulation study 
%using a virtual reality environment" the EMG signals were ampli-
%fied, band-pass filtered between 15 and 450 Hz. This is also supported by
%plotting the Hilber Huang transform of the EMG signals


%The band pass filter

[bemg,aemg]=butter(4,[15 450]./fnEMG,'bandpass');
aEPLf=filtfilt(bemg,aemg,aEPL);
aFPLf=filtfilt(bemg,aemg,aFPL);
clear bemg aemg

%Removing the 50Hz power line frequency
[bemg,aemg]=butter(4,[48 52]./fnEMG,'stop');
aEPLf=filtfilt(bemg,aemg,aEPLf);
aFPLf=filtfilt(bemg,aemg,aFPLf);
clear bemg aemg

%Rectify by taking the absolute of the EMG data as suggested
%in "Finger joint coordination during tapping"
aEPLf=abs(aEPLf);
aFPLf=abs(aFPLf);

%Normalize
if isempty(aEPLMax)==1
aEPLf=aEPLf./max(aEPLf);
else
    aEPLf=aEPLf./aEPLMax;
end

if isempty(aFPLMax)==1
aFPLf=aFPLf./max(aFPLf);
else
    aFPLf=aFPLf./aFPLMax;
end


%Obtain the envelope of the normalised original EMG signal
yepl=envelope(aEPLf,15,'peak');
yfpl=envelope(aFPLf,15,'peak');

%Use the activation dynamics differential equation to obtain the
%activations from integration using ode45 and spline interpolation of the
%normalised muscle excitations

pepl=spline(tEMG,yepl);
pfpl=spline(tEMG,yfpl);

ppepl=mkpp(pepl.breaks,pepl.coefs);
ppfpl=mkpp(pfpl.breaks,pfpl.coefs);

%Solve the differential equation shown in 
[~,a_epl_ode]=ode45(@(t,a) activation(t,a,ppepl),tEMG,0);
[~,a_fpl_ode]=ode45(@(t,a) activation(t,a,ppfpl),tEMG,0);


ac_epl=downsample(a_epl_ode,ceil(fsEMG/fs));
ac_fpl=downsample(a_fpl_ode,ceil(fsEMG/fs));

figure
plot(tEMG,a_epl_ode);
legend("EPL muscle activation");
xlabel("Time (s)");
ylabel("Muscle activation");
title("Muscle activation for the EPL muscle");

figure
plot(tEMG,a_fpl_ode);
legend("FPL muscle activation");
xlabel("Time (s)");
ylabel("Muscle activation");
title("Muscle activation for the FPL muscle");

%From Zajac the muscle force production is a*Fact+Fpass=Fm and the muscle
%moment is M=Fm*r

M_EPL_CMC=rEPL_CMC.*(ac_epl.*FA_EPL_CMC+PA_EPL_CMC);
M_FPL_CMC=rFPL_CMC.*(ac_fpl.*FA_FPL_CMC+PA_FPL_CMC);


M_EPL_MCP=rEPL_MCP.*(ac_epl.*FA_EPL_MCP+PA_EPL_MCP);
M_FPL_MCP=rFPL_MCP.*(ac_fpl.*FA_FPL_MCP+PA_FPL_MCP);


M_EPL_IP=rEPL_IP.*(ac_epl.*FA_EPL_IP+PA_EPL_IP);
M_FPL_IP=rFPL_IP.*(ac_fpl.*FA_FPL_IP+PA_FPL_IP);

Res=[M_EPL_CMC M_FPL_CMC M_EPL_MCP M_FPL_MCP M_EPL_IP M_FPL_IP];