%This code will be used for determining the spring and damper constants
%during dynamic movement of the fingers. A simple flexion/extension movement will be used
%and the muscles that will be considered are the FDP,FDS, EDC, EPL and FPL. Data for the Muscle moments are taken from
%the "A Musculoskeletal Model of the Hand and Wrist
%Capable of Simulating Functional Tasks". The idea is to minimze the
%difference between the observed angles and those determined from the
%itterative solution of the extended IBK model/Lagrangian. User defined inputs are the
%selection of minimising the IBK approximation or the complete Lagrangian
%function. If the experiment is performed in the direction of gravity this
%can also be specified. The parameters of the Maximum Voluntary Contraction
%(MVC) for the EDC, FDP, FDS,EPL and FPL muscles can  also be given if the
%normalisation of the signal is done under the MVC criterion. Otherwise the
%criterion of the maximum peak of the signal is used to normalise the EMG
%recordings.
global mn1 mn2 mn3;
close all;
clc;
%Select if you want to minimise the IBK or the Lagrangian. Set either to 1
%if you want to minimise the respective equations
IBK=1;
Lag=[];

%If the motion is in the direction of gravity change the following to 1
grav=[];

%mean Voltage of MVC for Normalisation could be either an empty matrix or a value
aEDCMax=0.38847;
aFDPMax=0.28946;
aFDSMax=0.18987;
aEPLMax=[];
aFPLMax=[];

%gravitational acceleration m/s^2
g=9.81;

%Set the sampling frequency of the MoCap
fs=125;
fn=fs/2;

%Human body density in Kg*/m^3
rho=1.16*10^3;


%Low pass filter the angular data.
[b,a]=butter(4,10/fn,'low');
th1f=filtfilt(b,a,th1).*pi/180;
th2f=filtfilt(b,a,th2).*pi/180;
th3f=filtfilt(b,a,th3).*pi/180;


%Set the time matrix
tim=zeros(length(th1),1);
for i=2:length(th1)
    tim(i)=tim(i-1)+1/fs;
end

figure
subplot(3,1,1)
plot(tim,th1f,tim,th1.*pi/180,'x');
title("MCP");
xlabel("Time (s)");
ylabel("Angle (rad)");
legend("Filtered","Raw data",'location','southeast');
subplot(3,1,2)
plot(tim,th2f,tim,th2.*pi/180,'x');
xlabel("Time (s)");
ylabel("Angle (rad)");
legend("Filtered","Raw data",'location','southeast');
title("PIP");
subplot(3,1,3)
plot(tim,th3f,tim,th3.*pi/180,'x');
xlabel("Time (s)");
ylabel("Angle (rad)");
legend("Filtered","Raw data",'location','southeast');
title("DIP");
syms t x

Par=input("Give the participant number: ");
in=input("Give the digit you are working with: ");

%Uncomment the following when I am working from my Laptop
I=importfile("Set file directory for participant segment length and moment of inertia data","");

%Import the parameters from the excel files.
I1=I.data(27+in,4);
I2=I.data(27+in,5);
I3=I.data(27+in,6);
Ia=I.data(27+in,7);
L=[I.data(in,6),I.data(in,4),I.data(in,2)].*10^-3;
m=[I.data(in,10),I.data(in,9),I.data(in,8)];
M1=m(1);
M2=m(2);
M3=m(3);
L1=L(1);
L2=L(2);
L3=L(3);
I=[I1 I2 I3];

%Function that calculates the muscle moments based on the
%M=Sum(r_i*(a_i*Fact,i+Fpass,i)) formula for each muscle
if in==1
   [tEMG,Res,eplf,fplf]=Thumb_EMG_2_Fit(th1f,th2f,th3f,aEPL,aFPL,aEPLMax,aFPLMax,fs); 
else

[SF_MCP,SF_PIP,SF_DIP,tEMG,Res,edcf,fdpf,fdsf]=EMG_2_Fit(th1f,th2f,th3f,aEDC,aFDP,aFDS,aEDCMax,aFDPMax,aFDSMax,tim,fs,in);
end


%Obtain the spline interpolation to each muscle moment for each degree of
%freedom
a=0;
while(a~=1)
    if isempty(IBK)==0
        N=[];
    end
    if isempty(Lag)==0
        N=input("Select a polynomial degree order: ");
    end

if in~=1
[EDC_MCP, FDP_MCP, FDS_MCP, EDC_PIP, FDP_PIP, FDS_PIP, EDC_DIP, FDP_DIP]=Moments(Res,tim,N,IBK,Lag,in);
else
    [EPL_CMC, FPL_CMC,~,EPL_MCP,FPL_MCP,~,EPL_IP,FPL_IP]=Moments(Res,tim,N,IBK,Lag,in);
end

a=input("Press 1 to continue: ");
end


close all;

%Initial conditions for integration. Movement starts from rest.
init=[th1f(1),0,th2f(1),0,th3f(1),0];

%Set the minimisation options
options=optimoptions('fmincon','MaxFunctionEvaluations',10^10,'MaxIterations',10^10,'ConstraintTolerance',10^-10,'FunValCheck','on','StepTolerance',1e-10);

if isempty(IBK)==0
if in~=1
%Set the initial conditions Damper, Spring, and Moment scaling parameters
r1=[0.0013 0.0482 1 1 1];
r2=[5.26E-04 0.0169 1 1 1];
r3=[8.23E-05 0.002814 1 1];

TH=importfile("Set file directory for equilibrium angles with no gravity","");
%Set the equilibrium angles in degrees determined from a static capture
if in==2
 theq=[mean(TH.data(5:end,59),'omitnan')  mean(TH.data(5:end,35),'omitnan') mean(TH.data(5:end,14),'omitnan')].*pi/180;
elseif in==3
 theq=[mean(TH.data(5:end,62),'omitnan')  mean(TH.data(5:end,41),'omitnan') mean(TH.data(5:end,17),'omitnan')].*pi/180;
 elseif in==4
 theq=[mean(TH.data(5:end,65),'omitnan')  mean(TH.data(5:end,47),'omitnan') mean(TH.data(5:end,20),'omitnan')].*pi/180;
else
 theq=[mean(TH.data(5:end,68),'omitnan')  mean(TH.data(5:end,53),'omitnan') mean(TH.data(5:end,23),'omitnan')].*pi/180;
end

%Proximal segment IBK optimisation
i=1;
[Min1,sum1]=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EDC_MCP,FDP_MCP,FDS_MCP,g,r,tim,init,th1f,i,I,theq(1),grav,in),r1,[],[],[],[],[r1(1) r1(2) 0.1 0.1 0.1],[2 10 2 2 2],[],options);

%Middle segment IBK optimisation
i=2;
[Min2,sum2]=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EDC_PIP,FDP_PIP,FDS_PIP,g,r,tim,init,th2f,i,I,theq(2),grav,in),r2,[],[],[],[],[r2(1) r2(2) 0.1 0.1 0.1],[2 10 1.8 1.8 1.8],[],options);

%Distal segment IBK optimisation
i=3;
[Min3,sum3]=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EDC_DIP,FDP_DIP,[],g,r,tim,init,th3f,i,I,theq(3),grav,in),r3,[],[],[],[],[r3(1) r3(2) 0.1 0.1],[2 5 1.8 1.8],[],options);


[~,Y1]=ode45(@(t,y) IBK_th1(M1,M2,M3,L1,EDC_MCP,FDP_MCP,FDS_MCP,g,Min1,I1,theq(1),grav,t,y),tim,[init(1) init(2)]);
figure
plot(tim,Y1(:,1).*180/pi,tim,th1f.*180/pi);
rmse=sqrt(mean((Y1(:,1).*180/pi-th1f.*180/pi).^2));
rsq=1-sum((th1f-Y1(:,1)).^2)/sum((th1f-mean(th1f)).^2);
legend("Minimisation","Gait Lab MCP data",'Location','southeast');
xlabel("Time (s)");
ylabel("Angle (degrees)");
title("Angular data and fit with R^2: "+rsq+" RMSE: "+rmse+" degrees");

[~,Y2]=ode45(@(t,y) IBK_th2(M2,M3,L2,EDC_PIP,FDP_PIP,FDS_PIP,g,Min2,I2,theq(2),grav,t,y),tim,[init(3) init(4)]);
figure
plot(tim,Y2(:,1).*180/pi,tim,th2f.*180/pi);
rmse=sqrt(mean((Y2(:,1).*180/pi-th2f.*180/pi).^2));
rsq=1-sum((th2f-Y2(:,1)).^2)/sum((th2f-mean(th2f)).^2);
legend("Minimisation","Gait Lab PIP data",'Location','southeast');
xlabel("Time (s)");
ylabel("Angle (degrees)");
title("Angular data and fit with R^2: "+rsq+" RMSE: "+rmse+" degrees");

[~,Y3]=ode45(@(t,y) IBK_th3(M3,L3,EDC_DIP,FDP_DIP,g,Min3,I3,theq(3),grav,t,y),tim,[init(5) init(6)]);
figure
plot(tim,Y3(:,1).*180/pi,tim,th3f.*180/pi);
rmse=sqrt(mean((Y3(:,1).*180/pi-th3f.*180/pi).^2));
rsq=1-sum((th3f-Y3(:,1)).^2)/sum((th3f-mean(th3f)).^2);
legend("Minimisation","Gait Lab DIP data",'Location','southeast');
xlabel("Time (s)");
ylabel("Angle (degrees)");
title("Angular data and fit with R^2: "+rsq+" RMSE: "+rmse+" degrees");

Write=input("Do you want the data for checking the Lagrangian equations?:Y/N ");

if(Write==Y || Write==y)
   writematrix([ppval(EDC_MCP,tim),ppval(FDP_MCP,tim),ppval(FDS_MCP,tim),ppval(EDC_PIP,tim),ppval(FDP_PIP,tim),ppval(FDS_PIP,tim),ppval(EDC_DIP,tim),ppval(FDP_DIP,tim),tim,th1f,th2f,th3f,Y1(:,1),Y2(:,1),Y3(:,1)],"data_digit_"+in+"_par_"+Par+".csv");
end


else
    %Set the initial parameter values
    rt1=[0.03 0.09 1 1];
    rt2=[0.001 0.05 1 1];
    rt3=[0.0005 0.03 1 1];
    
    %Set the equilibrium angles determined from a static capture
 theq=[40.93 43.77 26.75].*pi/180;
 
    i=1;
    [Mincmc,sumcmc]=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EPL_CMC,FPL_CMC,[],g,r,tim,init,th1f,i,I,theq(1),grav,in),rt1,[],[],[],[],[rt1(1) rt1(2) 0.1 0.1],[2 10 1 1],[],options);

    i=2;
    
    [Minmcp,summcp]=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EPL_MCP,FPL_MCP,[],g,r,tim,init,th2f,i,I,theq(2),grav,in),rt2,[],[],[],[],[rt2(1) rt2(2) 0.1 0.1],[2 10 1 1],[],options);

    i=3;
    [Minip,sumip]=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EPL_IP,FPL_IP,[],g,r,tim,init,th3f,i,I,theq(3),grav,in),rt2,[],[],[],[],[rt3(1) rt3(2) 0.1 0.1],[2 10 1 1],[],options);

    [~,Y1]=ode45(@(t,y) IBK_Thumb(M1,M2,M3,L1,L2,L3,I1,I2,I3,EPL_CMC,FPL_CMC,1,g,Mincmc,theq(1),grav,t,y),tim,[init(1) init(2)]);
figure
plot(tim,Y1(:,1),tim,th1f);
legend("Minimisation","Gait Lab CMC data",'Location','southeast');
xlabel("Time (s)");
ylabel("Angle (rad)");
title("Angular data and fit with optimised values");

[~,Y2]=ode45(@(t,y) IBK_Thumb(M1,M2,M3,L1,L2,L3,I1,I2,I3,EPL_MCP,FPL_MCP,2,g,Minmcp,theq(2),grav,t,y),tim,[init(3) init(4)]);
figure
plot(tim,Y2(:,1),tim,th2f);
legend("Minimisation","Gait Lab MCP data",'Location','southeast');
xlabel("Time (s)");
ylabel("Angle (rad)");
title("Angular data and fit with optimised values");

[~,Y3]=ode45(@(t,y) IBK_Thumb(M1,M2,M3,L1,L2,L3,I1,I2,I3,EPL_IP,FPL_IP,3,g,Minip,theq(3),grav,t,y),tim,[init(5) init(6)]);
figure
plot(tim,Y3(:,1),tim,th3f);
legend("Minimisation","Gait Lab IP data",'Location','southeast');
xlabel("Time (s)");
ylabel("Angle (rad)");
title("Angular data and fit with optimised values");

end
end
