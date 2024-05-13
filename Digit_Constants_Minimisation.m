%This code will be used for determining the spring and damper constants
%during dynamic movement of the fingers. A simple flexion/extension movement will be used
%and the muscles that will be considered are the FDP,FDS, EDC, EPL and FPL. Data for the Muscle 
%moments are taken from the "A Musculoskeletal Model of the Hand and Wrist
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

%If you want the data to use the linear envelope change the following to 1
env=[];

%mean Voltage of MVC for Normalisation could be either an empty matrix or a value
aEDCMax=0.30772;
aFDPMax=0.2008;
aFDSMax=0.15926;

trial=input("Which Cylindrical trial number you are doing: ");
if trial==0
    trial=" ";
end
prom=input("is it a repeat trial? y/n: ",'s');

if (prom=='Y' || prom=='y')
    cyl="repeat";
else
    cyl=" ";
end

%gravitational acceleration m/s^2
g=9.81;

%Set the sampling frequency of the MoCap
fs=125;
fn=fs/2;

%Human body density in Kg/m^3
rho=1.16*10^3;

[coll,rows]=size(th1);

th1f=zeros(coll,rows);
th2f=th1f;
th3f=th2f;

wn_mcp=zeros(1,rows);
wn_dip=wn_mcp;
wn_pip=wn_mcp;
n_mcp=wn_mcp;
n_pip=n_mcp;
n_dip=n_pip;

%Set the time matrix
tim=zeros(length(th1),1);
for i=2:length(th1)
    tim(i)=tim(i-1)+1/fs;
end


for i=1:4
%Low pass filter the angular data.
wn_mcp(i)=Residual_analysis_for_filtering(th1(:,i),fs,'MCP',i);
[b,a]=butter(4,wn_mcp(i)/fn,'low');
th1f(:,i)=filtfilt(b,a,th1(:,i)).*pi/180;
clear b a
wn_pip(i)=Residual_analysis_for_filtering(th2(:,i),fs,'PIP',i);
[b,a]=butter(4,wn_pip(i)/fn,'low');
th2f(:,i)=filtfilt(b,a,th2(:,i)).*pi/180;
clear b a
wn_dip(i)=Residual_analysis_for_filtering(th3(:,i),fs,'DIP',i);
[b,a]=butter(4,wn_dip(i)/fn,'low');
th3f(:,i)=filtfilt(b,a,th3(:,i)).*pi/180;
clear b a



disp("Important reminder: When MCP angular data have a sharp edge close to zero this is interpreted as finger motion starting from an extended joint position. Convert the angles before zero to negative. OpenSim muscle data take into account this negative sign and corresponds to an extension of the joint");

figure
subplot(3,1,1)
hold on;
plot(tim,th1f(:,i),tim,th1(:,i).*pi/180,'x');
title("MCP");
xlabel("Time (s)");
ylabel("Angle (rad)");
sh=input("The data exhibit a sharp edge close to zero? Y/N: ",'s');
if (sh=='Y' || sh=='y')
    [th1m,indm]=min(th1f(:,i));
    
    th1f(:,i)=[-th1f(1:indm-1,i); th1f(indm:end,i)];
    plot(tim,th1f(:,i),'o');
    legend("Filtered at "+wn_mcp(i)+" Hz","Raw data","Filtered data with OpenSim convention",'location','southeast');
else
legend("Filtered at "+wn_mcp(i)+" Hz","Raw data",'location','southeast');
end
hold off;
subplot(3,1,2)
plot(tim,th2f(:,i),tim,th2(:,i).*pi/180,'x');
xlabel("Time (s)");
ylabel("Angle (rad)");
legend("Filtered at "+wn_pip(i)+" Hz","Raw data",'location','southeast');
title("PIP");
subplot(3,1,3)
plot(tim,th3f(:,i),tim,th3(:,i).*pi/180,'x');
xlabel("Time (s)");
ylabel("Angle (rad)");
legend("Filtered at "+wn_dip(i)+" Hz","Raw data",'location','southeast');
title("DIP");

syms t x

%To avoid fitting when no movement occurs in the begining set the index
%point for each trial when the digit actually moves

figure
plot(th1f(:,i))
title("MCP filtered data");
xlabel("Index");
ylabel("Angle (rad)");
n_mcp(i)=input("Give the index where the motion begins for the MCP joint: ");
figure
plot(th2f(:,i))
title("PIP filtered data");
xlabel("Index");
ylabel("Angle (rad)");
n_pip(i)=input("Give the index where the motion begins for the PIP joint: ");
figure
plot(th3f(:,i))
title("DIP filtered data");
xlabel("Index");
ylabel("Angle (rad)");
n_dip(i)=input("Give the index where the motion begins for the DIP joint: ");

end

D=["Thumb","Index","Middle","Ring","Little"];
Par=input("Give the participant number: ");

for in=2:5

I=importfile("C:\Users\............Moment of Inertia excel file for each participant.xlsx","");

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

[SF_MCP,SF_PIP,SF_DIP,tEMG,Res,edcf,fdpf,fdsf]=EMG_2_Fit(th1f(:,in-1),th2f(:,in-1),th3f(:,in-1),aEDC,aFDP,aFDS,aEDCMax,aFDPMax,aFDSMax,tim,fs,in,env);



%Obtain the spline interpolation to each muscle moment for each degree of
%freedom
a=0;
while(a~=1)

    [EPL_CMC, FPL_CMC,~,EPL_MCP,FPL_MCP,~,EPL_IP,FPL_IP]=Moments(Res,tim,N,IBK,Lag,in);
    a=1;
    pause(20);
end


%close all;

%Initial conditions for integration. Movement starts from rest.
init=[th1f(n_mcp(in-1),in-1),0,th2f(n_pip(in-1),in-1),0,th3f(n_dip(in-1),in-1),0];

%Set the minimisation options
options=optimoptions('fmincon','MaxFunctionEvaluations',10^3,'MaxIterations',10^4,'ConstraintTolerance',10^-4,'FunValCheck','on','StepTolerance',1e-4,'UseParallel',true);

if isempty(IBK)==0
if in~=1
%Set the initial conditions Damper, Spring, and Moment scaling parameters
r1=[0.0013 0.0482 1 1 1];
r2=[5.26E-04 0.0169 1 1 1];
r3=[8.23E-05 0.002814 1 1];

theq=[th1f(end,in-1) th2f(end,in-1) th3f(end,in-1)];
i=1;
% Define the number of bootstrap samples
B1=5;

% Define the number of parameters
P=length(r1);
P2=length(r3);

% Initialize matrices to store parameter estimates and bootstrap datasets
par_mcp=zeros(B1,P);
par_pip=zeros(B1,P);
par_dip=zeros(B1,P2);

Min1=[];
Min2=[];
Min3=[];

Boot=0;

if isempty(Min1)==1
[Min1,sum1,~,out_mcp,lamda_mcp,~,hess_mcp]=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EDC_MCP,FDP_MCP,FDS_MCP,g,r,tim(n_mcp(in-1):end),init,th1f(n_mcp(in-1):end,in-1),i,I,theq(1),grav,in,Boot),r1,[],[],[],[],[r1(1) r1(2) 0.1 0.1 0.1],[2 10 2 2 2],[],options);
end

%Solve the IBK model for the minimised parameters
[~,Y1]=ode45(@(t,y) IBK_th1(M1,M2,M3,L1,EDC_MCP,FDP_MCP,FDS_MCP,g,Min1,I1,theq(1),grav,t,y),tim(n_mcp(in-1):end),[init(1) init(2)]);
low=[r1(1) r1(2) 0.1 0.1 0.1];
up=[2 10 2 2 2];
Boot=1;
residuals=th1f(n_mcp(in-1):end,in-1)-Y1(:,1);  % Calculate residuals from initial parameter estimates
for i1=1:B1
    disp("Bootstrap Method for MCP joint "+i1+"/"+B1);
        bootstrap_residuals=datasample(residuals, length(residuals), 'Replace', true); 
        bootstrap_y=Y1(:,1) + bootstrap_residuals;
        par_mcp(i1,:)=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EDC_MCP,FDP_MCP,FDS_MCP,g,r,tim(n_mcp(in-1):end),init,bootstrap_y,i,I,theq(1),grav,in,Boot),Min1,[],[],[],[],low,up,[],options);
        
end
% Calculate the parameter uncertainties as the standard deviation of parameter estimates
sd_mcp = std(par_mcp);

% Calculate the 95% bootstrap confidence intervals
bootstrap_CIs = prctile(par_mcp, [2.5, 97.5]);

%Middle segment IBK optimisation
i=2;
Boot=0;
if isempty(Min2)==1
[Min2,sum2,~,out_pip,lamda_pip,~,hess_pip]=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EDC_PIP,FDP_PIP,FDS_PIP,g,r,tim(n_pip(in-1):end),init,th2f(n_pip(in-1):end,in-1),i,I,theq(2),grav,in,Boot),r2,[],[],[],[],[r2(1) r2(2) 0.1 0.1 0.1],[2 10 2 2 2],[],options);
end

clear low up residuals
%Solve the IBK model for the minimised parameters
[~,Y2]=ode45(@(t,y) IBK_th2(M2,M3,L2,EDC_PIP,FDP_PIP,FDS_PIP,g,Min2,I2,theq(2),grav,t,y),tim(n_pip(in-1):end),[init(3) init(4)]);
low=[r2(1) r2(2) 0.1 0.1 0.1];
up=[2 10 2 2 2];
Boot=1;
residuals=th2f(n_pip(in-1):end,in-1)-Y2(:,1);% Calculate residuals from initial parameter estimates
for i1=1:B1
    disp("Bootstrap Method for PIP joint "+i1+"/"+B1);
        bootstrap_residuals=datasample(residuals, length(residuals), 'Replace', true);
        bootstrap_y=Y2(:,1) + bootstrap_residuals;
        par_pip(i1,:)=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EDC_PIP,FDP_PIP,FDS_PIP,g,r,tim(n_pip(in-1):end),init,bootstrap_y,i,I,theq(2),grav,in,Boot),Min2,[],[],[],[],low,up,[],options);
        
end
% Calculate the parameter uncertainties as the standard deviation of parameter estimates
sd_pip = std(par_pip);

% Calculate the 95% bootstrap confidence intervals
bootstrap_CIs_pip = prctile(par_pip, [2.5, 97.5]);

%Distal segment IBK optimisation
i=3;
Boot=0;

if isempty(Min3)==1
[Min3,sum3,~,out_dip,lamda_dip,~,hess_dip]=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EDC_DIP,FDP_DIP,[],g,r,tim(n_dip(in-1):end),init,th3f(n_dip(in-1):end,in-1),i,I,theq(3),grav,in,Boot),r3,[],[],[],[],[r3(1) r3(2) 0.1 0.1],[2 10 2 2],[],options);
end

clear low up residuals
[~,Y3]=ode45(@(t,y) IBK_th3(M3,L3,EDC_DIP,FDP_DIP,g,Min3,I3,theq(3),grav,t,y),tim(n_dip(in-1):end),[init(5) init(6)]);
low=[r3(1) r3(2) 0.1 0.1];
up=[2 10 2 2];
Boot=1;
residuals=th3f(n_dip(in-1):end,in-1)-Y3(:,1);  % Calculate residuals from initial parameter estimates
for i1=1:B1
    disp("Bootstrap Method for DIP joint "+i1+"/"+B1);
        bootstrap_residuals=datasample(residuals, length(residuals), 'Replace', true);
        bootstrap_y=Y3(:,1) + bootstrap_residuals;
        par_dip(i1,:)=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EDC_DIP,FDP_DIP,[],g,r,tim(n_dip(in-1):end),init,bootstrap_y,i,I,theq(3),grav,in,Boot),Min3,[],[],[],[],low,up,[],options);
        
end
% Calculate the parameter uncertainties as the standard deviation of parameter estimates
sd_dip = std(par_dip);

% Calculate the 95% bootstrap confidence intervals
bootstrap_CIs_dip = prctile(par_dip, [2.5, 97.5]);


%Plot results
figure
plot(tim(n_mcp(in-1):end),Y1(:,1).*180/pi,tim(n_mcp(in-1):end),th1f(n_mcp(in-1):end,in-1).*180/pi,'x');
rmse1=sqrt(mean((Y1(:,1).*180/pi-th1f(n_mcp(in-1):end,in-1).*180/pi).^2));
rsq1=1-sum((th1f(n_mcp(in-1):end,in-1)-Y1(:,1)).^2)/sum((th1f(n_mcp(in-1):end,in-1)-mean(th1f(n_mcp(in-1):end,in-1))).^2);
legend("Minimisation","Gait Lab MCP data",'Location','southeast');
xlabel("Time (s)");
ylabel("Angle (degrees)");
title("Angular data and fit "+D(in)+" with R^2: "+rsq1+" RMSE: "+rmse1+" degrees");


figure
plot(tim(n_pip(in-1):end),Y2(:,1).*180/pi,tim(n_pip(in-1):end),th2f(n_pip(in-1):end,in-1).*180/pi,'x');
rmse2=sqrt(mean((Y2(:,1).*180/pi-th2f(n_pip(in-1):end,in-1).*180/pi).^2));
rsq2=1-sum((th2f(n_pip(in-1):end,in-1)-Y2(:,1)).^2)/sum((th2f(n_pip(in-1):end,in-1)-mean(th2f(n_pip(in-1):end,in-1))).^2);
legend("Minimisation","Gait Lab PIP data",'Location','southeast');
xlabel("Time (s)");
ylabel("Angle (degrees)");
title("Angular data and fit "+D(in)+" with R^2: "+rsq2+" RMSE: "+rmse2+" degrees");


figure
plot(tim(n_dip(in-1):end),Y3(:,1).*180/pi,tim(n_dip(in-1):end),th3f(n_dip(in-1):end,in-1).*180/pi,'x');
rmse3=sqrt(mean((Y3(:,1).*180/pi-th3f(n_dip(in-1):end,in-1).*180/pi).^2));
rsq3=1-sum((th3f(n_dip(in-1):end,in-1)-Y3(:,1)).^2)/sum((th3f(n_dip(in-1):end,in-1)-mean(th3f(n_dip(in-1):end,in-1))).^2);
legend("Minimisation","Gait Lab DIP data",'Location','southeast');
xlabel("Time (s)");
ylabel("Angle (degrees)");
title("Angular data and fit "+D(in)+" with R^2: "+rsq3+" RMSE: "+rmse3+" degrees");

rsq=[rsq1 rsq2 rsq3];
rmse=[rmse1 rmse2 rmse3];
save_resultsM(Min1,Min2,Min3,mn1,mn2,mn3,in,Par,sd_mcp,sd_pip,sd_dip,trial,rsq,rmse,cyl,wn_mcp(in-1),wn_pip(in-1),wn_dip(in-1),env);
if (length(Y1(:,1))==length(th1) && length(Y2(:,1))==length(th2) && length(Y3(:,1))==length(th3))
    if isempty(env)==1
        if LL=='l'
        ffsave="C:\Users\panos\OneDrive - University of Warwick\PhD\Hand Trials\Results\Cylindrical Grasp\P"+Par+"\data_digit_"+in+"par_"+Par+"_Cylindical_"+cyl+"_trial_"+trial+".xlsx";
        else
            ffsave="C:\Users\u1857308\OneDrive - University of Warwick\PhD\Hand Trials\Results\Cylindrical Grasp\P"+Par+"\data_digit_"+in+"par_"+Par+"_Cylindical_"+cyl+"_trial_"+trial+".xlsx";
        end
    else
            if LL=='l'
        
             ffsave="C:\Users\panos\OneDrive - University of Warwick\PhD\Hand Trials\Results\Cylindrical Grasp\P"+Par+"\data_digit_"+in+"par_"+Par+"_Cylindical_"+cyl+"_trial_"+trial+"Linear_Envelope.xlsx";
            else
                ffsave="C:\Users\u1857308\OneDrive - University of Warwick\PhD\Hand Trials\Results\Cylindrical Grasp\P"+Par+"\data_digit_"+in+"par_"+Par+"_Cylindical_"+cyl+"_trial_"+trial+"Linear_Envelope.xlsx"; 
            end
    end
        writematrix([ppval(EDC_MCP,tim(n_mcp(in-1):end)),ppval(FDP_MCP,tim(n_mcp(in-1):end)),ppval(FDS_MCP,tim(n_mcp(in-1):end)),ppval(EDC_PIP,tim(n_pip(in-1):end)),ppval(FDP_PIP,tim(n_pip(in-1):end)),ppval(FDS_PIP,tim(n_pip(in-1):end)),ppval(EDC_DIP,tim(n_dip(in-1):end)),ppval(FDP_DIP,tim(n_dip(in-1):end)),tim(n_mcp(in-1):end),th1f(n_mcp(in-1):end,in-1),tim(n_pip(in-1):end),th2f(n_pip(in-1):end,in-1),tim(n_dip(in-1):end),th3f(n_dip(in-1):end,in-1),Y1(:,1),Y2(:,1),Y3(:,1)],ffsave,'Sheet',in);
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
pause(15);
end
