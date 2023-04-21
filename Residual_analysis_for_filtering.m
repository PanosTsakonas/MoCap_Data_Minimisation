function wnm=Residual_analysis_for_filtering(th,fs)

%This code performs the residual analysis suggested in the Biomechanics and
%motor control of human movement in section 3.4.4.3 from Winter. This
%function returns the cutoff frequency for the low pass filter determined
%using this method.The intercept a on the ordinate (at 0 Hz) is
%nothing more than the rms value of the noise, because Xˆi for a 0-Hz filter is
%nothing more than the mean of the noise over the N samples. When the data
%consist of true signal plus noise, the residual will be seen to rise above 
%the straight (dashed) line as the cutoff frequency is reduced. 
%This rise above the dashed line represents the signal distortion that is 
%taking place as the cutoff is reduced more and more. The final decision is where fc 
%should be chosen. The compromise is always a balance between the signal distortion and 
%the amount of noise allowed through. If we decide that both should be equal, 
%then we simply project a line horizontally from a to intersect the residual line at b. 
%The frequency chosen is fc', and at this frequency the signal distortion is represented by
%bc. 


clc;
wn=1;
fn=fs/2;
Wn=[1:floor(fn)];
ph=zeros(length(th),floor(fn));
z=zeros(1,floor(fn));

i=1;
while wn<=floor(fn)
    [b,a]=butter(4,wn/fn,'low');
    ph(:,i)=filtfilt(b,a,th);
    wn=wn+1;
    i=i+1;
end

for i=1:floor(fn)
    z(i)=sqrt(mean((th-ph(:,i)).^2));
end

figure
hold on;
plot(Wn,z)
ff=input("Select the frequency point where you want to calculate the linear model: ");
xlabel("Cutoff frequency (Hz)");
ylabel("Residuals (degrees)");
title("Residuals versus cutoff frequency plot");

%Linear fit of the straight portion of the residuals
Model=fitlm([ff:floor(fn)],z(ff:end));

plot([0:fn],feval(Model,[0:fn]),'--');
intercept=feval(Model,0);
plot([0:floor(fn)+1],ones(1,length([0:floor(fn)+1])).*intercept,'-.');
legend("Residuals","Linear fit of the noise residual","RMS of noise "+intercept+" degrees");
hold off;

FF=z<=intercept;

pp=find(FF,1);

a1=abs(intercept-z(pp));
a2=abs(intercept-z(pp-1));

mn=min([a1 a2]);

if mn==a1
    wnm=Wn(pp);
else
    wnm=Wn(pp-1);
end
end



