%This code is used to analyse the MVC data
clc;
close all;
fs=1000;
fn=fs/2;
aedcM=aedcM(1:end-1);
afdsM=afdsM(1:end-1);
afdpM=afdpM(1:end-1);
n=1:length(aedcM);

for i=1:3
    if i==1
        hht(emd(aedcM),fs);
        low=input("Give the low frequency for the band pass filter: ");
        high=input("Give the high frequency for the band pass filter: ");
    
    [b,a]=butter(4,[low high]./fn,'bandpass');
    aedcMax=filtfilt(b,a,aedcM);
    [yedc,~]=envelope(aedcMax,50,'peak');
    
    elseif i==2
        hht(emd(afdpM),fs);
        low=input("Give the low frequency for the band pass filter: ");
        high=input("Give the high frequency for the band pass filter: ");
    [b,a]=butter(4,[low high]./fn,'bandpass');
    afdpMax=filtfilt(b,a,afdpM);
    [yfdp,~]=envelope(afdpMax,50,'peak');
    else
        hht(emd(afdsM),fs);
        low=input("Give the low frequency for the band pass filter: ");
        high=input("Give the high frequency for the band pass filter: ");
    [b,a]=butter(4,[low high]./fn,'bandpass');
    afdsMax=filtfilt(b,a,afdsM);
    [yfds,~]=envelope(afdsMax,50,'peak');
    end
    close all;
end

    
   figure
plot(n,aedcMax,n,yedc);
legend("EMG Raw Data","EMG envelope");
xlabel("Index");
ylabel("Electric Potential (Volt)");
title("Raw EMG data and envelope for EDC muscle"); 

figure
plot(n,afdpMax,n,yfdp);
legend("EMG Raw Data","EMG envelope");
xlabel("Index");
ylabel("Electric Potential (Volt)");
title("Raw EMG data and envelope for FDP muscle");

figure
plot(n,afdsMax,n,yfds);
legend("EMG Raw Data","EMG envelope");
xlabel("Index");
ylabel("Electric Potential (Volt)");
title("Raw EMG data and envelope for FDS muscle");
    
    
nl=input("Give the starting index: ");
nh=input("Give the ending index: ");

aEDCMax=mean(yedc(nl:nh));
aFDPMax=mean(yfdp(nl:nh));
aFDSMax=mean(yfds(nl:nh));


disp("EDC Max voltage "+aEDCMax+" +/- "+std(yedc(nl:nh)));
disp("FDP Max voltage "+aFDPMax+" +/- "+std(yfdp(nl:nh)));
disp("FDS Max voltage "+aFDSMax+" +/- "+std(yfds(nl:nh)));