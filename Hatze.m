function Fv=Hatze(Fa,a,F,th,in)
%This function computes the normalised Force-Velocity interaction for each
%muscle using the method in Thelen 2003 Adjustment of muscle mechanics
%model parameters to simulate dynamic contractions in older adults.



syms h x
fl1=matlabFunction(exp(-(h-1)^2/0.45)); %Active force-length scale factor from Thelen 2003.

%Fits are taken from the Data in OpenSim model
if in==2
    if F=='EDC_MCP'
        Lm=feval(matlabFunction(0.1224*x+0.8781),th);
    elseif F=='FDP_MCP'
        Lm=feval(matlabFunction(-0.1308*x+1.001),th);
    elseif F=='FDS_MCP'
        Lm=feval(matlabFunction(-0.1365*x+1.051),th);
    elseif F=='EDC_PIP'
        Lm=feval(matlabFunction(-0.004322*x^2+0.0627*x+0.8814),th);
    elseif F=='FDP_PIP'
        Lm=feval(matlabFunction(-0.1012*x+1.003),th);
    elseif F=='FDS_PIP'
        Lm=feval(matlabFunction(0.009727*x^2-0.1034*x+1.052),th);
    elseif F=='EDC_DIP'
        Lm=feval(matlabFunction(-0.01155*x^2+0.05426*x+0.8797),th);
    elseif F=='FDP_DIP'
        Lm=feval(matlabFunction(-0.008109*x^2-0.04448*x+1.001),th);
    end
elseif in==3
    if F=='EDC_MCP'
        Lm=feval(matlabFunction(-0.01752*x^3+0.03269*x^2+0.1035*x+0.9992),th);
    elseif F=='FDP_MCP'
        Lm=feval(matlabFunction(-0.01717*x^2-0.09995*x+0.9989),th);
    elseif F=='FDS_MCP'
        Lm=feval(matlabFunction(-0.02252*x^2-0.1262*x+0.999),th);
    elseif F=='EDC_PIP'
        Lm=feval(matlabFunction(0.06029*x+1.003),th);
    elseif F=='FDP_PIP'
        Lm=feval(matlabFunction(-0.009293*x^2-0.08988*x+1.003),th);
    elseif F=='FDS_PIP'
        Lm=feval(matlabFunction(-0.01212*x^2-0.105*x+1.002),th);
    elseif F=='EDC_DIP'
        Lm=feval(matlabFunction(-0.003304*x^2+0.04484*x+1.001),th);
    elseif F=='FDP_DIP'
        Lm=feval(matlabFunction(-0.008319*x^2-0.05544*x+1.001),th);
    end
elseif in==4
    if F=='EDC_MCP'
        Lm=feval(matlabFunction(0.009867*x^4-0.02677*x^3+0.02975*x^2+0.08636*x+0.9998),th);
    elseif F=='FDP_MCP'
        Lm=feval(matlabFunction(-0.01169*x^2-0.1102*x+0.9996),th);
    elseif F=='FDS_MCP'
        Lm=feval(matlabFunction(-0.01658*x^2-0.1263*x+0.999),th);
    elseif F=='EDC_PIP'
        Lm=feval(matlabFunction(0.01199*x^3-0.03524*x^2+0.07608*x+0.9992),th);
    elseif F=='FDP_PIP'
        Lm=feval(matlabFunction(0.007089*x^3-0.0264*x^2-0.05486*x+1.001),th);
    elseif F=='FDS_PIP'
        Lm=feval(matlabFunction(0.008742*x^3-0.02808*x^2-0.03884*x+1),th);
    elseif F=='EDC_DIP'
        Lm=feval(matlabFunction(-0.01294*x^2+0.05507*x+0.9999),th);
    elseif F=='FDP_DIP'
        Lm=feval(matlabFunction(-0.008374*x^2-0.03833*x+1.001),th);
    end
elseif in==5
    if F=='EDC_MCP'
        Lm=feval(matlabFunction(0.004721*x^2+0.09731*x+0.9978),th);
    elseif F=='FDP_MCP'
        Lm=feval(matlabFunction(-0.00884*x^2-0.08256*x+1.002),th);
    elseif F=='FDS_MCP'
        Lm=feval(matlabFunction(-0.03104*x^2-0.1325*x+0.999),th);
    elseif F=='EDC_PIP'
        Lm=feval(matlabFunction(0.01807*x^3-0.04932*x^2+0.09257*x+0.9992),th);
    elseif F=='FDP_PIP'
        Lm=feval(matlabFunction(0.01998*x^3-0.04947*x^2-0.06765*x+0.9996),th);
    elseif F=='FDS_PIP'
        Lm=feval(matlabFunction(-0.01528*x^2-0.1082*x+1.002),th);
    elseif F=='EDC_DIP'
        Lm=feval(matlabFunction(-0.01914*x^2+0.03838*x+0.999),th);
    elseif F=='FDP_DIP'
        Lm=feval(matlabFunction(0.008211*x^3-0.02538*x^2-0.04179*x+1.001),th);
    end
end

fl=feval(fl1,Lm);

fv=matlabFunction(0.1433/(0.1074+exp(-1.409*sinh(3.2*h+1.6)))); %Hatze function from A 3-D dynamic model 
%of human finger for studying free movements

Fm=Fa./max(Fa);
Af=0.25;%Value taken from Thelen 2003
Flen=1.4;%Value taken from Thelen 2003
Fv=zeros(1,length(Fa));
b=0;
for i=1:length(Fa)
    if Fm(i)<=a(i)*fl(i)
        b=a(i)*fl(i)+Fm(i)/Af;
    else
        b=(2+2/Af)*(a(i)*fl(i)*Flen-Fm(i))/(Flen-1);
    end
    v=(0.25+0.75*a(i))*(Fm(i)-a(i)*fl(i))/b;
    Fv(i)=feval(fv,v);
end
end


