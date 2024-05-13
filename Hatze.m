function Fv=Hatze(fl,a)
syms h
fv=matlabFunction(0.1433/(0.1074+exp(-1.409*sinh(3.2*h+1.6)))); %Hatze function 
Fm=max(fl);
Af=0.25;%Value taken from Thelen 2003
Flen=1.4;%Value taken from Thelen 2003
Fv=zeros(1,length(fl));
b=0;
for i=1:length(fl)
    if Fm<=a(i)*fl(i)
        b=a(i)*fl(i)+Fm/Af;
    else
        b=(2+2/Af)*(a(i)*fl(i)*Flen-Fm)/(Flen-1);
    end
    v=(0.25+0.75*a(i))*(Fm-a(i)*fl(i))/b;
    Fv(i)=feval(fv,v);
end
end

