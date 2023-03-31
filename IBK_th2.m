function dydt=IBK_th2(M2,M3,L2,EDC,FDP,FDS,g,r,I2,theq,grav,t,y)
dydt=zeros(2,1);
%syms t2(t) t 
b2=r(1);
k2=r(2);
rEDC=r(3);
rFDP=r(4);
rFDS=r(5);

if isempty(grav)==0
    %Equation when gravity is in the plane of movement
    % eqn=rEDC*ppval(EDC,tim)+rFDP*ppval(FDP,tim)+rFDS*ppval(FDS,tim)==(I2*diff(t2(t),t,2)+b2*diff(t2(t),t)+k2*(t2(t)-theq)-0.5*(M2+2*M3)*g*L2*cos(t2(t)));
    dydt(1)=y(2);
    dydt(2)=(rEDC*ppval(EDC,t)+rFDP*ppval(FDP,t)+rFDS*ppval(FDS,t)-(b2*y(2)+k2*(y(1)-theq)-0.5*(M2+2*M3)*g*L2*cos(y(1))))/I2;
else
    %Equation when gravity is not in the plane of movement  
    %eqn=rEDC*ppval(EDC,tim)+rFDP*ppval(FDP,tim)+rFDS*ppval(FDS,tim)==I2*diff(t2(t),t,2)+b2*diff(t2(t),t)+k2*(t2(t)-theq);
    dydt(1)=y(2);
    dydt(2)=(rEDC*ppval(EDC,t)+rFDP*ppval(FDP,t)+rFDS*ppval(FDS,t)-(b2*y(2)+k2*(y(1)-theq)))/I2;
end 
 %[De,~]=odeToVectorField(eqn);  
 %dydt=matlabFunction(De,'Vars',{'t','Y'});
end