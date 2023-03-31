function dydt=IBK_th1(M1,M2,M3,L1,EDC,FDP,FDS,g,r,I1,theq,grav,t,y)
%syms t1(t) t
dydt=zeros(2,1);
B=r(1);
K=r(2);
rEDC=r(3);
rFDP=r(4);
rFDS=r(5);

if isempty(grav)==0
    %Equation when gravity is in the plane of movement
    dydt(1)=y(2);
    dydt(2)=(rEDC*ppval(EDC,t)+rFDP*ppval(FDP,t)+rFDS*ppval(FDS,t))/I1-(B*y(2)+K*(y(1)-theq)-0.5*(M1+2*(M2+M3))*g*L1*cos(y(1)))/I1;
     %eqn=rEDC*ppval(EDC,t)+rFDP*ppval(FDP,t)+rFDS*ppval(FDS,t)==I1*diff(t1(t),t,t)+B*diff(t1(t),t)+K*(t1(t)-theq)-0.5*(M1+2*(M2+M3))*g*L1*cos(t1(t));
else
     %Equation when gravity is not in the plane of movement
     dydt(1)=y(2);
     dydt(2)=(rEDC*ppval(EDC,t)+rFDP*ppval(FDP,t)+rFDS*ppval(FDS,t))/I1-(B*y(2)+K*(y(1)-theq))/I1;
     %eqn=rEDC*EDC+rFDP*FDP+rFDS*FDS==I1*diff(t1(t),t,2)+B*diff(t1(t),t)+K*(t1(t)-theq);     
end

    %[De,~]=odeToVectorField(eqn);
    %dydt=matlabFunction(De,'Vars',{'t','Y'});
end

    

