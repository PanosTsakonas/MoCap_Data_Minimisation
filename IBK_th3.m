function dydt=IBK_th3(M3,L3,EDC,FDP,g,r,I3,theq,grav,t,y)

dydt=zeros(2,1);
%syms  t3(t) t 
b3=r(1);
k3=r(2);
rEDC=r(3);
rFDP=r(4);

if isempty(grav)==0
  %Equation when gravity is in the plane of movement  
  %eqn=rEDC*ppval(EDC,tim)+rFDP*ppval(FDP,tim)==(I3*diff(t3(t),t,2)+b3*diff(t3(t),t)+k3*(t3(t)-theq)-M3*g*L3*cos(t3(t))/2);
  dydt(1)=y(2);
  dydt(2)=(rEDC*ppval(EDC,t)+rFDP*ppval(FDP,t)-(b3*y(2)+k3*(y(1)-theq)-M3*g*L3*cos(y(1))/2))/I3;
else
    %Equation when gravity is not in the plane of movement 
    %eqn=rEDC*ppval(EDC,tim)+rFDP*ppval(FDP,tim)==I3*diff(t3(t),t,2)+b3*diff(t3(t),t)+k3*(t3(t)-theq);
   dydt(1)=y(2);
   dydt(2)=(rEDC*ppval(EDC,t)+rFDP*ppval(FDP,t)-(b3*y(2)+k3*(y(1)-theq)))/I3;
end  
    %[De,~]=odeToVectorField(eqn);
    %dydt=matlabFunction(De,'Vars',{'t','Y'});
end