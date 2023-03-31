function dydt=IBK_Thumb(Icmc,Imcp,Iip,Ia,i,t,r,theq,y)

dydt=zeros(2,1);
B=r(1);
K=r(2);


dydt(1)=y(2);
    if i==1
dydt(2)=-(B*y(2)+K*(y(1)-theq))/Icmc;
    elseif i==2
     dydt(2)=-(B*y(2)+K*(y(1)-theq))/Imcp; 
    elseif i==3
        dydt(2)=-(B*y(2)+K*(y(1)-theq))/Iip; 
    else
        dydt(2)=-(B*y(2)+K*(y(1)-theq))/Ia;
    end
end
    


