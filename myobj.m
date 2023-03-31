function [obj,th_exp]=myobj(M1,M2,M3,L1,L2,L3,EXT,FL1,FL2,g,r,tim,init,th_G,i,I,theq,grav,in)
%This functions minimises the sum of the square differences between the IBK
%approximation for each degree of freedom and the filtered angular data obtained from the
%motion capture laboratory.
if in~=1
if i==1
    disp("Minimising the parameters for the movement of the proximal segment");
[~,Y]=ode45(@(t,y) IBK_th1(M1,M2,M3,L1,EXT,FL1,FL2,g,r,I(1),theq,grav,t,y),tim,[init(1) init(2)]);
%[~,Y]=ode45(IBK_th1(M1,M2,M3,L1,EXT,FL1,FL2,g,r,I(1),theq,grav,[],[]),tim,[init(1) init(2)]);
th_exp=Y(:,1);

elseif i==2
    disp("Minimising the parameters for the movement of the middle segment");
    [~,Y]=ode45(@(t,y) IBK_th2(M2,M3,L2,EXT,FL1,FL2,g,r,I(2),theq,grav,t,y),tim,[init(3) init(4)]);
    %[~,Y]=ode45(IBK_th2(M2,M3,L2,EXT,FL1,FL2,g,r,I(2),theq,grav,[],[]),tim,[init(1) init(2)]);
    th_exp=Y(:,1);
    

else
    disp("Minimising the parameters for the movement of the distal segment");
    [~,Y]=ode45(@(t,y) IBK_th3(M3,L3,EXT,FL1,g,r,I(3),theq,grav,t,y),tim,[init(5) init(6)]);
    %[~,Y]=ode45(IBK_th3(M3,L3,EXT,FL1,g,r,I(3),theq,grav,[],[]),tim,[init(1) init(2)]);
th_exp=Y(:,1);
end
else
    
    if i==1
        disp("Minimising the parameters for the CMC joint movement of the thumb");
    [~,Y]=ode45(@(t,y) IBK_Thumb(M1,M2,M3,L1,L2,L3,I(1),I(2),I(3),EXT,FL1,i,g,r,theq,grav,t,y),tim,[init(1) init(2)]);
    elseif i==2
        disp("Minimising the parameters for the MCP joint movement of the thumb");
     [~,Y]=ode45(@(t,y) IBK_Thumb(M1,M2,M3,L1,L2,L3,I(1),I(2),I(3),EXT,FL1,i,g,r,theq,grav,t,y),tim,[init(3) init(4)]); 
    else
        disp("Minimising the parameters for the IP joint movement of the thumb");
        [~,Y]=ode45(@(t,y) IBK_Thumb(M1,M2,M3,L1,L2,L3,I(1),I(2),I(3),EXT,FL1,i,g,r,theq,grav,t,y),tim,[init(5) init(6)]);
    end
    
    th_exp=Y(:,1);
end

%Least squares
obj=sum(sum((th_exp-th_G).^2));

%Root mean square error
%obj=sqrt(mean((th_exp-th_G).^2));
end


