function [EXT_P,FL1_P,FL2_P, EXT_M, FL1_M, FL2_M, EXT_D, FL_D]=Moments(Res,tim,N,IBK,Lag,in)
%This function calculates the spline interpolation to the muscle
%moments at each degree of freedom.
if in~=1
syms t
D=["EDC MCP","FDP MCP","FDS MCP", "EDC PIP", "FDP PIP", "FDS PIP", "EDC DIP", "FDP DIP"];
for i=1:8
    if isempty(IBK)==0
    P=spline(tim,Res(:,i));
    pp=mkpp(P.breaks,P.coefs);
    end
    
    if isempty(Lag)==0
   P=polyfit(tim,Res(:,i),N);
    pp=0;
   for j=0:N
        pp=pp+P(j+1)*t^(N-j);
    end
    end
    if i==1
    EXT_P=pp; 
    elseif i==2
FL1_P=pp;
elseif i==3
FL2_P=pp;
    elseif i==4
EXT_M=pp;
    elseif i==5
FL1_M=pp;
    elseif i==6
FL2_M=pp;
    elseif i==7
EXT_D=pp;
    else
FL_D=pp;
    end
    if isempty(Lag)==0
    Rs=1-sum((feval(matlabFunction(pp,tim))-Res(:,i)).^2)/sum((Res(:,i)-mean(Res(:,i))).^2);
    figure
    plot(tim,feval(matlabFunction(pp,tim)),tim,Res(:,i),'x');
    xlabel("Time (s)");
    ylabel("Moment of "+D(i)+" (N*m)");
    title("Polynomial fit of "+D(i)+" muslce with R^2 "+num2str(Rs));
    end
    
    if isempty(IBK)==0
     figure
    plot(tim,ppval(pp,tim),tim,Res(:,i),'x');
    xlabel("Time (s)");
    ylabel("Moment of "+D(i)+" (N*m)");
    title("Spline fit of "+D(i)+" muslce");
    end   
        
clear Rs P pp; 
end
else
    disp("I am doing the thumb moment fits");
    syms t
    D=["EPL CMC" "FPL CMC" "EPL MCP" "FPL MCP" "EPL IP" "FPL IP"];
    for i=1:6
      
      if isempty(IBK)==0
    P=spline(tim,Res(:,i));
    pp=mkpp(P.breaks,P.coefs);
      end
    
    if isempty(Lag)==0
   P=polyfit(tim,Res(:,i),N);
    pp=0;
   for j=0:N
        pp=pp+P(j+1)*t^(N-j);
    end
    end
      if i==1
          EXT_P=pp;
      elseif i==2
          FL1_P=pp;
      elseif i==3
          EXT_M=pp;
      elseif i==4
          FL1_M=pp;
      elseif i==5
          EXT_D=pp;
      else
          FL_D=pp;
      end
      
    if isempty(Lag)==0
    Rs=1-sum((feval(matlabFunction(pp,tim))-Res(:,i)).^2)/sum((Res(:,i)-mean(Res(:,i))).^2);
    figure
    plot(tim,feval(matlabFunction(pp,tim)),tim,Res(:,i),'x');
    xlabel("Time (s)");
    ylabel("Moment of "+D(i)+" (N*m)");
    title("Polynomial fit of "+D(i)+" muslce with R^2 "+num2str(Rs));
    end
    
    if isempty(IBK)==0
     figure
    plot(tim,ppval(pp,tim),tim,Res(:,i),'x');
    xlabel("Time (s)");
    ylabel("Moment of "+D(i)+" (N*m)");
    title("Spline fit of "+D(i)+" muslce");
    end   
        
clear Rs P pp;   
    
    end
FL2_P=[];
FL2_M=[];
end
end


