function Res=Search_For_Global(M1,M2,M3,L1,L2,L3,EDC_MCP,FDP_MCP,FDS_MCP,g,r,t,init,th,i,I,theq,grav,in,Boot,Itter)

%Set the minimisation options
options=optimoptions('fmincon','FunValCheck','on','UseParallel',true);
R=zeros(Itter,length(r));

for j=1:Itter
    R=r+rand(1,5).*r;
    
    if i<3
        Res(j,:)=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EDC_MCP,FDP_MCP,FDS_MCP,g,r,t,init,th,i,I,theq,grav,in,Boot),R,[],[],[],[],[r(1) r(2) 0.1 0.1 0.1],[2 10 2 2 2],[],options);
    else
        Res(j,:)=fmincon(@(r) myobj(M1,M2,M3,L1,L2,L3,EDC_MCP,FDP_MCP,FDS_MCP,g,r,t,init,th,i,I,theq,grav,in,Boot),R,[],[],[],[],[r(1) r(2) 0.1 0.1 ],[2 10 2 2 ],[],options);
    end
end
end

