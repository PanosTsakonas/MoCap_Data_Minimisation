function save_resultsM(Min1,Min2,Min3,mn1,mn2,mn3,in,Par,sd_mcp,sd_pip,sd_dip,trial,rsq,rmse,cyl,wn_mcp,wn_pip,wn_dip)

param={'B (Nms/rad)','sdB (Nms/rad)','K (Nm/rad)','sdK (Nm/rad)','rho_1','sdrho_1','rho_2','sdrho_2','rho_3','sdrho_3','Comments'};
D=["Thumb","Index","Middle","Ring","Little"];
fing=transpose({D(in),'MCP','PIP','DIP'});
m={Min1(1) sd_mcp(1) Min1(2) sd_mcp(2) Min1(3) sd_mcp(3) Min1(4) sd_mcp(4) Min1(5) sd_mcp(5) "R^2 :"+rsq(1)+" RMSE: "+rmse(1)+" deg"};
p={Min2(1) sd_pip(1) Min2(2) sd_pip(2) Min2(3) sd_pip(3) Min2(4) sd_pip(4) Min2(5) sd_pip(5) "R^2 :"+rsq(2)+" RMSE: "+rmse(2)+" deg"};
d={Min3(1) sd_dip(1) Min3(2) sd_dip(2) Min3(3) sd_dip(3) Min3(4) sd_dip(4) "" ""  "R^2 :"+rsq(3)+" RMSE: "+rmse(3)+" deg"};
M=[{'MN1';mn1} {'MN2';mn2} {'MN3';mn3} {'MCP_cutoff (Hz)';wn_mcp} {'PIP_cutoff (Hz)';wn_pip} {'DIP_cutoff (Hz)';wn_dip}];
file="P"+Par+"\Cylindrical_"+cyl+"_"+trial+"_Par_"+Par+".xlsx";

if in==2
writecell(param,file,'Sheet',1,'Range','B1');
writecell(fing,file,'Sheet',1,'Range','A2');
writecell(m,file,'Sheet',1,'Range','B3');
writecell(p,file,'Sheet',1,'Range','B4');
writecell(d,file,'Sheet',1,'Range','B5');
writecell(M,file,'Sheet',1,'Range','N3');
elseif in==3   
writecell(fing,file,'Sheet',1,'Range','A7');
writecell(m,file,'Sheet',1,'Range','B8');
writecell(p,file,'Sheet',1,'Range','B9');
writecell(d,file,'Sheet',1,'Range','B10');
writecell(M,file,'Sheet',1,'Range','N8');
elseif in==4
writecell(fing,file,'Sheet',1,'Range','A12');
writecell(m,file,'Sheet',1,'Range','B13');
writecell(p,file,'Sheet',1,'Range','B14');
writecell(d,file,'Sheet',1,'Range','B15');
writecell(M,file,'Sheet',1,'Range','N13');
else
writecell(fing,file,'Sheet',1,'Range','A17');
writecell(m,file,'Sheet',1,'Range','B18');
writecell(p,file,'Sheet',1,'Range','B19');
writecell(d,file,'Sheet',1,'Range','B20');
writecell(M,file,'Sheet',1,'Range','N18');
end

    
