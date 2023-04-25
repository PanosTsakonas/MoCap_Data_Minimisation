function dadt=activation(t,a,u)

dadt=0;

% Tact and Tdeact are taken from the "Real-time simulation of hand motion for prosthesis control". The initial
%muscle activation since there is no motion before the experiment is assumed to be equal to zero. The differential
%equation is also taken from the same paper.

Tact=15*10^-3;
Tdeact=50*10^-3;

dadt=(ppval(u,t)./Tact+(1-ppval(u,t))./Tdeact).*(ppval(u,t)-a);

end
