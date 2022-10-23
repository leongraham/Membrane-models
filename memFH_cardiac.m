clear all
FH_Vt=5.0;
FH_Vp=100.0;
C1 = 0.2%0.5;
C2 = 0.005%0.02;
C3 = 0.015;
C4 = 0.005;

dt=0.1;
Vm=6
w=0;

for i=1:1000,
    Vm = Vm + dt*(C1*Vm*(Vm/FH_Vt-1.0)*(1.0-Vm/FH_Vp) - C2*Vm*w);
    w = w + dt*C3*(Vm - C4*w);
    V(i)=Vm;
    W(i)=w;
end

plot(V)
grid on
%hold on
%plot(W,'r')