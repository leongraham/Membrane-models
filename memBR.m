%Beeler-Reuter membrane model
clear all

%membrane parameters
Gna = 4%2.5%4.0;
Gnc = 0.003%0.0077%0.003;
Ena = 50;
Gs = 0.09;%0.09

%intial membrane potential
Vm=-84.0;
Cm=1.0;
dt=0.01;
t=2000;
stimulus=-500; %-500
start_stim=20000;

%initial conditions for gating variables
ax1 = (0.0005*exp(0.083*(Vm + 50.0)))/(exp(0.057*(Vm + 50.0)) + 1.0);
bx1 = (0.0013*exp(-0.06*(Vm + 20.0)))/(exp(-0.04*(Vm + 20.0)) + 1.0);
am = (-1.0*(Vm + 47.0))/(exp(-0.1*(Vm + 47.0)) - 1.0);
bm = (40.0*exp(-0.056*(Vm + 72.0)));
ah = (0.126*exp(-0.25*(Vm + 77.0)));
bh = (1.7)/(exp(-0.082*(Vm + 22.5)) + 1.0);
aj = (0.055*exp(-0.25*(Vm + 78.0)))/(exp(-0.2*(Vm + 78.0)) + 1.0);
bj = (0.3)/(exp(-0.1*(Vm + 32.0)) + 1.0);
ad = (0.095*exp(-0.01*(Vm -5.0)))/(exp(-0.072*(Vm - 5.0)) + 1.0);
bd = (0.07*exp(-0.017*(Vm + 44.0)))/(exp(0.05*(Vm + 44.0)) + 1.0);
af = (0.012*exp(-0.008*(Vm + 28.0)))/(exp(0.15*(Vm + 28.0)) + 1.0);
bf = (0.0065*exp(-0.02*(Vm + 30.0)))/(exp(-0.2*(Vm + 30.0)) + 1.0);

%steady state
m = am/(am + bm);
h = ah/(ah + bh);
j = aj/(aj + bj);
d = ad/(ad + bd);
f = af/(af + bf);
x1 = ax1/(ax1 + bx1);

%initial resting intracellular Ca concentration
cai = 2e-4;


for i=1:t,
    
    if (i<=5)%start_stim)
        stim=stimulus;
    else
        stim=0;
    end
    
    ax1 = (0.0005*exp(0.083*(Vm + 50.0)))/(exp(0.057*(Vm + 50.0)) + 1.0);
    bx1 = (0.0013*exp(-0.06*(Vm + 20.0)))/(exp(-0.04*(Vm + 20.0)) + 1.0);
    am = (-1.0*(Vm + 47.0))/(exp(-0.1*(Vm + 47.0)) - 1.0);%(-1.0*(Vm + 47.0))/(exp(-0.1*(Vm + 47.0)) - 1.0);
    bm = (40.0*exp(-0.056*(Vm + 72.0)));%(40.0*exp(-0.056*(Vm + 72.0)));
    ah = (0.126*exp(-0.25*(Vm + 77.0)));
    bh = (1.7)/(exp(-0.082*(Vm + 22.5)) + 1.0);
    aj = (0.055*exp(-0.25*(Vm + 78.0)))/(exp(-0.2*(Vm + 78.0)) + 1.0);
    bj = (0.3)/(exp(-0.1*(Vm + 32.0)) + 1.0);
    ad = (0.095*exp(-0.01*(Vm -5.0)))/(exp(-0.072*(Vm - 5.0)) + 1.0);
    bd = (0.07*exp(-0.017*(Vm + 44.0)))/(exp(0.05*(Vm + 44.0)) + 1.0);
    af = (0.012*exp(-0.008*(Vm + 28.0)))/(exp(0.15*(Vm + 28.0)) + 1.0);
    bf = (0.0065*exp(-0.02*(Vm + 30.0)))/(exp(-0.2*(Vm + 30.0)) + 1.0);

    Es = -82.3 - 13.0287*log(cai);
    Is = Gs*d*f*(Vm - Es);
    Ik1 = 0.35*(4.0*(exp(0.04*(Vm + 85.0)) - 1.0)/(exp(0.08*(Vm + 53.0)) + exp(0.04*(Vm + 53.0))) + 0.2*(Vm + 23.0)/(1.0 -exp(-0.04*(Vm + 23.0))));
    Ix1 = x1*0.8*(exp(0.04*(Vm + 77.0)) - 1.0)/exp(0.04*(Vm + 35.0));
    Ina = (Gna*m*m*m*h*j + Gnc)*(Vm - Ena);
    
    Vm = Vm - dt*(1/Cm)*(Ik1 + Ix1 + Ina + Is + stim);

    m = m + dt*(am*(1.0-m) - bm*m);
    h = h + dt*(ah*(1.0-h) - bh*h);
    j = j + dt*(aj*(1.0-j) - bj*j);
    d = d + dt*(ad*(1.0-d) - bd*d);
    f = f + dt*(af*(1.0-f) - bf*f);
    x1 = x1 + dt*(ax1*(1.0-x1) - bx1*x1);

    %calcium uptake
    cai = cai + dt*((-10e-7)*Is + 0.07*((10e-7) - cai));
      
    V(i)=Vm;
end

plot(V)
grid on
