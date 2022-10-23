%Luo-Rudy membrane model
clear all

Gna = 23.0;
Gsi = 0.09;
Gk = 0.282;
Gk1 = 0.6047;
Gkp = 0.0183;
Gb = 0.03921;

Ena = 54.4;
Ek = -77.0;
Ek1 = -87.26; %???
Ekp = Ek1;
Eb = 59.87;

Ccai = 2.0e-2;%1.0e-9;

%intial membrane potential
Vm=-84.0;

dt=0.01;
t=40000;
stimulus=-1000;
start_stim=1000;

%initial gating parameters
am = 0.32*(Vm + 47.13)/(1.0 - exp(-0.1*(Vm + 47.13)));
bm = 0.08*exp(-Vm/11.0);

ah = 0.135*exp(-(80.0 + Vm)/6.8);
bh = 3.56*exp(0.079*Vm) + (3.1e5)*exp(0.35*Vm);

aj = ((-1.2714e5)*exp(0.2444*Vm) - (3.474e-5)*exp(-0.04391*Vm))*(Vm + 37.78)/(1.0 + exp(0.311*(Vm + 79.23)));
bj = 0.1212*exp(-0.01052*Vm)/(1.0 + exp(-0.1378*(Vm + 40.14)));

ad = 0.095*exp(-0.01*(Vm - 5.0))/(1.0 + exp(-0.072*(Vm - 5.0)));
bd = 0.07*exp(-0.017*(Vm + 44.0))/(1.0 + exp(0.05*(Vm + 44.0)));
    
af = 0.012*exp(-0.008*(Vm + 28))/(1.0 + exp(0.15*(Vm + 28.0)));
bf = 0.0065*exp(-0.02*(Vm + 30.0))/(1.0 + exp(-0.2*(Vm + 30.0)));

ax = 0.0005*exp(0.083*(Vm + 50.0))/(1.0 + exp(0.057*(Vm + 50.0)));
bx = 0.0013*exp(-0.06*(Vm + 20.0))/(1.0 + exp(-0.04*(Vm + 20.0)));

ak1 = 1.02/(1.0 + exp(0.2385*(Vm - Ek1 - 59.215)));
bk1 = (0.49124*exp(0.08032*(Vm - Ek1 + 5.476)) + exp(0.06175*(Vm - Ek1 - 594.31)))/(1.0 + exp(-0.5143*(Vm - Ek1 + 4.753)));
    
    
m = am/(am+bm);
h = ah/(ah+bh);
j = aj/(aj+bj);
d = ad/(ad+bd);
f = af/(af+bf);
X = ax/(ax+bx);
K1 = ak1/(ak1+bk1);

for i=1:t,
    
    if (i<=10)%start_stim)
        stim=stimulus;
    else
        stim=0;
    end
    
    %fast sodium current
    if (Vm >= -40.0)
        ah = 0.0;
        aj = 0.0;
        bh = 1.0/(0.13*(1.0 + exp((Vm + 10.66)/(-11.1))));
        bj = 0.3*exp((-2.535e-7)*Vm)/(1.0 + exp(-0.1*(Vm + 32)));
    end
      
    if (Vm < -40.0)
        ah = 0.135*exp(-(80.0 + Vm)/6.8);
        bh = 3.56*exp(0.079*Vm) + (3.1e5)*exp(0.35*Vm);
        aj = ((-1.2714e5)*exp(0.2444*Vm) - (3.474e-5)*exp(-0.04391*Vm))*(Vm + 37.78)/(1.0 + exp(0.311*(Vm + 79.23)));
        bj = 0.1212*exp(-0.01052*Vm)/(1.0 + exp(-0.1378*(Vm + 40.14)));
    end
      
    am = 0.32*(Vm + 47.13)/(1.0 - exp(-0.1*(Vm + 47.13)));
    bm = 0.08*exp(-Vm/11.0);

    %slow inward current
    %cai = Ccai;
    Esi = 7.7 - 13.0287*log(Ccai);
    ad = 0.095*exp(-0.01*(Vm - 5.0))/(1.0 + exp(-0.072*(Vm - 5.0)));
    bd = 0.07*exp(-0.017*(Vm + 44.0))/(1.0 + exp(0.05*(Vm + 44.0)));
    af = 0.012*exp(-0.008*(Vm + 28))/(1.0 + exp(0.15*(Vm + 28.0)));
    bf = 0.0065*exp(-0.02*(Vm + 30.0))/(1.0 + exp(-0.2*(Vm + 30.0)));
    

    %outward currents
    %time-dependent potassium current
    if (Vm > -100.0)
        Xi = 2.837*(exp(0.04*(Vm + 77.0)) - 1)/((Vm + 77.0)*exp(0.04*(Vm + 35.0)));
    end
    
    if (Vm <= -100)
       Xi = 1.0;
    end
      
    ax = 0.0005*exp(0.083*(Vm + 50.0))/(1.0 + exp(0.057*(Vm + 50.0)));
    bx = 0.0013*exp(-0.06*(Vm + 20.0))/(1.0 + exp(-0.04*(Vm + 20.0)));

    %time-independent potassium current
    ak1 = 1.02/(1.0 + exp(0.2385*(Vm - Ek1 - 59.215)));
    bk1 = (0.49124*exp(0.08032*(Vm - Ek1 + 5.476)) + exp(0.06175*(Vm - Ek1 - 594.31)))/(1.0 + exp(-0.5143*(Vm - Ek1 + 4.753)));
    K1  = ak1/(ak1+bk1);

    %plateau potassium current
    Kp = 1.0/(1.0 + exp((7.488 - Vm)/5.98));

    Ina = Gna*m*m*m*h*j*(Vm - Ena);
    Isi = Gsi*d*f*(Vm - Esi);
    Ik = Gk*X*Xi*(Vm - Ek);
    Ik1 = Gk1*K1*(Vm - Ek1);
    Ikp = Gkp*Kp*(Vm - Ekp);
    Ib = Gb*(Vm + Eb);

    Vm = Vm - dt*( Ina + Isi + Ik + Ik1 + Ikp + Ib + stim);

    m = m + dt*(am*(1.0-m) - bm*m);
    h = h + dt*(ah*(1.0-h) - bh*h);
    j = j + dt*(aj*(1.0-j) - bj*j);
    d = d + dt*(ad*(1.0-d) - bd*d);
    f = f + dt*(af*(1.0-f) - bf*f);
    X = X + dt*(ax*(1.0-X) - bx*X);

    %calcium uptake
    Ccai = Ccai + dt*((-10e-4)*Isi + 0.07*((10e-4) - Ccai));
      
    V(i)=Vm;
    stim=0;
end

plot(V)
grid on
