function outfunc=newzsmrecyclefunc(p,t,f,bh)
clear all;
clf;
clc;
%% legend
%c(1)= cemeoh, c(2)=cedme, c(3)=ceh2o, c(4)=cec2h4, c(5)=cec3h6,
%c(6)=cec4h8, c(7)=cec5h10, c(8)=cec6h12, c(9)=cecx, c(10)= cec, c(11)=
%cech4
%c(12:22) are bubble phase molar concentrations (kmol/m3)

%r1=-k2*c(1)*z-k15*c(1)*z-req1+2*req2; 
%r2=-k3*c(2)*z+0.5*req1-req2;
%r3=k2*c(1)*z+k3*c(2)*z+k15*c(1)*z+0.5*req1-req2;
%r4=5*k4*c(9)*z-k5*c(4)*z;
%r5=(10/3)*k6*c(9)*z-k7*c(5)*z;
%r6=2.5*k8*c(9)*z-k9*c(6)*z;
%r7=2*k10*c(9)*z-k11*c(7)*z;
%r8=1.666*k12*c(9)*z-k13*c(8)*z;
%r9=0.1*k2*c(1)*z+0.2*k3*c(2)*z-k4*c(9)*z+0.2*k5*c(4)*z-k6*c(9)*z+0.3*k7*c(5)*z-k8*c(9)*z+0.4*k9*c(6)*z-k10*c(9)*z...
    %+0.5*k11*c(7)*z-k12*c(9)*z+0.6*k13*c(8)*z;
%r10=0.5*k15*c(1)*z;
%r11=0.5*k15*c(1)*z;


%% mass balance to be solved
function F=massbalance(c)
z=1/(1+Kw*c(3));
%eqm rxn
%Keqb=21400/Keq;
req1=(21400*((c(1))^2))/((1+15000*c(1)+47000*c(3))^2);%21400*c(1)^2;
req2=((21400/4.369)*c(2)*c(3))/((1+15000*c(1)+47000*c(3))^2);%Keqb*c(2)*c(3);
%rxn rates 
      r=[-k2*c(1)*z-k15*c(1)*z-req1+2*req2,...
         -k3*c(2)*z+0.5*req1-req2,...
         k2*c(1)*z+k3*c(2)*z+k15*c(1)*z+0.5*req1-req2,...
        5*k4*c(9)*z-k5*c(4)*z,...
         (10/3)*k6*c(9)*z-k7*c(5)*z,...
         2.5*k8*c(9)*z-k9*c(6)*z,...
     	2*k10*c(9)*z-k11*c(7)*z,...
        1.666*k12*c(9)*z-k13*c(8)*z,...
         0.1*k2*c(1)*z+0.2*k3*c(2)*z-k4*c(9)*z+0.2*k5*c(4)*z-k6*c(9)*z+0.3*k7*c(5)*z-k8*c(9)*z+0.4*k9*c(6)*z-k10*c(9)*z...
        +0.5*k11*c(7)*z-k12*c(9)*z+0.6*k13*c(8)*z,...
        0.5*k15*c(1)*z,...
        0.5*k15*c(1)*z];
%mass eqns
for i=1:11
    %odd=emulsion, even=bubble
    F(2*i-1)=-Vnew(1)*c(i)+V(a-1,1)*C(a-1,i)+Vnew(3)*K(a)*(c(i+11)-c(i))-Wn*sum(r)*(c(i)/(sum(c(1:11))-c(9)))+Wn*r(i); %emulsion
    F(2*i)=-Vnew(2)*c(i+11)+V(a-1,2)*C(a-1,i+11)-Vnew(3)*K(a)*(c(i+11)-c(i))+Wn*sum(r)*(c(i)/(sum(c(1:11))-c(9))); %bubble
  
end
%mass eqns for cx
F(17)=-Vnew(1)*c(9)+Vnew(1)*C(a-1,9)+Wn*r(9);
F(18)=c(20);

end
%% solving for evolflow and bvolflow 
function G = volumechange (Q)    
    G(1)= sum(mw.*conc(1:11))*Q(1)+sum(mw.*conc(12:22))*Q(2)-(sum(mw.*C(a-1,1:11))*V(a-1,1)+sum(mw.*C(a-1,12:22))*V(a-1,2));% mass balance between stages
end
%% constants 
p=evalin('base', 'p');
t=evalin('base', 't');
f=evalin('base', 'f');
bh=evalin('base', 'bh');

P=p; %3030000; %pressure constant of reactor, Pa, set at 30 atm to decrease volume flow
T=t;%500+273.15; %temperature constant of reactor, K
n=200; %no of cstr stages reduced from 500 to 200 to speed up excel generation
H=bh;%10; %m, total height of reactor, increase height, increase catalyst,
%more methanol reacted (but does not mean more C2C3, can be less)
beddensity=744;%744 from lit;%kg/m3 %using 200 instead for now to reduce catalyst used
rollp=1270;%kg/m^3, particledensity
rollg=15.44; %kg/m^3, gasdensity from hysys
miu =2.003e-005;%si units, viscosity from hysys
dp=308*10^(-6);%m, particle diameter, x3 from literature because to improve
Ff=0.004; %fine fraction
g=9.81;%m^2/s
D=0.4*10^(-4);%m^2/s, molar diffusion coefficient
emf=0.405;%voidage at min fluid

R=1.986*10^(-3);%kcal/mol.K

molarflow0=f; %3.73 kmol/s for 28500 m3/h
q0=molarflow0*8314*T/P; %m3/s
%q0=28500/3600; %m3/s %stage 0 volumetric flow rate from packed bed hysys =24380 m3/h, 36500 is w/o recycle for 880kt/a capacity
%molarflow0=P*q0/(8314*T); %kmol/s , assuming ideal gas, initial molar flow into reactor
molarflowmeth=0.9979*molarflow0; %kmol/s, mol fraction from problem statement
molarflowh20=0.0071*molarflow0; % kmol/s
%% values calculation
%hydrodynamic parameters that stay relatively constant throughout reactor
umf =(miu/(dp*rollg))*((33.7^2+0.0408*(dp^3*rollg*((rollp-rollg)*g)/miu^2))^(1/2)-33.7); % m/s
umb=umf*((2300*rollg^0.126*miu^0.523*exp(0.716*Ff))/(dp^0.8*g^0.934*(rollp-rollg)^0.934));%m/s
%test to compare initial umf and umb to that of final stage
lol=umf;
lol2=umb;
%endtest
%calculation of ut 
gravitationterm=g*(rollp-rollg);
ut=((((1.78*10^-2)*gravitationterm^2)/(rollg*miu))^0.3333)*...
    dp; %ms-1, (u0<ut for bubbling fluidization to work)

%calculation of internal diameter, volume of reactor, and weight of
%catalyst per stage
u0=8.5*umf; %ms-1, superficial velocity of inlet stage
area=q0/u0;%area of tube, m2
d=2*((area/pi)^0.5);%m internal diameter
reactorvolume=(pi/4)*d^2*H;%m3, bed volume
catalystweight=beddensity*reactorvolume; %kg
Wn=catalystweight/n;%kg 

%calculation of WHSV (massflow/catalystweight)
massflow0=(molarflowmeth*32.0429+molarflowh20*18.0159)*3600; %kg/h
WHSV=massflow0/catalystweight; %h-1

%calculation of stage 0 parameters
h0=0; %height of bed is 0 at inlet

db0=0.00376*(u0*100-umb*100)^2; %cm
dbm=0.652*(area*10000*(u0*100-umb*100))^0.4; %cm
db=(dbm/100)-((dbm/100)-(db0/100))*exp((-0.3*h0)/d); %m
ub=u0-umb+0.711*(g*db)^(0.5); %m
emulfrac=(u0-umb)/(ub-umb); %dimensionless, sigma is emulsion frac?
ubbar=ub;dbbar=db;
kbe=(umb/3)+((4*D*emf*ubbar)/(pi*dbbar))^(1/2);
Kbe=6*kbe/(db); %s-1

qb0=(1-emulfrac)*q0;%m3/s, bubble volume flow rate in stage 0
qe0=emulfrac*q0;%m3/s
Vb=(1-emulfrac)*reactorvolume;%m3, volume of bubble phase in reactor
Ve=emulfrac*reactorvolume;%m3
Vbn=Vb/n; %m3, volume of bubble phase per stage
Ven=Ve/n;
cemeohfeed=molarflowmeth/(qe0); %kmol/m3
ceh20feed=molarflowh20/(qe0); %kmol/m3
%% rxn constants
%constants
k2=1.7*exp(-6.28/(R*T));% kmol/kg s
k3=9*exp(-15.2/(R*T));
k4=2.7*10^(5)*exp(-20.2/(R*T));
k5=14.6*exp(-21/(R*T));
k6=2.76*10^(3)*exp(-11.46/(R*T));
k7=65*exp(-19.41/(R*T));
k8=87.1*exp(-6.94/(R*T));
k9=0.23*exp(-14.12/(R*T));
k10=0.73*exp(-0.69/(R*T));
k11=1.33*10^(-4)*exp(-3.92/(R*T));
k12=2.52*exp(-1.59/(R*T));
k13=13.5*exp(-5.23/(R*T));
k15=0.5*exp(-8.80/(R*T));
Kw=9520.7*exp(-12.1/(R*T));
%Keq=exp(-13.36+2835/T+1.675*log(T)-2.39*10^(-4)*T-2.1*10^(-7)*T^2);


%% molecular weight and viscosity of compounds
mw = [32.0429 , 46.0699, 18.0159, 28.0540, 42.0810, 56.1080, 70.1350 , 84.1620, 140.2700, 12.01, 16.04];
%mw(1)=meohmw, mw(2)=dmemw, mw(3)=h20mw, mw(4)=c2h4mw, mw(5)=c3h6mw, mw(6)=c4h8mw, mw(7)=c5h10mw,
%mw(8)=c6h12mw, mw(9)=cxmw, mw(10)=cmw, mw(11)=ch4mw %kg/kmol
vis= [0.00001949, 0.00002377, 0.00002841, 0.00002361, 0.00002221, 0.00002053, 0.00001891, 0.00001772, 0.00001386, 0, 0.00002326];
%vis in Pa.s at 500oC from hysys

%% creation of matrix required
%row 1= stage 0, row ntotal = stage ntotal -1
C = zeros(n,22); %Concentrations 
V = zeros(n, 4); %Volumetric flow rates(V1,V2) and Volume (V3,V4), emulsion bubble respectively 
efs = zeros(n, 1); %emulsionfraction, sigma is not bubble fraction apparently
K = zeros(n,1); %mass transf coefficient
h=zeros(n,1); %height of bed at each stage
dbs=zeros(n,1);%bubblediameter 
ubs=zeros(n,1);%bubblevelocity 
u0s=zeros(n,1);%superficial velocity 
uts=zeros(n,1); %terminal velocity
derp=zeros(n,1); %no. of iteration of while loop for each stage
%% setting stage 0 values to matrix
C(1,1)=cemeohfeed; %emulsion feed methanol conc
C(1,3)=ceh20feed; %emulsion feed h20 conc
V(1,1)=qe0;V(1,2)=qb0;V(1,3)=Ven;V(1,4)=Vbn;
efs(1)=emulfrac;
K(1)=Kbe;
h(1)=h0;
dbs(1)=db;
ubs(1)=ub;
u0s(1)=u0;
uts(1)=ut;
%% iterations
herp2=0;
C0=C(1,:);
V0=V(1,:);
Vr=1; %a guess to kickstart recycle loop
Vrnew=3; %a guess to kickstart recycle loop

while abs(Vr-Vrnew)>=0.001 %recycle loop
    herp2=herp2+1;
    Vr=Vrnew %for troubleshoot, place this at end
    if herp2>=2
        V(1,1)=V0(1)+Vr*emulfrac;
        V(1,2)=V0(2)+Vr*(1-emulfrac);
        V(1,3)= Ven;
        V(1,4)= Vbn;
        C(1,1:11)=(C0(1:11)*V0(1)+Cr(1:11)*Vr*emulfrac)/V(1);
        C(1,12:22)=(C0(12:22)*V0(2)+Cr(1:11)*Vr*(1-emulfrac))/V(2);
    else
        Vr=1;
    end
for a = 2:n; %no. of stages
       h(a)=h(a-1)+H/(n-1);
       Vnew= V(a-1,:);
       u0=0.5;      
       u0new=u0s(a-1);
       herp=0;
       
       while abs(u0new-u0) > 0.0012 %cant solve <500 iter if tolerance <0.001 at 300 cstr stgs
           %if fsolve stalls, result will be inaccurate, loosen tolerance
           %more no. of cstr stages, can tighten tolerance more
      
       herp=herp+1;
       
       u0=u0new;     
       %Concentration calculations
       c0 = C(a-1,:);
       conc = fsolve(@massbalance,c0);
       
       %Volume calculations
       Q0=Vnew;
       vol=fsolve(@volumechange,Q0);
       
       Vnew = vol;
       qnew = Vnew(1)+Vnew(2); %new qe qb
       u0new = qnew/area;
       
       %calc of new gas density rollgnew and visnew
       emassfracnew=(conc(1:11).*mw*vol(1))./sum(conc(1:11).*mw*vol(1));
       bmassfracnew=(conc(12:22).*mw*vol(2))./sum(conc(12:22).*mw*vol(2));
       rollgnew=sum(mw.*conc(1:11).*emassfracnew)*(1-efs(a))+sum(mw.*conc(12:22).*bmassfracnew)*efs(a);
       miunew=sum(vis.*emassfracnew)*(1-efs(a))+sum(vis.*bmassfracnew)*efs(a);
       %calc of ut
       gravitationterm=g*(rollp-rollgnew);
       uts(a)=((((1.78*10^-2)*gravitationterm^2)/(rollgnew*miunew))^0.3333)*...
       dp;
       %calc of umfnew and umbnew
       umf =(miunew/(dp*rollgnew))*((33.7^2+0.0408*(dp^3*rollgnew*((rollp-rollgnew)*g)/miunew^2))^(1/2)-33.7); % m/s
       umb=umf*((2300*rollgnew^0.126*miunew^0.523*exp(0.716*Ff))/(dp^0.8*g^0.934*(rollp-rollgnew)^0.934));%m/s
  
       
       db0=0.00376*(u0new *100-umb*100)^2; %cm
       dbm=0.652*(area*10000*(u0new *100-umb*100))^0.4; %cm
       dbs(a)=(dbm/100)-((dbm/100)-(db0/100))*exp((-0.3*h(a)/d)); %m
       ubs(a)=u0new -umb+0.711*(g*dbs(a))^(0.5);
       efs(a)=(u0new -umb)/(ubs(a)-umb); 
       kbe=(umb/3)+((4*D*emf*sum(ubs(1:a))/a)/(pi*sum(dbs(1:a))/a))^(1/2);

       K(a)=6*kbe/(dbs(a));
       
 
       
       %Modifying volume matrix using emulfraction
       Vnew(3)=efs(a)*(V(a-1,3)+V(a-1,4)); %new Ve
       Vnew(4)=(1-efs(a))*(V(a-1,3)+V(a-1,4)); %new Vb
       
       if herp==600
           break
       end
       
       end
       derp(a)=herp;
       u0s(a)=u0new;
       C(a,:)= conc; %Assigning a-1 stage solution to conc matrix
       V(a,:)= Vnew; %Assigning a-1 stage solution to vol matrix
       if herp==600
           break
       end
       
end
%assigning V2 and C2
V2=sum(V(n,1:2)); %combined qe and qb, 1 element
C2=(C(n,1:11)*V(n,1)+C(n,12:22)*V(n,2))/V2;
F2=C2.*V2;
%defining Cout and calc Vout
Fout=zeros(1,11); Fout(3:11)=F2(3:11);
Fr=zeros(1,11); Fr(1:2)=F2(1:2);
Vout=sum(Fout)*6250*T/P; Vrnew=sum(Fr)*6250*T/P; %ideal gas correlation causing problem ...
%(Vout>V2, as such 8314 is reduced to ensure V2=Vout+Vr)
Cout=Fout./Vout; Cr=Fr./Vrnew;%combined conc of emul and bubble, 11 elements

%if abs(Vr-Vrnew)<=1
    %fprintf('great success')
    %break
%else
    %Vr=Vrnew
%end
if herp2==100
    outfunc=zeros(1,15);
    outfunc(1:15)=404;
    fprintf('recycle failed')
    break
end

end

%% finding molar, mass flow and component fractions
   if herp==600
       outfunc=zeros(1,15);
    outfunc(1:15)=202;
        derp
        fprintf('Cant solve within 600 iterations')
        %parameters to call for troubleshoot
        u0
        u0new
        Vr
        V2
        C2
        Vout
        Cout
        Vrnew
        Cr
        return
   end 
   if herp2==100
       return
   end
   
   for ii=1:n
    emolarflow(ii,1:11)=C(ii,1:11).*V(ii,1);
    bmolarflow(ii,1:11)=C(ii,12:22).*V(ii,2);
   end

totalmolarflow = emolarflow + bmolarflow;

   for jj=1:11
    totalmassflow (:,jj)= mw(jj).*totalmolarflow(:,jj);
   end
   
   
 %out of reactor  
rxtorexitTmolarflow=sum(totalmolarflow(n,:));
rxtorexitTmassflow=sum(totalmassflow(n,:));

rxtorexitcompmolefrac = totalmolarflow(n,:)/sum(totalmolarflow(n,:));
rxtorexitcompmassfrac=totalmassflow(n,:)/sum(totalmassflow(n,:));

rxtorexitcompmolarflow=totalmolarflow(n,:);
rxtorexitcompmassflow=totalmassflow(n,:);

%after separation of methanol
Tmolarflowout=sum(Fout);
Tmassflowout=sum(Fout.*mw);

compmolarfracout=Fout/sum(Fout);
compmassfracout=(Fout.*mw)/sum(Fout.*mw);

compmolarflowout=Fout;
compmassflowout=Fout.*mw;
%% conversion calculation
meohconversion=(molarflowmeth-totalmolarflow(n,1))*100/molarflowmeth; % percent, conversion after reactor

c2c3massflow=sum(totalmassflow(n,4:5)); %kg/s
plantcapacity=c2c3massflow*31.536; %kt/a, assuming 365 days

inletmethanolmassflow=molarflowmeth*mw(1)*31.536; %kt/a

%% matrix for storing in outside function
outfunc=zeros(1,18);
outfunc(1:11)=rxtorexitcompmolefrac;%totalmassflow(n,:).*31.536;
outfunc(12)=rxtorexitTmassflow.*31.536;
outfunc(13)=plantcapacity;
outfunc(14)=d;
outfunc(15)=catalystweight;
outfunc(16)=inletmethanolmassflow;
outfunc(17)=WHSV;
outfunc(18)=meohconversion;



%% misc parameters for calling
molarflow0;
d;
catalystweight;
rollgnew;
miunew;
umf;
umb;
V2;
C2;
Vout;
Cout;
Vrnew;
Cr;
end
