%Mirror Information
%Coordinates origin as drawing 13783 rev D, pg 4/8
%Note: Drawing orient X axis accoring to the back face of the mirror
%Therefore X dimensions shall be sign reversed. For everything
%Z axis is also reversed
%T matrix is introduced to rotate the drawing coordinates to the 
%parent figure coordinate system and from inches to meters
T=0.0254*[-1 0 0;0 1 0;0 0 -1];

% units in meters

%Crude Approximation base on drawings
%CG coord Center of gravity coordinates, respect to the glassCoord
CG_Mirror=[0.02; -0.026;0.295];
%Hardpoints Puck positions respect to the glassCoord
W=zeros(3,6);
W(:,1)=[-62.001;-104.938;0];
W(:,2)=[62.001;-104.938;0];
W(:,3)=[121.879;-1.225;0];
W(:,4)=[59.878;106.163;0];
W(:,5)=[-59.878;106.163;0];
W(:,6)=[-121.879;-1.225;0];

W=T*W;

%Vertex positions respect to the cell considering the origin of the cell
%Coordinates system with a 0,0,0 at the same position that the mirror
%Coordinate system in its nominal position
V=zeros(3,6);
V(:,1)=[0;-175.797;99.327];
V(:,2)=[0;-175.797;99.327];
V(:,3)=[152.245;87.899;99.327];
V(:,4)=[152.245;87.899;99.327];
V(:,5)=[-152.245;87.899;99.327];
V(:,6)=[-152.245;87.899;99.327];

V=T*V;
%Hp vectors
Wcg=zeros(3,6);
%Nominal length of HPs
Lhp=zeros(1,6);
for i1=1:6 
   Wcg(:,i1)=W(:,i1)-CG_Mirror; 
   Vcg(:,i1)=V(:,i1)-CG_Mirror; 
end
Hp=W-V;
Hpn=Hp;
for i1=1:6
    Hpn(:,i1)=Hp(:,i1)/norm(Hp(:,i1));
    Lhp(i1)=norm(Hp(:,i1));
end


%Nominal Stiffness of the Hardpoint
%a vector Accounts for individual variations
%units in N/m
HpK=[170e6; 170e6; 170e6; 170e6; 170e6; 170e6];

%Stiffness Matrix
dF_dx=zeros(3,1);
dF_dy=zeros(3,1);
dF_dz=zeros(3,1);
dF_drx=zeros(3,1);
dF_dry=zeros(3,1);
dF_drz=zeros(3,1);
dM_dx=zeros(3,1);
dM_dy=zeros(3,1);
dM_dz=zeros(3,1);
dM_drx=zeros(3,1);
dM_dry=zeros(3,1);
dM_drz=zeros(3,1);

for i1=1:6
    dF_dx=dF_dx+Hpn(1,i1)*(-1)*HpK(i1)*Hpn(:,i1);
    dF_dy=dF_dy+Hpn(2,i1)*(-1)*HpK(i1)*Hpn(:,i1);
    dF_dz=dF_dz+Hpn(3,i1)*(-1)*HpK(i1)*Hpn(:,i1);
    dF_drx=dF_drx+((Wcg(3,i1)*Vcg(2,i1)-Wcg(2,i1)*Vcg(3,i1))/norm(Hp(:,i1)))*(-1)*HpK(i1)*Hpn(:,i1);
    dF_dry=dF_dry+((Wcg(1,i1)*Vcg(3,i1)-Wcg(3,i1)*Vcg(1,i1))/norm(Hp(:,i1)))*(-1)*HpK(i1)*Hpn(:,i1);
    dF_drz=dF_drz+((Wcg(2,i1)*Vcg(1,i1)-Wcg(1,i1)*Vcg(2,i1))/norm(Hp(:,i1)))*(-1)*HpK(i1)*Hpn(:,i1);
    dM_dx=dM_dx+(Hp(1,i1)/norm(Hp(:,i1)))*(-1)*HpK(i1)*cross(Wcg(:,i1),Hpn(:,i1));
    dM_dy=dM_dy+(Hp(2,i1)/norm(Hp(:,i1)))*(-1)*HpK(i1)*cross(Wcg(:,i1),Hpn(:,i1));
    dM_dz=dM_dz+(Hp(3,i1)/norm(Hp(:,i1)))*(-1)*HpK(i1)*cross(Wcg(:,i1),Hpn(:,i1));
    dM_drx=dM_drx+((Wcg(3,i1)*Vcg(2,i1)-Wcg(2,i1)*Vcg(3,i1))/norm(Hp(:,i1)))*(-1)*HpK(i1)*cross(Wcg(:,i1),Hpn(:,i1));
    dM_dry=dM_dry+((Wcg(1,i1)*Vcg(3,i1)-Wcg(3,i1)*Vcg(1,i1))/norm(Hp(:,i1)))*(-1)*HpK(i1)*cross(Wcg(:,i1),Hpn(:,i1));
    dM_drz=dM_drz+((Wcg(2,i1)*Vcg(1,i1)-Wcg(1,i1)*Vcg(2,i1))/norm(Hp(:,i1)))*(-1)*HpK(i1)*cross(Wcg(:,i1),Hpn(:,i1));
end
%Stiffness Matrix K finally;
Km=[dF_dx dF_dy dF_dz dF_drx dF_dry dF_drz; dM_dx dM_dy dM_dz dM_drx dM_dry dM_drz];
clear('dF_dx','dF_dy','dF_dz','dF_drx','dF_dry','dF_drz');
clear('dM_dx','dM_dy','dM_dz','dM_drx','dM_dry','dM_drz');
%Inertial matrix
Mmass=17000; % Kg
Ix=17000*(4.2^2)/4;
Iy=17000*(4.2^2)/4;
Iz=17000*(4.2^2)/2;
Mm=diag([Mmass Mmass Mmass Ix Iy Iz]);
clear('Mmass','Ix','Iy','Iz');
A=inv(Mm)*Km;
[Mvect Mfrec]=eig(A);
GainM=inv(Mvect)*inv(Mm);
Mdamp=0.1*[1 1 1 1 1 1];

%Frequency of sampled model of the mirror dynamic response and Harpoint
%response
Fsmm=500; %Hz
Mode1MirrorTf=c2d(tf([1],[1 2*Mdamp(1)*abs(sqrt(Mfrec(1,1))) abs(Mfrec(1,1))]),1/Fsmm,'zoh');
Mode2MirrorTf=c2d(tf([1],[1 2*Mdamp(2)*abs(sqrt(Mfrec(2,2))) abs(Mfrec(2,2))]),1/Fsmm,'zoh');
Mode3MirrorTf=c2d(tf([1],[1 2*Mdamp(3)*abs(sqrt(Mfrec(3,3))) abs(Mfrec(3,3))]),1/Fsmm,'zoh');
Mode4MirrorTf=c2d(tf([1],[1 2*Mdamp(4)*abs(sqrt(Mfrec(4,4))) abs(Mfrec(4,4))]),1/Fsmm,'zoh');
Mode5MirrorTf=c2d(tf([1],[1 2*Mdamp(5)*abs(sqrt(Mfrec(5,5))) abs(Mfrec(5,5))]),1/Fsmm,'zoh');
Mode6MirrorTf=c2d(tf([1],[1 2*Mdamp(6)*abs(sqrt(Mfrec(6,6))) abs(Mfrec(6,6))]),1/Fsmm,'zoh');

clear('Mfrec','A','Mdamp');
%Inverse Kinematix linear matrix.
%This matrix is used to offset the forces on the mirror
%dynamics block
dL_dx=zeros(6,1);
dL_dy=zeros(6,1);
dL_dz=zeros(6,1);
dL_drx=zeros(6,1);
dL_dry=zeros(6,1);
dL_drz=zeros(6,1);
for i1=1:6
    dL_dx(i1)=Hpn(1,i1);
    dL_dy(i1)=Hpn(2,i1);
    dL_dz(i1)=Hpn(3,i1);
    dL_drx(i1)=(Wcg(3,i1)*Vcg(2,i1)-Wcg(2,i1)*Vcg(3,i1))/norm(Hp(:,i1));
    dL_dry(i1)=(Wcg(1,i1)*Vcg(3,i1)-Wcg(3,i1)*Vcg(1,i1))/norm(Hp(:,i1));
    dL_drz(i1)=(Wcg(2,i1)*Vcg(1,i1)-Wcg(1,i1)*Vcg(2,i1))/norm(Hp(:,i1));
end
KinematicsM=[dL_dx dL_dy dL_dz dL_drx dL_dry dL_drz];
clear('dL_dx','dL_dy','dL_dz','dL_drx','dL_dry','dL_drz');
InvKinematicsM=inv(KinematicsM);

%inverse Kinematics *stifness matrix
KmXInvK=Km*InvKinematicsM;

%CG motion to Surface motion Transformation matrix
% This matrix transforms coordinates from Mirror CG motion to Optical
% surface origin point as defined in drawing 14558 page 2/6 and drawing 13784 page
% 2/2
%Surface origin point coordinates,  units in meters
So=[0;0;0.4783];
Socg=So-CG_Mirror;
CG2So=[1 0 0 0 Socg(3) -Socg(2); 0 1 0 -Socg(3) 0 Socg(1);0 0 1 Socg(2) -Socg(1) 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
SoKinematics=KinematicsM*inv(CG2So);
%MHP Matrix as UA-95 document
MHP=zeros(6);
MHP(1:3,:)=Hpn;
for i1=1:6
    MHP(4:6,i1)=cross(Wcg(:,i1),Hpn(:,i1));
end


%Hard point position system
%Electric Motor Inductance per phase
HDLph=0.0014; %H

%Electric Motor Resistance per phase
HDRph=0.33; %Ohm

%Electric Motor Torque Constant
HDKm=66; %Nm/A

%Electric Motor Moment of Inertia 
HDIm=8.9e-4; %Kg m^2

%Electric motor Back EMF constant
HDKv=7.4/3*161/(2*pi*60); %Vs

%Motor shaft to Lineal speed
HPKw2Lin=1e-3/(161*(2*pi));%m/rads

%Torque Controller
Wn=2*pi*10;
damp=1;
HDKpTor=Wn/(HDKm*HDRph/HDLph);
HDKiTor=Wn*(1-HDLph/HDRph)*HDKpTor;


%Speed control
%Plant 1/Is
%PI controller Kpv+Kiv/s
%H model
%(wn1+wn2)*(s+wn1*wn2/(wn1+wn2)
%------------------------------
%(s+wn1)(s+wn2)

wn1=10*2*pi;
wn2=1*2*pi;
Kpv=(wn1+wn2)*HDIm;
Kiv=(wn1*wn2)*HDIm;
Pv=tf([1],[HDIm 0]);
Cv=tf([Kpv Kiv],[1 0]);
Hv=minreal(Pv*Cv/(1+Pv*Cv));

%position Controller
wm=2.5*2*pi;
dm=0.9;
w0=5*wm;
Kpp=(wm^2+2*dm*wm*w0)/(wn1*HPKw2Lin);
Kip=Kpp*w0*wm/(wm+2*dm*w0);
Pp=tf([HPKw2Lin],[1 0]);
Cp=tf([Kpp Kip],[1 0]);
Hp=minreal(Pp*Cp*Hv/(1+Pp*Cp*Hv));
if 0
    bode(Hv,Hp)
    grid
    figure(2)
    step(Hv,Hp)
    grid on
end
HPdtf=c2d(Hp,1/Fsmm,'zoh');
clear('HDLph','HDRph','HDKm','HDIm','HDKv','HPKw2Lin','HDKpTor','HDKiTor');
clear('Wn','damp','wn1','wn2','Kpv','Kiv','Pv','Cv','Hv','wm','dm','w0');
clear('Kpp','Kip','Pp','Cp','Hp')


%Global Force Controllers
% Made up by 6 independent discrete PID+n controllers
% updated at Fsfc Hz
Fscf=100;, % Hz
Cfkp=0.0*ones(1,6);
Cfki=0.2*ones(1,6);
Cfkd=0*ones(1,6);
Cfn=1*ones(1,6);

