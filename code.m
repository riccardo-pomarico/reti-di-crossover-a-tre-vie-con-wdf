%% Sound Analysis, Synthesis and Processing 
% Homework #4

%% Group components:
% Brusca Alfredo 	10936149
% Pomarico Riccardo 10661306

clear all
close all
clc

%% Import Input Audio Signal
[Vin,~] = audioread('ExpSweep.wav');

%% LTSpice Files for Ground-Truth
[OutLowSpice,~]=audioread('outlowsweep.wav');
[OutMidSpice,~]=audioread('outmidsweep.wav');
[OutHighSpice,FsLTSpice]=audioread('outhighsweep.wav');
TsLTSpice=1/FsLTSpice;

%% Sampling frequency (to be varied: FsLTSpice/downSampFact, with downSampFact={4,3,2})
downSampFact=2;
Fs =FsLTSpice/downSampFact; 

%% Downsample Input Signal
Vin=Vin([1:downSampFact:end]);

%% Sampling Period
Ts=1/Fs;
%% Number of Samples
Nsamp=length(Vin);
%% Simulated time
tstop=Nsamp*Ts;
%% Parameters of Dynamic Element
L1=0.35*10^(-3);
L2=0.35*10^(-3);
L3=3.5*10^(-3);
L4=3.5*10^(-3);
C1= 2.8*10^(-6);
C2= 2.8*10^(-6);
C3= 28*10^(-6);
C4= 4.7*10^(-6);
C5=28*10^(-6);
C6=47*10^(-6);
%% Resistive Parameters

R1=10;
RspkLow=8;
R2=10;
RspkMid=8;
RspkHigh=8;

%% WDF setting of free parameters (adaptation conditions)

Z3 = Ts/(2*C1);
Z5 = RspkHigh;
Z6 = 2*L1/Ts;
Z4 = Z5*Z6/(Z5+Z6);
Z2 = Z4;
Z1 = Z2+Z3;

Z9 = 2*L4/Ts;
Z11 = Ts/(2*C5);
Z15 = RspkLow;
Z17 = R2;
Z18 = Ts/(2*C6);
Z16 = Z17+Z18;
Z14 = Z16;
Z13 = Z14*Z15/(Z14+Z15);
Z12 = Z13;
Z10 = Z11*Z12/(Z11+Z12);
Z8 = Z10;
Z7 = Z8+Z9;

Z19 = Ts/(2*C2);
Z20 = 2*L3/Ts;
Z21 = RspkMid;
Z23 = Ts/(2*C3);
Z24 = R1;
Z25 = Ts/(2*C4);
Z26 = 2*L2/Ts;


%% Computation of Scattering Matrices

Bser = [1,1,1];
Qpar = [1,1,1];

% Diagonal impedance matrices
Zser1 = diag([Z1, Z2, Z3]);
Zpar1 = diag([Z4, Z5, Z6]);
Zser2 = diag([Z7, Z8, Z9]);
Zser3 = diag([Z16, Z17, Z18]);
Zpar2 = diag([Z10, Z11, Z12]);
Zpar3 = diag([Z13, Z14, Z15]);

% Calculate S matrices
sS1 = eye(3) - 2*Zser1*Bser'*(Bser*Zser1*Bser'\Bser);
sS2 = eye(3) - 2*Zser2*Bser'*(Bser*Zser2*Bser'\Bser);
sS3 = eye(3) - 2*Zser3*Bser'*(Bser*Zser3*Bser'\Bser);
sP1 = 2*Qpar'*(Qpar*(Zpar1\Qpar')\(Qpar/Zpar1)) - eye(3);
sP2 = 2*Qpar'*(Qpar*(Zpar2\Qpar')\(Qpar/Zpar2)) - eye(3);
sP3 = 2*Qpar'*(Qpar*(Zpar3\Qpar')\(Qpar/Zpar3)) - eye(3);


syms Z22 real
Z = diag([Z19, Z20, Z21, Z22, Z23, Z24, Z25, Z26]);
B = [1, 0, 0, 0; 
    0, 1, 0, 0; 
    0, 0, 1, 0;
    0, 0, 0, 1;
    -1, 0, 0, 1;
    -1, -1, -1, 1;
    -1, -1, -1, 1;
    0, 0, 0, 1]';
S = eye(8) - 2*Z*B'*((B*Z*B')\B);
Z22 = double(solve(S(4,4)==0));
Z = diag([Z19, Z20, Z21, Z22, Z23, Z24, Z25, Z26]);
S = eye(8) - 2*Z*B'*((B*Z*B')\B);

%% Initialization of Waves

a3 = 0;
b5 = 0;
a6 = 0;

a9 = 0;
a11 = 0;
b15 = 0;
b17 = 0;
a18 = 0;

a19 = 0;
b21 = 0;
a20 = 0;
a23 = 0;
b24 = 0;
a25 = 0;
a26 = 0;


%% Initialize Output Signals
% Low
VoutLow=zeros(size(Vin));
% Mid
VoutMid=zeros(size(Vin));
% High
VoutHigh=zeros(size(Vin));

ii=0;
while (ii<Nsamp)
    ii=ii+1;

    %% Manage Dynamic Elements
    b3 = a3;
    b6 = -a6;

    b9 = -a9;
    b11 = a11;
    b18 = a18;
    
    b19 = a19;
    b20 = -a20;
    b23 = a23;
    b25 = a25;
    b26 = -a26;        
    
    %% Forward Scan
    a4 = sP1(1, :) * [0; b5; b6];
    b2 = a4;
    a1 = sS1(1, :) * [0; b2; b3];

    a16 = sS3(1, :) * [0; b17; b18];
    b14 = a16;
    a13 = sP3(1, :) * [0; b14; b15];
    b12 = a13;
    a10 = sP2(1, :) * [0; b11; b12];
    b8 = a10;
    a7 = sS2(1, :) * [0; b8; b9];

    a22 = S(4, :) * [b19; b20; b21; 0; b23; b24; b25; b26];

    %% Local Root Scattering

    b1 = 2*Vin(ii) - a1;
    b7 = 2*Vin(ii) - a7;
    b22 = 2*Vin(ii) - a22;

    %% Backward Scan
    a2 = sS1(2, :) * [b1; b2; b3];
    a3 = sS1(3, :) * [b1; b2; b3];
    b4 = a2;
    a5 = sP1(2, :) * [b4; b5; b6];
    a6 = sP1(3, :) * [b4; b5; b6];

    a8 = sS2(2, :) * [b7; b8; b9];
    a9 = sS2(3, :) * [b7; b8; b9];
    b10 = a8;
    a11 = sP2(2, :) * [b10; b11; b12];
    a12 = sP2(3, :) * [b10; b11; b12];
    b13 = a12;
    a14 = sP3(2, :) * [b13; b14; b15];
    a15 = sP3(3, :) * [b13; b14; b15];
    b16 = a14;
    a18 = sS3(3, :) * [b16; b17; b18];  

    a19 = S(1, :) * [b19; b20; b21; b22; b23; b24; b25; b26];
    a20 = S(2, :) * [b19; b20; b21; b22; b23; b24; b25; b26];
    a21 = S(3, :) * [b19; b20; b21; b22; b23; b24; b25; b26];
    a23 = S(5, :) * [b19; b20; b21; b22; b23; b24; b25; b26];
    a25 = S(7, :) * [b19; b20; b21; b22; b23; b24; b25; b26];
    a26 = S(8, :) * [b19; b20; b21; b22; b23; b24; b25; b26];

    %% Read Output
    VoutLow(ii)=-(a15+b15)/2;
    VoutMid(ii)=-(a21+b21)/2;
    VoutHigh(ii)=-(a5+b5)/2;
    
end


%% Output Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(TsLTSpice*[1:length(OutLowSpice)],OutLowSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutLow,'b--','Linewidth',1); grid on; xlim([0,tstop]); 
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
title('Output Signals','Fontsize',18,'interpreter','latex');
subplot(312)
plot(TsLTSpice*[1:length(OutMidSpice)],OutMidSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutMid,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(TsLTSpice*[1:length(OutHighSpice)],OutHighSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutHigh,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

%% Error Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(Ts*[1:Nsamp],OutLowSpice([1:downSampFact:end])-VoutLow,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
title(['Error Signals. $F_{\mathrm{s}}=$ ',num2str(Fs),' Hz'],'Fontsize',18,'interpreter','latex');
subplot(312)
plot(Ts*[1:Nsamp],OutMidSpice([1:downSampFact:end])-VoutMid,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(Ts*[1:Nsamp],OutHighSpice([1:downSampFact:end])-VoutHigh,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

