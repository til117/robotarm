% Project/Lab3 - Computer Aided Control Design
% Students: Ioannis Til

%% Transfer function [Assignment 1]

% Independent cell
clear all; close all; clc;

% Parameters
Lm = 2;
Rm = 21;
b = 1;
Kt = 38;
Km = 0.5;
n= 1/20;
[J,umax] = lab3robot(930901); % J=5, umax=90

% Transfer function of robot arm
s = tf('s');
G = (n*Kt)/(((J*s+b)*(s*Lm+Rm)+Km*Kt)*s)

% Verification
lab3robot(G,930901)

%% P-Controller [Assignment 2]

% Independent cell
clear all; close all; clc;

% Transfer function of robot arm
s = tf('s');
G = tf([0 0 0 1.9],[10 107 40 0]);

% Find K that gives <5% overshoot and short rise time on the
% unit step response of the closed loop system
figure('Name','Different K values','NumberTitle','off')
hold on
K = linspace(4,4.2,5); % Test with values around 4 (gives %5)
legendindex = {zeros(1,length(K))}; 
for i=1:length(K);
    Gc = feedback(K(i)*G,1); % The closed loop system
    step(Gc,50); % The step response of Gc
    legendindex{i} = ['K = ' num2str(K(i))];
end
legend(legendindex)
% From the above loop we note a value just below 4.13

% P-Controller
figure('Name','Unit step response of the closed loop system','NumberTitle','off')
Kp = 4.12; % Choose K just under 4.13
Gc = feedback(Kp*G,1); % The closed loop system
step(Gc) % The step response of Gc
legend('Kp = 4.12')
Info=stepinfo(Gc)

%% System Parameters [Assignment 3]

% Independent cell
clear all; close all; clc;

% Transfer function of robot arm
s = tf('s');
G = tf([0 0 0 1.9],[10 107 40 0]);

% The open loop system
figure('Name','Open loop system','NumberTitle','off')
Kp = 4.12; % Unstable for K>225.2
Go = Kp*G;
margin(Go)
[Gm,Fm,Wp,Wc]=margin(Go)

% The closed loop system
figure('Name','Closed loop system','NumberTitle','off')
Gc = feedback(Kp*G,1);
margin(Gc)
BW = bandwidth(Gc)

%% Lead-Lag Compensator [Assignment 6]

% Independent cell
clear all; close all; clc;

% Transfer function of robot arm
s = tf('s');
G = tf([0 0 0 1.9],[10 107 40 0]);

% P-Controller
Kp = 4.12;
Go = Kp*G; % The open loop system
[Gm,Fm,Wp,Wc]=margin(Go);
Gc = feedback(Kp*G,1); % The closed loop system

% We study different parameters to fulfill the specifications
Wc2 = Wc*4 % We increase the rise time with a factor 4
[~,Fm2] = bode(Go,Wc2) % The new phase margin
Fm_new = Fm-(180+Fm2)+5.45 % Compensate for phase margin
beta = (1-sind(Fm_new))/(1+sind(Fm_new)) % Find beta

% Find the appropriate K
Wc2 = Wc*4;
tD = 1/(Wc2*sqrt(beta));
tmpFlead = (tD*s+1)/(beta*tD*s+1); % Neutralised Lead-Controller (K=1)
[Am,~]=bode(tmpFlead*G,Wc2);
K = 1/Am % The value that satisfies: K*|Flead|*|G|=1

% Lead-Controller
Flead = (tD*s+1)/(beta*tD*s+1)

% Find the appropriate gamma
%{
syms z
G = ((1/20)*38)/(((5*z+1)*(z*2+21)+0.5*38)*z);
e1=limit(1/(G*K*z),z,0); % e1=1.6282
%}
gamma = 0.05/1.6282 % gamma*e1=(Stat.control err.)=0.05
tI=14/Wc2; % 14 trial and error

% Lag-Controller
Flag = (tI*s+1)/(tI*s+gamma);

% Lead-Lag Compensator
F = K*Flead*Flag

% Comparing the open/closed loop systems
figure('Name','Open loop systems','NumberTitle','off')
Go2 = F*G;
bode(Go2,Go,'-.')
legend('Lead-Lag Controller','P-Controller')
figure('Name','Step response of closed loop systems','NumberTitle','off')
Gc2 = feedback(Go2,1);
step(Gc2,Gc,'-.')
legend('Lead-Lag Controller','P-Controller')

%{
%Checking if |u|<umax for all t
figure
u=F*(1-Gc2)
%step(u)
%
t=0:0.1:100;
S2=1/(1+Go2);
[Y,T]=lsim(S2,t,t);%checking stationary error smaller than 0.05
plot(T,Y), xlabel('t'), ylabel('error'), title('Error with ramp')
%}

%% Sensitivity Functions [Assignment 8]

% Dependent cell - The previous cell needs to run first to obtain the 
% parameters for this cell (G,Kp,F)
close all; clc;

% The Sensitivity Functions
S1 = 1/(1+G*Kp);
S2 = 1/(1+G*F);
bodemag(S1,'-.',S2), grid on
legend('P-Controller','Lead-Lag Controller')

%% Robustness [Assignment 9]

% Dependent cell - The previous cell needs to run first to obtain the 
% parameters for this cell (S2)
close all; clc;

% Robustness Criterion: |dGx| < 1/|T|
T=1-S2;
dG1 = (s+10)/40; 
dG2 = (s+10)/(4*(s+0.01));
figure('Name','T(s) & 1/dG','NumberTitle','off')
bodemag(T^-1,dG1,'--',dG2,':')
legend('1/T(s)','DG_1','DG_2')
% From the plot we note that dG2 does not fulfill the Robustness Criterion

%% State Space Model [Assignment 11]

% Independent cell
clear all; close all; clc;


% Parameters
Lm = 2;
Rm = 21;
b = 1;
Kt = 38;
Km = 0.5;
n= 1/20;
[J,umax] = lab3robot(930901); % J=5, umax=90

% Static Space Model
A=[0 n 0; 0 -b/J Kt/J; 0 -Km/Lm -Rm/Lm];
B=[0; 0; 1/Lm];
C=[1 0 0];

S=[B A*B A^2*B]; % Controllability Matrix
O=[C; C*A; C*A^2]; % Observability Matrix

% Check if the system is controllable
if det(S)==0  
    disp('The system is not controllable!')
else
    disp('The system is controllable!')
end

% Check if the system is observable
if det(O)==0
    disp('The system is not observable!')
else
    disp('The system is observable!')
end

%{
% Static space using MATLAB (wrong?)
s = tf('s');
G = ((1/20)*38)/(((5*s+1)*(s*2+21)+0.5*38)*s);
Static_Space = ss(G)
%}

%% Pole Placement [Assignment 12]

% Dependent cell - The previous cell needs to run first to obtain the 
% parameters for this cell (A,B,C)
close all; clc;

s = tf('s');
G = tf([0 0 0 1.9],[10 107 40 0]);
Gs = ss(G)
poles = eig(Gs.a)
p1=0; p2=-0.3879; p3=-10.3121;
% The task of deciding the correct pole-placement was assigned to a 
% group of well-trained dolphins. After some (not that much) time, the
% following results were obtained:       (p.184 - 9.2)
p1=-2+2i; p2=-2-2i; p3=-1.5;
L=place(A,B,[p1 p2 p3]) % Placing the poles
Gc0=ss(A-B*L,B,C,0)
L0=1/dcgain(Gc0)

%% Evaluation

% Independent cell
clear all; close all; clc;

% Everything needed for validation of the derived controllers/models:
s = tf('s');
% Assignment 1 (OK)
G = tf([0 0 0 1.9],[10 107 40 0]);
% Assignment 2 (OK)
Kp = 4.12;
% Assignment 6 (OK)
beta=0.1701; K=12.9351; tD=3.4081; tI=19.678; gamma=0.0307;
Flead = (tD*s+1)/(beta*tD*s+1)
Flag = (tI*s+1)/(tI*s+gamma);3
F = K*Flead*Flag
% Assignment 11 (OK)
A=[0 0.05 0; 0 -0.2 7.6; 0 -0.25 -10.5];
B=[0; 0; 0.5];
C=[1 0 0];
% Assignment 12 (OK)
p1=-2+2i; p2=-2-2i; p3=-1.5;
L=place(A,B,[p1 p2 p3]);
Gc0=ss(A-B*L,B,C,0);
L0=1/dcgain(Gc0);

% Validation
lab3robot(G,Kp,F,A,B,C,L,L0,930901)
