clear;

%% Predefined Parameters
% the mass of the carrige (kg)
M = 1;      

% the effective pendulum length (m)
L = 0.842;    

% the friction coefficient (kg/s)
F = 1;      

% the gravitational acceleration (m/s)
g = 9.8093;             

B1_2 = 0; B3 = 0; B4_5 = 0; B6_7 = 0; B8 = 0; B9 = 0; B10 = 0;

%% State Space Representations
% system matrix
A = [0     1      0   0;
     0   -F/M     0   0;
     0     0      0   1;
     0   F/(L*M) g/L  0];           
 
% input matrix 
B = [    0;
        1/M;
         0;
      -1/(L*M)];                    

% output matrix  
C = [1 0 0 0;
     0 0 1 0];                      

% feedthrough matrix 
D = [0;
     0];                            
 
%% Reachability (Controllability)
% reachability matrix
R = ctrb(A,B); 

rank_R = rank(R);

%% Observability
% observability matrix
O = obsv(A,C);  

rank_O = rank(O);

%% B1 & B2
% state-space model
states = {'s(t)' 's_dot' 'phi(t)' 'phi_dot'}; inputs = {'u(t)'}; outputs = {'s(t)'; 'phi(t)'};
sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);                  

% transfer function model
sys_tf = tf(sys_ss);                
sys_zpk = zpk(sys_tf);            % pole-zero format

poles = pole(sys_tf);
isstable(sys_tf);                 % verify stability

% desired closed loop eigenvalues
P = [-0.5+1i -0.5-1i -8 -9];

% solve for K using Ackerman Formula instead of Place to remove mismatching
K = acker(sys_ss.A,-sys_ss.B,P);
% disp('Feedback Gain Matrix: '); disp(K);

% create closed loop system
Acl = A+(B*K);  
sys_cl = ss(Acl,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

% four various initial states
x0_1 = [-0.5 0 0 0];
x0_2 = [0 -0.5 0 0];
x0_3 = [0 0 -0.7 0];
x0_4 = [0 0 0 -0.5];

t = 0:.1:12; equ = zeros(1,121);

% simulate time response of dynamic system
u = 0*t;    % set reference equal to zero 
[yl1,t,xl1] = lsim(sys_cl,u,t,x0_1); [yl2,t,xl2] = lsim(sys_cl,u,t,x0_2);
[yl3,t,xl3] = lsim(sys_cl,u,t,x0_3); [yl4,t,xl4] = lsim(sys_cl,u,t,x0_4);
ul1 = K*xl1'; ul2 = K*xl2'; ul3 = K*xl3'; ul4 = K*xl4';

% initial response of the closed loop system
if (B1_2 == 1)
    analysis1(xl1,xl2,xl3,xl4,t,equ,1);
    
    figure(2)
    subplot(2,2,1); plot(t,ul1); hold on; plot(t,equ,'k--','LineWidth',0.05);
    xlabel('Time (seconds)'); ylabel('u(t)');
    title('Reponse of u(t) to Initial Condition: [-0.5 0 0 0]');
    
    subplot(2,2,2); plot(t,ul2); hold on; plot(t,equ,'k--','LineWidth',0.05);
    xlabel('Time (seconds)'); ylabel('u(t)');
    title('Reponse of u(t) to Initial Condition: [0 -0.5 0 0]');
    
    subplot(2,2,3); plot(t,ul3); hold on; plot(t,equ,'k--','LineWidth',0.05);
    xlabel('Time (seconds)'); ylabel('u(t)');
    title('Reponse of u(t) to Initial Condition: [0 0 -0.7 0]');
    
    subplot(2,2,4); plot(t,ul4); hold on; plot(t,equ,'k--','LineWidth',0.05);
    xlabel('Time (seconds)'); ylabel('u(t)');
    title('Reponse of u(t) to Initial Condition: [0 0 0 -0.5]');    
end

%% B3
% simulate time response of dynamic system
[t,x1] = ode45(@(t,x) plant(t,x,K,M,L,F,g),t,x0_1'); [t,x2] = ode45(@(t,x) plant(t,x,K,M,L,F,g),t,x0_2');
[t,x3] = ode45(@(t,x) plant(t,x,K,M,L,F,g),t,x0_3'); [t,x4] = ode45(@(t,x) plant(t,x,K,M,L,F,g),t,x0_4');

% initial response of the closed loop system without linearlization
if B3 == 1
    analysis1(x1,x2,x3,x4,t,equ,3);
    
    x3_infos = lsiminfo(x3(:,1),t,0); x3_infophi = lsiminfo(x3(:,3),t,0);
    xl3_infos = lsiminfo(xl3(:,1),t,0); xl3_infophi = lsiminfo(xl3(:,3),t,0);
    
    figure(4)
    subplot(2,1,1); yyaxis left; plot(t,xl3(:,1)); hold on; plot(t,equ,'k--','LineWidth',0.05);
    title('Linearized system: reponse of y(t) to Initial Condition: [0 0 -0.7 0]'); ylabel('y_1: s(t) (m)');
    p=find(xl3(:,1)==max(xl3(:,1))); 
    text(t(p),xl3(p,1),['(',num2str(t(p)),',',num2str(xl3(p,1)),')'],'color','b');
    p=find(xl3(:,1)==min(xl3(:,1))); 
    text(t(p),xl3(p,1),['(',num2str(t(p)),',',num2str(xl3(p,1)),')'],'color','b');
    text(xl3_infos.SettlingTime,0,(num2str(xl3_infos.SettlingTime)),'color','b');
    
    yyaxis right; plot(t,xl3(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; %plot(t,equ,'k--','LineWidth',0.05);
    p=find(xl3(:,3)==max(xl3(:,3))); 
    text(t(p),xl3(p,3),['(',num2str(t(p)),',',num2str(xl3(p,3)),')'],'color','r');
    p=find(xl3(:,3)==min(xl3(:,3))); 
    text(t(p),xl3(p,3),['(',num2str(t(p)),',',num2str(xl3(p,3)),')'],'color','r');
    text(xl3_infophi.SettlingTime,0,(num2str(xl3_infophi.SettlingTime)),'color','r');
    
    subplot(2,1,2); yyaxis left; plot(t,x3(:,1)); hold on; plot(t,equ,'k--','LineWidth',0.05);
    title('Nonlinear system: reponse of y(t) to Initial Condition: [0 0 -0.7 0]'); ylabel('y_1: s(t) (m)');
    p=find(x3(:,1)==max(x3(:,1))); 
    text(t(p),x3(p,1),['(',num2str(t(p)),',',num2str(x3(p,1)),')'],'color','b');
    p=find(x3(:,1)==min(x3(:,1))); 
    text(t(p),x3(p,1),['(',num2str(t(p)),',',num2str(x3(p,1)),')'],'color','b');
    text(x3_infos.SettlingTime,0,(num2str(x3_infos.SettlingTime)),'color','b');
    
    yyaxis right; plot(t,x3(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; plot(t,equ,'k--','LineWidth',0.05);
    p=find(x3(:,3)==max(x3(:,3))); 
    text(t(p),x3(p,3),['(',num2str(t(p)),',',num2str(x3(p,3)),')'],'color','r');
    p=find(x3(:,3)==min(x3(:,3))); 
    text(t(p),x3(p,3),['(',num2str(t(p)),',',num2str(x3(p,3)),')'],'color','r');
    text(x3_infophi.SettlingTime,0,(num2str(x3_infophi.SettlingTime)),'color','r');    
end

%% B4 & B5
% sampling time and period
% completely unstable
% T = 0.15; period = 3;
% T = 0.14; period = 4.48;
% T = 0.136; period = 5.44;
% T = 0.135; period = 6.75;
% T = 0.134; period = 6.7;
% T = 0.133; period = 7.98; 

% oscilation
% T = 0.132; period = 13.2*4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = 0.13; period = 13*4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = 0.129; period = 12.9*3;
% T = 0.127; period = 19.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = 0.126; period = 12.6*4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = 0.125; period = 12.5;
% T = 0.123; period = 12.3;
% T = 0.12; period = 12;
% T = 0.119; period = 11.9; 

% asymptoyiclly stable 
T = 0.08; period = 12;


t = 0:T:period; times = round(T\period); equ = zeros(1,times+1); 

s1 = zeros(times+1,4); s1(1,:) = x0_1;
s2 = zeros(times+1,4); s2(1,:) = x0_2;
s3 = zeros(times+1,4); s3(1,:) = x0_3;
s4 = zeros(times+1,4); s4(1,:) = x0_4;

% simulate time response of dynamic system
for i = 1:times
    t_ = 0:0.001:T;
    zoh1 = K*x0_1'; zoh2 = K*x0_2'; zoh3 = K*x0_3'; zoh4 = K*x0_4';
    [t1_,xs1] = ode45(@(t,x) plant_s(t,x,zoh1,M,L,F,g),t_,x0_1');
    [t2_,xs2] = ode45(@(t,x) plant_s(t,x,zoh2,M,L,F,g),t_,x0_2');
    [t3_,xs3] = ode45(@(t,x) plant_s(t,x,zoh3,M,L,F,g),t_,x0_3');
    [t4_,xs4] = ode45(@(t,x) plant_s(t,x,zoh4,M,L,F,g),t_,x0_4');
    x0_1 = xs1(end,:); x0_2 = xs2(end,:); x0_3 = xs3(end,:); x0_4 = xs4(end,:);
    s1(i+1,:) = x0_1; s2(i+1,:) = x0_2; s3(i+1,:) = x0_3; s4(i+1,:) = x0_4;
end

% initial response of the closed loop system without linearlization
if B4_5 == 1 
    analysis3(s1,s2,s3,s4,t,equ,5);
end

%% B6 
sys_ssd = c2d(sys_ss,T);

%% B7
x0_1 = [-0.5 0 0 0];
x0_2 = [0 -0.5 0 0];
x0_3 = [0 0 -0.7 0];
x0_4 = [0 0 0 -0.5];

% desired closed loop eigenvalues
Pd = [exp(P(1)*T) exp(P(2)*T) exp(P(3)*T) exp(P(4)*T)];

% solve for K using Ackerman Formula instead of Place to remove mismatching
Kd = acker(sys_ssd.A,-sys_ssd.B,Pd);
% disp('Feedback Gain Matrix: '); disp(K);

% create closed loop system
Acld = sys_ssd.A+(sys_ssd.B*Kd);  
sys_cld = ss(Acld,sys_ssd.B,sys_ssd.C,sys_ssd.D,T,'statename',states,'inputname',inputs,'outputname',outputs);

% simulate time response of dynamic system
u = 0*t;    % set reference equal to zero
[yl1d,t,xl1d] = lsim(sys_cld,u,t,x0_1); [yl2d,t,xl2d] = lsim(sys_cld,u,t,x0_2);
[yl3d,t,xl3d] = lsim(sys_cld,u,t,x0_3); [yl4d,t,xl4d] = lsim(sys_cld,u,t,x0_4);

% initial response of the closed loop system
if (B6_7 == 1)
    analysis3(xl1d,xl2d,xl3d,xl4d,t,equ,6);
end

%% B8 
% sampling time and period
% completely unstable
% Td = 0.23; periodd = Td*30;
% Td = 0.209; periodd = Td*50;
% Td = 0.208; periodd = Td*200;

% oscilation
% Td = 0.207; periodd = Td*400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Td = 0.2; periodd = 20*2;
% Td = 0.198; periodd = 19.8*3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Td = 0.197; periodd = 19.7;
% Td = 0.194; periodd = 19.4;
% Td = 0.189; periodd = 18.9;
% Td = 0.186; periodd = 18.6;
% Td = 0.184; periodd = 18.4;
% Td = 0.183; periodd = 18.3;

% asymptoyiclly stable 
% Td = 0.182; periodd = 18.2;
% Td = 0.18; periodd = 12.6;
Td = 0.08; periodd = 12;

td = 0:Td:periodd; timesd = round(Td\periodd); equd = zeros(1,timesd+1); 

s1d = zeros(timesd+1,4); s1d(1,:) = x0_1;
s2d = zeros(timesd+1,4); s2d(1,:) = x0_2;
s3d = zeros(timesd+1,4); s3d(1,:) = x0_3;
s4d = zeros(timesd+1,4); s4d(1,:) = x0_4;

% simulate time response of dynamic system
for i = 1:timesd
    t_ = [0 Td];
    zoh1d = Kd*x0_1'; zoh2d = Kd*x0_2'; zoh3d = Kd*x0_3'; zoh4d = Kd*x0_4';
    [t1_,xs1d] = ode45(@(t,x) plant_s(t,x,zoh1d,M,L,F,g),t_,x0_1');
    [t2_,xs2d] = ode45(@(t,x) plant_s(t,x,zoh2d,M,L,F,g),t_,x0_2');
    [t3_,xs3d] = ode45(@(t,x) plant_s(t,x,zoh3d,M,L,F,g),t_,x0_3');
    [t4_,xs4d] = ode45(@(t,x) plant_s(t,x,zoh4d,M,L,F,g),t_,x0_4');
    x0_1 = xs1d(end,:); x0_2 = xs2d(end,:); x0_3 = xs3d(end,:); x0_4 = xs4d(end,:);
    s1d(i+1,:) = x0_1; s2d(i+1,:) = x0_2; s3d(i+1,:) = x0_3; s4d(i+1,:) = x0_4;
end

% initial response of the closed loop system without linearlization
if B8 == 1
    analysis3(s1d,s2d,s3d,s4d,td,equd,7);
end

%% B9
x0_1 = [-0.5 0 0 0];
x0_2 = [0 -0.5 0 0];
x0_3 = [0 0 -0.7 0];
x0_4 = [0 0 0 -0.5];

% control law 
Q = [10   0   0   0;     % penalize s(t) error
      0  0.1  0   0;     % penalize s_dot(t) error
      0   0  10   0;     % penalize phi(t) error
      0   0   0   0.1];  % penalize phi_dot(t) error

R = 1;                   % penalize actuator effort

% solve for K using LQR 
[Kr,Sr,er] = dlqr(sys_ssd.A,-sys_ssd.B,Q,R);

% calculate the optimal cost
J1 = x0_1*Sr*x0_1'; J2 = x0_2*Sr*x0_2'; J3 = x0_3*Sr*x0_3'; J4 = x0_4*Sr*x0_4';
J = [J1 J2 J3 J4];

% create closed loop system
Aclr = sys_ssd.A+(sys_ssd.B*Kr);  
sys_clr = ss(Aclr,sys_ssd.B,sys_ssd.C,sys_ssd.D,T,'statename',states,'inputname',inputs,'outputname',outputs);

% simulate time response of dynamic system
u = 0*t;    % set reference equal to zero
[yl1r,t,xl1r] = lsim(sys_clr,u,t,x0_1); [yl2r,t,xl2r] = lsim(sys_clr,u,t,x0_2);
[yl3r,t,xl3r] = lsim(sys_clr,u,t,x0_3); [yl4r,t,xl4r] = lsim(sys_clr,u,t,x0_4);

% initial response of the closed loop system
if (B9 == 1)
    analysis2(xl1r,xl2r,xl3r,xl4r,t,equ,8)
end

%% B10
% sampling time and period
% completely unstable
% Tr = 1; periodr = Tr*5;
% Tr = 0.36; periodr = Tr*30;
% Tr = 0.347; periodr = Tr*120;

% oscilation
% Tr = 0.346; periodr = Tr*200;
% Tr = 0.34; periodr = Tr*200;
% Tr = 0.339; periodr = Tr*200;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tr = 0.334; periodr = Tr*200;
% Tr = 0.326; periodr = Tr*200;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tr = 0.325; periodr = Tr*200;
% Tr = 0.324; periodr = Tr*100;
% Tr = 0.3; periodr = Tr*50;
% Tr = 0.297; periodr = Tr*80;

% asymptoyiclly stable 
% Tr = 0.296; periodr = Tr*80;
% Tr = 0.18; periodr = 12.6;
Tr = 0.08; periodr = 12;

tr = 0:Tr:periodr; timesr = round(Tr\periodr); equr = zeros(1,timesr+1); 

s1r = zeros(timesr+1,4); s1r(1,:) = x0_1;
s2r = zeros(timesr+1,4); s2r(1,:) = x0_2;
s3r = zeros(timesr+1,4); s3r(1,:) = x0_3;
s4r = zeros(timesr+1,4); s4r(1,:) = x0_4;

% simulate time response of dynamic system
for i = 1:timesr
    t_ = [0 Tr];
    zoh1r = Kr*x0_1'; zoh2r = Kr*x0_2'; zoh3r = Kr*x0_3'; zoh4r = Kr*x0_4';
    [t1_,xs1r] = ode45(@(t,x) plant_s(t,x,zoh1r,M,L,F,g),t_,x0_1');
    [t2_,xs2r] = ode45(@(t,x) plant_s(t,x,zoh2r,M,L,F,g),t_,x0_2');
    [t3_,xs3r] = ode45(@(t,x) plant_s(t,x,zoh3r,M,L,F,g),t_,x0_3');
    [t4_,xs4r] = ode45(@(t,x) plant_s(t,x,zoh4r,M,L,F,g),t_,x0_4');
    x0_1 = xs1r(end,:); x0_2 = xs2r(end,:); x0_3 = xs3r(end,:); x0_4 = xs4r(end,:);
    s1r(i+1,:) = x0_1; s2r(i+1,:) = x0_2; s3r(i+1,:) = x0_3; s4r(i+1,:) = x0_4;
end

% initial response of the closed loop system without linearlization
if B10 == 1  
    analysis2(s1r,s2r,s3r,s4r,tr,equr,9);
end