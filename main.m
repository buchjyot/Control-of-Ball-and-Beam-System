%% Clear Workspace
% MiniProject - Modern Control Theory
clear;close all; clc;

%% Part-3: Linearized System Dynamics
L = 1;
J = 0.0676;
m = 0.9048;
r = 0.03;
Jb = 0.000326;
g = 9.81;

% Initial conditions
p0 = 0.25;
x1 = p0;
v0 = 0;
u0 = m*g*p0;
x0 = [x1;v0;0;0];

A = [0 1 0 0;
    0 0 -m*g*r^2/(Jb + m*r^2) 0;
    0 0 0 1;
    -m*g/(m*x1^2+J+Jb) 0 0 -2*m*x1*v0/(m*x1^2+J+Jb)];
B = [0;0;0;1/(m*x1^2+J+Jb)];
C = [1 0 0 0];
D = 0;

% Eigenvalues
[V,Di]=jordan(A);
fprintf('Eigenvalues: \n\n');
disp(diag(Di));

% Impulse Response
tfinal = 10;
sysC = ss(A,B,C,D);
[~,t,x] = impulse(sysC,tfinal);
x(:,1) = x(:,1) + p0*ones(size(x(:,1)));
figure(1);

ylabelArray = {'$p(t)$','$\dot{p}(t)$','$\theta(t)$','$\dot{\theta}(t)$'};
titleArray = {'Ball Position','Ball Velocity','Beam Angle','Beam Angular Velocity'};
for i = 1:4
    sp = subplot(3,2,i);
    plt1 = plot(t,x(:,i),'LineStyle','-','Color','b');
    set(gca,'XLim',[0 0.6]);
    ylabel(ylabelArray{i});
    title(titleArray{i});
end

% Legend Plot
sh=subplot(3,2,5);
p=get(sh,'position');
lh=legend(sh,plt1,'Open Loop');
set(lh,'position',p);
axis(sh,'off');

%% Part-4: System Analysis
sysObsv = obsv(A,C);
sysCtrb = ctrb(A,B);
fprintf('Rank of Observability Matrix: %d \n\n',rank(sysObsv));
fprintf('Rank of Controllability Matrix: %d \n\n',rank(sysCtrb));
fprintf('Since the system is both controllable and observable, it is minimal.\n\n');

%% Part-5: State Feedback Controller Design
P = [-5 -4 -3 -2];
K = acker(A,B,P);
Aaug = A-B*K;
sysAug = ss(Aaug,B,C,D);
[~,t,x] = impulse(sysAug,tfinal);
x(:,1) = x(:,1) + p0*ones(size(x(:,1)));
figure(1);

% Verify Eigenvalues
fprintf('Eigenvalues of augmented system with full state feedback:\n\n');
disp(eig(Aaug));

%% Part-6: State Feedback Controller Design Plot
for i = 1:4
    subplot(3,2,i)
    hold on;
    plt2 = plot(t,x(:,i));
    set(gca,'XLim',[0 0.6]);
    ylabel(ylabelArray{i});
    title(titleArray{i});
end
hold off;

% Actuation effort
subplot(3,2,6);
u = -K*x';
tau_plt1 = plot(t,u,'LineStyle','-.');
title('Applied Torque');ylabel('$\tau(t)$');

% Legend Plot
sh=subplot(3,2,5);
p=get(sh,'position');
lh=legend(sh,[plt1 plt2 tau_plt1],'Open Loop','Closed Loop','$\tau(t)$ Full State Feedback');
set(lh,'position',p);
axis(sh,'off');

%% Part-7: Observer based Feedback Control Design

% Observer poles 10 times the controller poles
O = 10*P;
G = place(A',C',O)';
Aaug = [...
    A -B*K;
    G*C A-G*C-B*K
    ];
Baug = [B;B];
Caug = [C zeros(size(C))];
x0_observer = 0.5*ones(4,1);

% Verify Eigenvalues
fprintf('### Eigenvalues of augmented system controller dynamics with observer based feedback:\n\n');
disp(eig(A-B*K));
fprintf('### Eigenvalues of augmented system controller dynamics with observer based feedback:\n\n');
disp(eig(A-G*C));

sysAug = ss(Aaug,Baug,Caug,D);
x0Aug = [0;0;0;0;x0_observer];
[~,timp,ximp] = impulse(sysAug,tfinal);
[~,tinit,xinit] = initial(sysAug,x0Aug,timp);
xadd = ximp + xinit;
xadd(:,1) = xadd(:,1) + p0*ones(size(xadd(:,1)));
xadd(:,5) = xadd(:,5) + p0*ones(size(xadd(:,1)));
figure(1);
for i = 1:4
    subplot(3,2,i)
    hold on;
    plt3 = plot(tinit,xadd(:,i),'Color','g');
    plt4 = plot(tinit,xadd(:,4+i),'--');
    set(gca,'XLim',[0 0.6]);
    ylabel(ylabelArray{i});
    title(titleArray{i});
end
hold off;

% Actuation effort
subplot(3,2,6);
hold on;
u = -K*xadd(:,5:8)';
tau_plt2 = plot(t,u,'LineStyle','-.');

% Legend Plot
sh=subplot(3,2,5);
p=get(sh,'position');
lh=legend(sh,[plt1 plt2 plt3 plt4 tau_plt1 tau_plt2],...
    '$x(t)$:Open Loop States',...
    '$x(t)$:Full State Feedback','$x(t)$:Observer based Feedback',...
    '$\hat{x}(t)$:Observer based Feedback','$\tau(t)$:Full State Feedback',...
    '$\tau(t)$:Observer based Feedback');
set(lh,'position',p);
axis(sh,'off');

%% Part-8: Optimal State Feedback Controller Design
Q = diag([100 10 10 10]);
R = 10;
[K,S,e] = lqr(sysC,Q,R);

Qest = 0.01*diag([1 1 1 1]);
Rest = 0.01;
G = eye(4);
[LGain,P,E] = lqe(A,G,C,Qest,Rest);

Aaug = [...
    A -B*K;
    LGain*C A-LGain*C-B*K
    ];
Baug = [B;B];
Caug = [C zeros(size(C))];
sys_aug = ss(Aaug,Baug,Caug,D);
x0_observer = 0.5*ones(4,1);

% Verify Eigenvalues
fprintf('*** Eigenvalues of augmented system controller dynamics with LQG:\n\n');
disp(eig(A-B*K));
fprintf('*** Eigenvalues of augmented system observer dynamics with LQG:\n\n');
disp(eig(A'-G'*C'));

sysAug = ss(Aaug,Baug,Caug,D);
x0Aug = [0;0;0;0;x0_observer];
[yimp,timp,ximp] = impulse(sysAug,tfinal);
[yinit,tinit,xinit] = initial(sysAug,x0Aug,timp);
xadd = ximp + xinit;
xadd(:,1) = xadd(:,1) + p0*ones(size(xadd(:,1)));
xadd(:,5) = xadd(:,5) + p0*ones(size(xadd(:,1)));
figure(1);
for i = 1:4
    subplot(3,2,i)
    hold on;
    plt5 = plot(tinit,xadd(:,i),'Color','m');
    plt6 = plot(tinit,xadd(:,4+i),'--');
    set(gca,'XLim',[0 1]);
    ylabel(ylabelArray{i});
    title(titleArray{i}); 
end
hold off;

% Actuation effort
subplot(3,2,6);
hold on;
u = -K*xadd(:,5:8)';
tau_plt3 = plot(tinit,u,'LineStyle','-.');
set(gca,'XLim',[0 0.5]);

% Legend Plot
sh=subplot(3,2,5);
p=get(sh,'position');
lh=legend(sh,[plt1 plt2 plt3 plt4 plt5 plt6 tau_plt1 tau_plt2 tau_plt3],...
    '$x(t)$:Open Loop States',...
    '$x(t)$:Full State Feedback','$x(t)$:Observer based Feedback',...
    '$\hat{x}(t)$:Observer based Feedback','$x(t)$:LQG States',...
    '$\hat{x}(t)$:LQG Estimates','$\tau(t)$:Full State Feedback',...
    '$\tau(t)$:Observer based Feedback','$\tau(t)$:LQG Input');
set(lh,'position',p);
axis(sh,'off');
print('-fillpage','Report/fig1','-dpdf')