%Dynamic matrix control with input saturation constrained
%CBA 2020 

%% 
clear, close, clc
R=4;                                          %Prediction horizon 
L=2;                                          %Control horizon 
alpha =0.7;                                   %Parameter determines how fast
                                              %the trajectory reaches the setpoint
T=2e-3;                                       %sampling time
t=0:T:10e-2;                                  %time vector
N = length(t);                                %Model horizon
l=0.29;                                       %Beam length (m)
g=-9.8;                                       %Gravity (m/s^2)
m = 0.00388;                                  %mass
d = 0.065;                                    %Motor arm
r = 0.015;                                    %Ball radius
I = ((2*m*r^2)/3);                            %moment of inertia
s = tf('s');
G = -m*g*d/l/(I/r^2+m)/s^2                    %Transfer function TF
umax=(1.5e+06)                                %Saturation
Tz = c2d(G, T)                                %Discretized TF
Ad = cell2mat(Tz.den);                        %A(z-¹)
Bd = cell2mat(Tz.num);                        %B(z-¹)
[Ac Bc Cc Dc] = tf2ss([1.318],[1 0 0])        %Discrete space-state
[A1,B1]=c2d(Ac,Bc,T);
C1 = Cc; D1 = Dc;

%% Anti-Windup
n=size(B1,1); m=size(B1,2);
ops = sdpsettings('solver','sedumi','sedumi.eps',1e-5);
Qa=sdpvar(n,n,'symmetric');
La =  sdpvar(m,n, 'full');
Xa=sdpvar(n,n, 'full');
Ua = sdpvar(m,m,  'symmetric');
xp = sdpvar(n,1);
mu=sdpvar(1);
x_aw=[0;0];
LMI=[[-Qa -La' zeros(n,m) (C1*Qa+D1*La)' (A1*Qa+B1*La)';
    -La -2*Ua eye(m) (D1*Ua)' (B1*Ua)';
    zeros(m,n) eye(m) -mu*eye(m) zeros(m,n) zeros(m,m);
    (C1*Qa+D1*La) (D1*Ua) zeros(m,m) -eye(m) zeros(m,n);
    (A1*Qa+B1*La) (B1*Ua) zeros(n,m) zeros(n,m) -Qa]<=0];
AW = optimizer(LMI,mu,ops,xp,{Qa,La,mu});
sol = AW{x_aw};
Fa = sol{2}*inv(sol{1});
fprintf('\n\n'); disp('AW-Ganho AW: '); fprintf('%f \t', Fa);
%%
% Open loop response
x_aw=[0 0]'; u=1; y=0;
for k=1:N
    x_aw(:,k+1)=A1*x_aw(:,k)+B1*u;
    y(k)=C1*x_aw(:,k)+D1*u;
end
%% DMC-AW
a=y(2:end);                                     %Step response coefficients
h(1)=a(1);
for j=2:length(a);
    h(j) = a(j) - a(j-1);
end
f=0;                                           %Motion suppression parameter
A = toeplitz(a(1:R), [a(1) zeros(1,L-1)]);
disp('Matriz dinâmica A'); disp(A);
Kc = inv(A'*A + f*eye(L))*A';
KT = Kc(1,:);
disp('Kc'); disp(Kc);
y_aw(1)=0; x_aw=[0 0]';
u_aw=0; r(1:20)=20; r(21:41)=15;
% Coprime Factors
a_=A1+B1*Fa;
b_=B1;
c_=C1+D1*Fa;
d_=D1;
Mz_Ia=a_;Mz_Ib=b_; Mz_Ic=Fa; Mz_Id=0;
Nza=a_;Nzb=b_; Nzc=c_; Nzd=d_;
ud=0; x1=[0;0]; yd=0; u_d=0;y_d=0; y_in=0;
for k=1:1+40
    time(k) = (k-1)*T;
    for m=1:R
        S(m) =0;
        for i=m+1:N-2
            if k+m-i>0
                S(m)=S(m) +h(i)*deltau(k+m-i);
            end
        end
    end
    for i=1:R
        P(i)=0;
        for m=1:i
            P(i)=P(i)+S(m);
        end
    end
    E(k) = r(k)-y_in;
    for i=1:R
        El(i)= (1-alpha^i)*E(k)-P(i);
    end
    deltau(k) = KT*El';
    if k==1
        u_aw(k) = deltau(k);
    else
        u_aw(k) = u_aw(k-1) +deltau(k)-ud;
    end
    u_ant=u_aw(k);
    if u_aw(k)<-umax, u_aw(k)=-umax;end
    if u_aw(k)>umax, u_aw(k)=umax;end
    u_dep=u_aw(k);
    %  Anti windup Application
    u_til(k)=u_ant-u_dep;
    % M(z)-I
    x1(:,k+1)=(A1+B1*Fa)*x1(:,k)+B1*u_til(k);
    ud=Fa*x1(:,k);
    % Nz(z) = G2M(z)
    yd=(C1+D1*Fa)*x1(:,k)+D1*u_til(k);
    %%
    x_aw(:,k+1)=A1*x_aw(:,k)+B1*u_aw(k);
    y_aw(k+1)=C1*x_aw(:,k)+D1*u_aw(k);
    %
   % y_aw=awgn(y_aw, 60,'measured'); 
    y_in=y_aw(k+1)+yd; 
end

%% Without AW 
x=[0 0]'; u=1; y=0;
for k=1:N
    x(:,k+1)=A1*x(:,k)+B1*u;
    y(k)=C1*x(:,k)+D1*u;
end
%DMC
a=y(2:end);                                     %Step response coefficients
h(1)=a(1);
for j=2:length(a);
    h(j) = a(j) - a(j-1);
end
f=0;
A = toeplitz(a(1:R), [a(1) zeros(1,L-1)]);
disp('Matriz dinâmica A'); disp(A);
Kc = inv(A'*A + f*eye(L))*A';
KT = Kc(1,:);
disp('Matriz Kc '); disp(Kc);
yr(1)=0; x=[1 5]';
u=0; r(1:20)=20; r(21:41)=15;
for k=1:1+40
    time(k) = (k-1)*T;
    for m=1:R
        S(m) =0;
        for i=m+1:N-2
            if k+m-i>0
                S(m)=S(m) +h(i)*deltau(k+m-i);
            end
        end
    end
    for i=1:R
        P(i)=0;
        for m=1:i
            P(i)=P(i)+S(m);
        end
    end
    E(k) = r(k)-yr(k);
    for i=1:R
        El(i)= (1-alpha^i)*E(k)-P(i);
    end
    deltau(k) = KT*El';
    if k==1
        u(k) = deltau(k);
    else
        u(k) = u(k-1) +deltau(k);
    end
    x(:,k+1)=A1*x(:,k)+B1*u(k);
    yr(k+1)=C1*x(:,k)+D1*u(k);   
end

%% figures 
figure(1); 
plot(y_aw(2:end), 'k', 'linewidth', 2.5); hold on; plot(yr(1:end-1), 'g', 'linewidth', 2.5); hold on; stairs(r,'b--','linewidth', 2); hold on; grid on; legend('DMC with AW', 'DMC', 'Reference')
set(gca,'fontsize',20,'fontname','times new roman');  xlim([0 40])
ylabel('Output (cm)'); xlabel('Time (s)'); 
figure (2); stairs(u_aw, 'k','linewidth',2); hold on; stairs(u, 'g--','linewidth',2);  grid on;
set(gca,'fontsize',20,'fontname','times new roman');   ylabel('Control signal'); xlabel('Time (s)'); xlim([0 40])
