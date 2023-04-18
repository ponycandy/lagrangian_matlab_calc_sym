%% 双节cartpole倒立摆例子
clc;clear all;
syms m_cart m1 m2 l1 l2 g F
syms x theta1 theta2 Dx Dtheta1 Dtheta2
%% 计算各质量的速度
%cart
T=1/2*m_cart*Dtheta1*Dtheta1;
%第一质量块
xiangduisudu_sudu=Dtheta1*l1*[cos(theta1),sin(theta1)];
liandaisudu=[Dx,0];
jueduisudu=xiangduisudu_sudu+liandaisudu;
jueduisudu_square=jueduisudu(1)^2+jueduisudu(2)^2;
T_pole_1=1/2*m1*jueduisudu_square;
%第二质量块
xiangduisudu_sudu=Dtheta2*l2*[cos(theta2),sin(theta2)];
jueduisudu=xiangduisudu_sudu+jueduisudu;
jueduisudu_square=jueduisudu(1)^2+jueduisudu(2)^2;
T_pole_2=1/2*m2*jueduisudu_square;
%% 势能计算
P_pole_1=m1*g*l1*(1-cos(theta1));
P_pole_2=-m2*g*(l2*cos(theta2)-l1*(1-cos(theta1)));

%% 拉格朗日函数计算

Larg=T+T_pole_1+T_pole_2-P_pole_1-P_pole_2;

q  = [x, theta1,theta2];
Dq = [Dx, Dtheta1,Dtheta2];
Eq = LagrangeDynamicEqDeriver(Larg, q, Dq);
%% 计算NMPC所需矩阵并生成文件

tau=[F;0;0];


DDq=solve_for_DDq(Eq,tau,q);
new_func=[Dx,DDq.DDx,Dtheta1,DDq.DDtheta1,Dtheta2,DDq.DDtheta2];
state=[x,Dx,theta1,Dtheta1,theta2,Dtheta2];
d_state=new_func;
control=[F];
jac_mat=solve_for_DDQ_jacobian(new_func,control,state,q);

matlabFunction(d_state,'File','state_file');
matlabFunction(jac_mat,'File','jac_mat_file');



