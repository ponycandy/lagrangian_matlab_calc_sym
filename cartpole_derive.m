clc;clear all;

%% cartpole
syms x theta  Dx Dtheta  %务必使用D而不是d
syms l m1 m2 t g F

%% Kinetic and Potential Energy

T_cart=1/2*m1*Dx^2;

xiangduisudu_sudu=Dtheta*l*[cos(theta),sin(theta)];
liandaisudu=[Dx,0];
jueduisudu=xiangduisudu_sudu+liandaisudu;
jueduisudu_square=jueduisudu(1)^2+jueduisudu(2)^2;
T_pole=1/2*m2*jueduisudu_square;

P_pole=m2*g*l*(1-cos(theta));

Larg=T_cart+T_pole-P_pole;
%%
q  = [x, theta];
Dq = [Dx, Dtheta];
Eq = LagrangeDynamicEqDeriver(Larg, q, Dq);

tau=[F;0];


DDq=solve_for_DDq(Eq,tau,q);
new_func=[Dtheta,DDq.DDtheta,Dx,DDq.DDx];
state=[theta,Dtheta,x,Dx];
d_state=new_func;
control=[F];
jac_mat=solve_for_DDQ_jacobian(new_func,control,state,q);

matlabFunction(d_state,'File','state_file');
matlabFunction(jac_mat,'File','jac_mat_file');