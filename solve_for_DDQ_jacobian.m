function jac_mat=solve_for_DDQ_jacobian(new_func,control,state,q)



N = length(q);

DDq = sym(zeros(1, N));
for ii = 1:N
    DDq(ii) = sym(['DD', char(q(ii))]);
end

Dq = sym(zeros(1, N));
for ii = 1:N
    Dq(ii) = sym(['D', char(q(ii))]);
end

up_head=[control,state];
left_head=new_func;

jac_mat=jacobian(left_head,up_head);
