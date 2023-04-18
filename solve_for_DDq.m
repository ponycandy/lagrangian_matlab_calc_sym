function output=solve_for_DDq(leftside,right_side,q)

N = length(leftside);

DDq = sym(zeros(1, N));
for ii = 1:N
    DDq(ii) = sym(['DD', char(q(ii))]);
end

Dq = sym(zeros(1, N));
for ii = 1:N
    Dq(ii) = sym(['D', char(q(ii))]);
end
eqns= leftside==right_side;

output=solve(eqns,DDq,"real",true);

