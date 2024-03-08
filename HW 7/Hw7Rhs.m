function out = Hw7Rhs(t,y)
if t==0
    A = [[0 1];[-1 0]];
else
    A = [[0 1];[-1 -3/t]];
end

out = A*y;
end