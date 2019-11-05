function op = opPCISampling(n,idx)

if (min(idx) < 1) || (max(idx) > n)
   error('Index parameter must be integer and match dimensions of the operator');
end

op = @(x,mode) op_internal(n,idx,x,mode);

function y = op_internal(n,idx,x,mode)

if mode == 1
   y = x(idx);
elseif mode == 2 
   y = zeros(n,1);
   y(idx) = x;
end