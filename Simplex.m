function [optval, xopt] = Simplex(c, A, b)
  xopt = FindBaseSolution(A,b);
  optval = c*xopt;
  s = FindBaseSolution(A, b);
  if(s.find(NaN))
    return
  end
  [optval, xopt] = SimpexHelp(c, A, b, s);
  
end

function [start] = FindBaseSolution(A, b)
  c = zeros([1, size(A,2)]);
  c = [c ones([1, size(A,1)])];
  start = NaN*ones([size(A,2), 1]);
  A = [A eye(size(A,1))];
  [ov ox] = SimpexHelp(c, A, b, [zeros(size(start)), ones(size(A,1))]);
  if(abs(ov) < 0.001)
    start = ox(1:size(start,2));
  end
end

function [optval, xopt] = SimpexHelp(c, A, b, B)
  B = logical(B);
  if(size(A(:,B)) != size(eye(size(A,1))))
    
  end
end