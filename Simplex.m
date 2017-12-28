function [xopt, optval] = Simplex(c, A, b, AddSlack = true, NonNegative = true, Minimize=true, LowerOrEqual = true)
  Settings = struct("AddSlack", AddSlack, "Negative", ~NonNegative, "Minimize", Minimize, "LowerOrEqual", LowerOrEqual);
  
  if ~Settings.Minimize
    c = -c;
  end
  
  if ~Settings.LowerOrEqual
    A = -A;
    b = -b;
  end
  
  if Settings.Negative
    A = [A -A];
    c = [c -c];
  end
  
  if Settings.AddSlack
    A = [A eye(size(A,1))];
    c = [c zeros([1 size(A,1)])];
    B = 1:size(A,2);
    B = B(size(A,2)-size(A,1)+1 : end);
  else
    [B, s, A, b] = findBase(A,b);
    if ~s
      optval = NaN;
      xopt = Nan*ones([size(A,2) 1]);
      printf("Problem is Unsolvable");
    end
  end
  
  [xopt optval] = SimplexAlgorithm(c, A, b, B);
  
  if optval == -inf
    printf("Problem is Indefinite");
  end
  if optval == NaN
    printf("Problem is Unsolvable");
  end
  
  if ~Settings.Minimize
    optval = -optval;
  end
  
  if Settings.Negative
    h = xopt(1:end-size(A,1));
    h = h(1:end/2) - h(end/2 +1 : end);
    xopt = [h; xopt(end-size(A,1)+1:end)];
  end
  
end

function [B, solvable, M, r] = findBase(A, b)
  solvable = true;
  I = find(b < 0);
  A(I,:) = -A(I,:);
  b(I) = -b(I);
  
  c = [zeros([1 size(A,2)]) ones([1 size(A, 1)])];
  A = [A eye(size(A,1))];
  B = 1 : size(A,2);
  B = B(size(A,2)-size(A,1)+1 : end);
  
  [xopt, optval, B, A, r] = SimplexAlgorithm(c, A, b, B);
  if abs(optval) >= 0.001
    solvable = false;
  end
  M = A(:,1:size(A,2)-size(A,1));
end

function [xopt, optval, Base, M, r] = SimplexAlgorithm(c, A, b, B)
  optval = 0;
  xopt = zeros([size(A,2) 1]);
  
  for i = 1:size(B,2)
    optval = optval + (-c(B(i))/A(i,B(i)))*b(i);
    c = c + (-c(B(i))/A(i,B(i)))*A(i,:);
  end
  
  M = A;
  r=b;
  if size(A,1) ~= size(b,1)
    optval = NaN;
    xopt = NaN*ones([size(A,2) 1]);
    Base = NaN;
    return;
  elseif size(A,2) ~= size(c,2)
    optval = NaN;
    xopt = NaN*ones([size(A,2) 1]);
    Base = NaN;
    return;
  end
  
  while size(find(c<0), 2) > 0
    I = find(c<0);
    i = I(1);
    
    a = A(:,i);
    N = find(a > 0);
    if size(N,1) == 0
      optval = -inf;
      xopt(B) = b;
      xopt(i) = inf;
      Base = B;
      M = A;
      r = b;
      return;
    end
    T = b(N)./a(N);
    n = N(find(T==min(T))(1));
    
    B(n) = i;
    
    a = A(n, :);
    x = b(n);
    I = 1 : size(A,1);
    I = setdiff(I, n);
    
    for k = I
      b(k) = b(k) + (-A(k,i)/a(i))*x;
      A(k,:) = A(k,:) + (-A(k,i)/a(i))*a;
    end
    
    optval = optval + (-c(i)/a(i))*x;
    c = c + (-c(i)/a(i))*a;
    
  end
  optval = -optval;
  xopt(B) = b;
  Base = B;
  M = A;
  r = b;
end