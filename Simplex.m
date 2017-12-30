function [xopt, optval] = Simplex(c, A, b, AddSlack = true, NonNegative = true, Minimize=true, LowerOrEqual = true)
  #usage: [xopt, optval] = Simplex(targetFunction, SystemMatrix, righthandSide, AddingSlackVariables, SolutionNonNegative, Minimization, LowerOrEqual)
  #
  #This is a function to Solve a Linear Problem via the Simplex Algorithm
  #It takes the target function (the c in c*x),
  #the System Matrix (the A from A*x)
  #and the right hand side of the Equation (the b from Ax <= b)
  #
  #Also the optional parameters:
  #@var{AddSlack}:
  #   A bool if the Algorithm needs to add Slackvariables or if they are already included
  #@var{NonNegative}:
  #   A bool if the solution should be all positive
  #@var{Minimize}:
  #   A bool if the Problem is a minimization problem
  #@var{LowerOrEqual}:
  #   A bool if the Problem is of the sort Ax<=b

  %Creating the Settings struct
  Settings = struct("AddSlack", AddSlack, "Negative", ~NonNegative, "Minimize", Minimize, "LowerOrEqual", LowerOrEqual);

  %if it is not a minimization make it to one  
  if ~Settings.Minimize
    c = -c;
  end
  
  %if the System is Ax>=b transform it to -Ax<=-b
  if ~Settings.LowerOrEqual
    A = -A;
    b = -b;
  end
  
  %if x can be negative create x=x_pos - x_neg with x_pos, x_neg >= 0
  if Settings.Negative
    A = [A -A];
    c = [c -c];
  end
  
  %if Slackvariables shall be added, add a slack variable for each equation
  if Settings.AddSlack
    A = [A eye(size(A,1))];
    c = [c zeros([1 size(A,1)])];
    %you get automaticaly a base (the slack variables)
    B = 1:size(A,2);
    B = B(size(A,2)-size(A,1)+1 : end);
  %if no Slackvariables need beeing added (Ax=b) then you need to find a base
  else
    %finding a Base also transform the right hand side and System Matrix
    [B, s, A, b] = findBase(A,b);
    %s is a bool if the original problem is solvable
    if ~s
      optval = NaN;
      xopt = Nan*ones([size(A,2) 1]);
      printf("Problem is Unsolvable");
      return;
    end
  end
  
  %Solve the Linear Problem of the kind min{cx : Ax=b, x>=0} with the Base B
  [xopt optval] = SimplexAlgorithm(c, A, b, B);
  
  %if the optimal value is -inf then the Problem is indefinite
  if optval == -inf
    printf("Problem is Indefinite");
  end
  
  if optval == NaN
    printf("Problem is unsolvable, check for Dimention Mismatch triggerd");
  end
  
  %if the problem was not originaly a Minimization, then the optimal value is actualy -optval
  if ~Settings.Minimize
    optval = -optval;
  end
  
  %if the original problem allowed negative solutions, put the solution together  
  if Settings.Negative
    h = xopt(1:end-size(A,1));
    h = h(1:end/2) - h(end/2 +1 : end);
    xopt = [h; xopt(end-size(A,1)+1:end)];
  end
  
end

function [B, solvable, M, r] = findBase(A, b)
  %initialize the solvability
  solvable = true;
  %find all entries of b where b is <0
  I = find(b < 0);
  %all equations where this is the case multiply by -1 on both sides
  A(I,:) = -A(I,:);
  b(I) = -b(I);
  
  %now: b>=0
  
  %Build the Help problem min{sum(d) : Ax + d = b, x,d>=0}
  %Per default d=b and x = 0 is a base solution of the problem
  c = [zeros([1 size(A,2)]) ones([1 size(A, 1)])];
  A = [A eye(size(A,1))];
  B = 1 : size(A,2);
  B = B(size(A,2)-size(A,1)+1 : end);
  
  %Run SimplexAlgorithm on this Problem  
  [xopt, optval, B, A, r] = SimplexAlgorithm(c, A, b, B);
  %if the optimal value is not near 0 (floating point arithmic), the original problem has no solution
  if abs(optval) >= 0.001
    solvable = false;
  end
  %remember the transformation to get the System Matrix in this form, no need to calculate it twice
  M = A(:,1:size(A,2)-size(A,1));
end

function [xopt, optval, Base, M, r] = SimplexAlgorithm(c, A, b, B)
  %optimal value initialize with 0 same for optimal x
  optval = 0;
  xopt = zeros([size(A,2) 1]);
  
  %bring c in the representation of the Base B
  for i = 1:size(B,2)
    optval = optval + (-c(B(i))/A(i,B(i)))*b(i);
    c = c + (-c(B(i))/A(i,B(i)))*A(i,:);
  end
  
  %remember the SystemMatrix and the right hand side, useful for
  M = A;
  r=b;

  %check for dimention mismatches  
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
  
  %as long as the target function can be minimized
  while size(find(c<0), 2) > 0
    %get indices of all x which can be increased so that the function gets smaller
    I = find(c<0);
    %choose the first one
    i = I(1);

    %get the collum of the system Matrix corresponding to that x value    
    a = A(:,i);
    %get indices of all entries which are positive
    N = find(a > 0);
    %if there are none, the problem is indefinite
    if size(N,1) == 0
      optval = -inf;
      xopt(B) = b;
      xopt(i) = inf;
      Base = B;
      M = A;
      r = b;
      return;
    end
    %divide the corresponding entries in b by the values in a
    T = b(N)./a(N);
    %find the minimum of these values and take the first if there are multiple
    n = N(find(T==min(T))(1));
    
    %The nth entry in B leaves the base and xi enters it
    B(n) = i;
  
    %we don't need the collum anymore we override it by the nth line  
    a = A(n, :);
    %also remember the nth entry in b
    x = b(n);
    %get all indices except n
    I = 1 : size(A,1);
    I = setdiff(I, n);
    
    %Transform the System Matrix and the right hand side so that xi is in the Base
    for k = I
      b(k) = b(k) + (-A(k,i)/a(i))*x;
      A(k,:) = A(k,:) + (-A(k,i)/a(i))*a;
    end
    
    %Transform c and evaluate new optimal value    
    optval = optval + (-c(i)/a(i))*x;
    c = c + (-c(i)/a(i))*a;
    
  end
  %the real optval is -optval
  optval = -optval;
  %set the xopt to equal b after all these transformations
  xopt(B) = b;
  %the Base is B
  Base = B;
  %remember the System Matrix and right hand side
  M = A;
  r = b;
end