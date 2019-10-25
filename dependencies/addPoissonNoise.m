function [b_sim_noisy] = addPoissonNoise(alpha, bkg, b_sim)

lambda      = zeros(size(b_sim));
lambda(:)   = b_sim*alpha + bkg; 
b_sim_noisy = poissrnd(lambda);
end


function rnd = local_poissrnd (lambda, varargin)
if (nargin < 1)
    print_usage ();
end

  if (nargin == 1)
    sz = size (lambda);
  elseif (nargin == 2)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error ("poissrnd: dimension vector must be row vector of non-negative integers");
    end
  elseif (nargin > 2)
    if (any (cellfun (@(x) (~ isscalar (x) || x < 0), varargin)))
      error ("poissrnd: dimensions must be non-negative integers");
    end
    sz = [varargin{:}];
  end

  if (~isscalar (lambda) && ~isequal (size (lambda), sz))
    error ("poissrnd: LAMBDA must be scalar or of size SZ");
  end

  if (~isreal (lambda))
    error ("poissrnd: LAMBDA must not be complex");
  end

  if (isa (lambda, "single"))
    cls = "single";
  else
    cls = "double";
  end

  if (isscalar (lambda))
    if (lambda >= 0 && lambda < Inf)
      rnd = randp (lambda, sz);
      if (strcmp (cls, "single"))
        rnd = single (rnd);
      end
    else
      rnd = NaN (sz, cls);
    end
  else
    rnd = NaN (sz, cls);

    k = (lambda >= 0) & (lambda < Inf);
    rnd(k) = randp (lambda(k));
  end

end