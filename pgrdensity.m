% File    : pgrdensity.m
% System  : MATLAB
% Purpose : Triangular PGR fabric plots.
% Author  : Frederick W. Vollmer
% Date    : Aug 6, 2020
% Notice  : Copyright (c) 2020 Frederick W. Vollmer 
% License : See LICENSE
%
% Citation
% --------
% The algorithms used in this code are presented in:
%
% Vollmer, F.W., 2020. Representing progressive fabric paths on a 
%   triangular plot using a fabric density index and crystal axes 
%   eigenvector barycenters. Geological Society of America Abstracts with
%   Proggrams, v. 51, submitted. 
% Vollmer F.W., 1990. An application of eigenvalue methods to structural 
%   domain analysis. Geological Society of America Bulletin, v. 102, n. 6,
%   786?791, ISSN 00167606
% Vollmer F.W., 1989. A triangular fabric plot with applications for 
%   structural analysis. EOS Transactions American Geophysical Union 
%   70:463
%
% One or more should be cited for usage of this or derivative code.
%-------------------------------------------------------------------------

function [x,y,z] = pgrdensity(n, err, opts, expect)
% PGRDENSITY  Returns a grid of density for a triangular fabric plot.
%   Returns [x,y,z] vectors with lengths of n. [x, y] values are in the 
%   range [-1,-1] to [1,1], the [x,y] coordinates are converted to [P,G,R] 
%   indexes, and then eigenvalues to calculate the fabric density index, 
%   D (Vollmer, 2020), in [z]. [x,y] points outside of the triangle return 
%   z = NaN. 
%   n = Number of grid nodes in x and y between -1 an 1.
%   err = Error allowed in PGR indexes at plot margins, this is used to 
%     avoid edge effects when contouring. Trial and error may be required. 
%     For n=150, err=0.005 works well.
%   opts 
%     0 = Return density index.
%     1 = Return intensity index.
%     2 = Return expected density index, expect must be protolith 
%         eigenvalues.
%     3 = Return expected intensity index, expect must be protolith 
%         eigenvalues.
%
  x = zeros(n,n);
  y = zeros(n,n);
  z = zeros(n,n);
  inc = 2.0/(n-1);
  xi = -1.0;
  for i = 1:n
    yi = -1.0;
    for j = 1:n
      zi = densityxy(xi, yi, err, opts, expect); 
      x(i,j) = xi;
      y(i,j) = yi;
      z(i,j) = zi;
      yi = yi + inc;
    end
    xi = xi + inc;
  end
end

function [d] = densityxy(x, y, err, opts, expect)
% DENSITYXY  Returns the fabric density index for [x,y] on unit triangle.
%   NaN is returned for points off the triangle more than +/- err.
%   [x,y] = Cartesian coordinates in the range [-1,-1] to [1,1] that 
%     include the triangle within a unit circumcircle.
%   err = Error allowed in PGR indexes at plot margins, this is used to 
%     avoid edge effects when contouring. Trial and error may be required. 
%     For n=150, err=0.005 works well.
%   opts 
%     0 = Return density index (Vollmer 2020).
%     1 = Return intensity index (Lisle 1985, Vollmer 1990).
%     2 = Return expected density index, expect must be protolith 
%         eigenvalues.
%     3 = Return expected intensity index, expect must be protolith 
%         eigenvalues (Hunter 2014).
%
  r = 0.0;
  p = 0.0;
  g = 0.0;
  [r, p, g] = xytotri(x, y, false, err); % apex down
  if (r + p + g) == 0.0 % not on plot
    d = NaN;
  else
    e = pgrtoeigen(p, g, r);
    if opts > 1 % expected protolith , 2 or 3
      e1 = e(1) - expect(1);
      e2 = e(2) - expect(2);
      e3 = e(3) - expect(3);
    else
      e1 = e(1) - 1/3;
      e2 = e(2) - 1/3;
      e3 = e(3) - 1/3;
    end;
    ss = e1*e1 + e2*e2 + e3*e3;
    if rem(opts, 2) % 1 or 3
      d = 7.5 * ss; % intensity 
    else    
      d = sqrt(1.5 * ss); % density 
    end
  end
end

function [e] = pgrtoeigen(p, g, r)
% PGRTOEIGEN  Converts PGR indexes to eigenvalues.
  e(3) = r/3.0;
  e(2) = g*0.5 + e(3);
  e(1) = p + e(2);
end

function [a,b,c] = xytotri(x, y, apexUp, err)
% XYTOTRI  Converts [x,y] to triangular coordinates.
%   Converts [x,y] in the range [-1,-1] to [1,1] to triangular coordinates, 
%   [a,b,c]. The plot is scaled to a height of 1.5 to give a unit 
%   circumcircle with the origin at [0,0]. [0,0,0] is returned for 
%   coordinates lie outside the triangle, allowed error is given by err. }
%   apexUp = Whether the plot has an apex at the top. For apexUp = true 
%     the apexes are:
%       a = [0,1]
%       b = [-sqrt(3)/2,-1/2]
%       c = [+sqrt(3)/2,-1/2]
%     For apexUp = false the y values are negated. 
%   err = Allowed error in a+b+c=1 while still returning a value. Used to 
%     prevents edge effects during, and may require adjustment for a given 
%     plot. 
%
  a = 0;
  b = 0;
  c = 0;
  rt3 = 1.732050807568877; % sqrt 3
  xt = x;
  if apexUp
    yt = y;
  else
    yt = -y;
  end
  xt = xt / 1.5;
  yt = (yt + 0.5) / 1.5;
  at = yt;
  if (at > 1.0+err) || (at < -err)
    return
  end
  bt = -(xt * rt3 - 1.0 + at) * 0.5;
  if (bt > 1.0+err) || (bt < -err) 
    return
  end
  ct = 1.0 - at - bt;
  if (ct > 1.0+err) || (ct < -err)
    return  
  end
  a = at;
  b = bt;
  c = ct;
end

