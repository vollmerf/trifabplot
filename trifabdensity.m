function [x,y,z] = trifabdensity(n, err, opts, expect)
% TRIFABDENSITY  Returns a grid of density for a triangular fabric plot.
%   Returns a grid as [x,y,z] where [x, y] values are in the range [-1,-1] 
%   to [1,1]. [x,y] are converted to [P,G,R] indexes (Vollmer 1989, 1990), 
%   and eigenvalues to calculate the fabric density index, D 
%   (Vollmer 2020), in [z]. 
%
% Input
%   n = Number of grid nodes calculated in [x,y] between [-1,-1 and [1.1].
%   err = Error allowed in PGR indexes at plot margins, is used to 
%     avoid edge effects when contouring. For n=150, err=0.005 works,
%     adjustment may be required.
%   opts = Sets returned index:
%     0 = Density
%     1 = Intensity
%     2 = Expected density, expect must be protolith eigenvalues
%     3 = Expected intensity, expect must be protolith eigenvalues
%
% Output
%   [x,y] = Values in the range [-1,-1] to [1,1]
%   [z]   = Fabric density or intensity index if [x,y] is on triangle, 
%           otherwise NaN. 
%
% Syntax
%   [x,y,z] = trifabdensity();
%   [x,y,z] = trifabdensity(150,0.005,1);
%   [x,y,z] = trifabdensity(159,0.005,2,[0.2,0.2,0.6]);

% END HELP
% File    : trifabdensity.m
% System  : MATLAB
% Purpose : Triangular fabric (PGR or Vollmer) plots.
% Author  : Frederick W. Vollmer
% Date    : Aug 7, 2020
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


  switch nargin
    % n, err, opts, expect
    case 0
      n = 150;
      err = 0.005;
      opts = 0;
      expect = [1/3,1/3,1/3];
    case 1
      err = 0.005;
      opts = 0;
      expect = [1/3,1/3,1/3];
    case 2
      opts = 0;
      expect = [1/3,1/3,1/3];
    case 3
      expect = [1/3,1/3,1/3];
    case 4
      opts = opts;
    otherwise
      return
  end  
  x = zeros(n,n);
  y = zeros(n,n);
  z = zeros(n,n);
  inc = 2.0/(n-1); % [-1,-1] to [1,1]
  xi = -1.0;
  for i = 1:n
    yi = -1.0;
    for j = 1:n
      zi = densityXY(xi, yi, err, opts, expect); 
      x(i,j) = xi;
      y(i,j) = yi;
      z(i,j) = zi;
      yi = yi + inc;
    end
    xi = xi + inc;
  end
end

function [d] = densityXY(x, y, err, opts, expect)
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
  r = 0.0;
  p = 0.0;
  g = 0.0;
  [r, p, g] = XYToTri(x, y, false, err); % apex down
  if (r + p + g) == 0.0 % not on plot
    d = NaN;
  else
    e = PGRToEigen(p, g, r);
    if opts > 1 % expected, 2 or 3
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
    else % 0 or 2   
      d = sqrt(1.5 * ss); % density 
    end
  end
end

function [e] = PGRToEigen(p, g, r)
% PGRTOEIGEN  Converts PGR indexes to eigenvalues.
  e(3) = r/3.0;
  e(2) = g*0.5 + e(3);
  e(1) = p + e(2);
end

function [a,b,c] = XYToTri(x, y, apexUp, err)
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
%   err = Allowed error in a+b+c=1 to return a value. Used to prevent edge 
%     effects during, and may require adjustment for a given plot. 
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
