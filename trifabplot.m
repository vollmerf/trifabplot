function [pgr, points, frame] = trifabplot(eig)
% TRIFABPLOT  Returns Cartesian coordinates of symbols for a triangular 
% fabric plot (Vollmer 1989, 1990).
%
% Input
%   eig    = Vector of normalized eigenvalues of fabric orientation 
%            matrixes, with [e1,e2,e3] as rows.
%
% Output
%   pgr    = PGR indexes of the eigenvalues.
%   points = [x,y] Cartesian cordinates of triangular plot symbols.
%   frame  = Coordinates of the triangular frame enclosed within a unit 
%            circumcircle centered at [0,0].
%
% Syntax
%    [pgr,points,frame] = trifabplot(eig);

% END HELP
% File    : trifabplot.m
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
%   Programs. Vol 52, No. 6, doi: 10.1130/abs/2020AM-358862.
% Vollmer F.W., 1990. An application of eigenvalue methods to structural 
%   domain analysis. Geological Society of America Bulletin, v. 102, n. 6,
%   786?791, ISSN 00167606
% Vollmer F.W., 1989. A triangular fabric plot with applications for 
%   structural analysis. EOS Transactions American Geophysical Union 
%   70:463
%
% One or more should be cited for usage of this or derivative code.
%-------------------------------------------------------------------------

  n = length(eig);
  for i = 1:n
    pgr(i,:) = eigentopgr(eig(i,:));
    points(i,:) = pgrtoxy(pgr(i,:));
  end
  frame = triframe(false);
end

function [pgr] = eigentopgr(eig)
% EIGENTOPGR  Converts orientation tensor eigenvalues to PGR indexes. 
%   The eigenvalues are assumed to be normalized.
  pgr(1) = eig(1)-eig(2);
  pgr(2) = 2.0 * (eig(2)-eig(3));
  pgr(3) = 3.0 * eig(3);
end

function [xy] = triframe(apexUp) 
% TRIFRAME  Returns the coordinates of the plot triangle. 
%   Returns the coordinates of a equilateral triangle that is 
%   circumscribed by a unit circle centered at the origin.
  r3 = 1.732050807568877; % sqrt 3
  xy(1,1) = -0.5 * r3; % top left
  xy(1,2) = 0.5;
  xy(2,1) = 0.5 * r3; % top right
  xy(2,2) = 0.5;
  xy(3,1) = 0.0; % bottom
  xy(3,2) = -1.0;
  if apexUp
    xy(:,2) = -xy(:,2)
  end
end
  
function [xy] = tritoxy(abc, apexUp)
% TRITOXY  Convert triangular coordinates to Cartesian. 
%   Plot is scaled to a height of 1.5 and centered in a unit circumcircle 
%   with the origin at [0,0]. a and b are assumed to be normalized, so 
%   a+b+c = 1 where [a,b,c] = abc. 
  irt3 = 0.577350269189626; % 1/(sqrt 3);
  xy(2) = abc(1); 
  xy(1) = (1.0 - abc(1) - 2.0 * abc(2)) * irt3;
  % scale to 1.5 height for unit circumcircle, origin at circumcenter
  xy(1) = xy(1) * 1.5;
  xy(2) = xy(2) * 1.5 - 0.5;
  if ~apexUp 
    xy(2) = -xy(2); 
  end
end

function [xy] = pgrtoxy(pgr)
% PGRTOXY  Convert triangular PGR coordinates to Cartesian. 
  t(1) = pgr(3); % a = r, bottom
  t(2) = pgr(1); % b = p, top left
  t(3) = pgr(2); % c = g, top right
  xy = tritoxy(t, false);
end

