% File    : pgrtest.m
% System  : MATLAB
% Purpose : Test program for pgrplot.m and pgrdensity.m
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
% One or more of should be cited for usage of this or derivative code.
%
% The example data sets are from:
%
% Hansen L.N., Zhao Y.H, Zimmerman M.E., Kohlstedt D.L, 2014. Protracted 
%   fabric evolution in olivine: Implications for the relationship among 
%   strain, crystallographic fabric, and seismic anisotropy. Earth and 
%   Planetary Science Letters 387:157?168, ISSN 0012821X
% Hunter N.J., Weinberg R., Wilson C., Luzin V., Misra S., 2019. Quartz 
%   deformation across interlayered monomineralic and polymineralic rocks: 
%   A comparative analysis. Journal of Structural Geology v. 119, p. 
%   118?134, ISSN 01918141
%
% Notes
% -----
% The input data file is a csv file consisting of rows:
%
%   e1,e2,e3,weight
%
% where [e1,e2,e3] are the normalized eigenvalues of the orientation
% tensor, and weight is a value used to determine the symbol color.
%-------------------------------------------------------------------------

% get comma delimited test file of e1,e2,e3,weight
% select one of two example files
eig = csvread(['olbary-Hansen2014.csv']); % olivine barycenters
%eig = csvread(['qtz-Hunter2019.csv']); % quartz c-axes

wt = eig(:,4); % weights for symbol color

% create a linear color map from white to red
lc = 256;
c1 = [1, 1, 1];
c2 = [1, 0, 0];
clin = [linspace(c1(1),c2(1),lc)', linspace(c1(2),c2(2),lc)', ... 
        linspace(c1(3),c2(3),lc)'];
% map weights to colors 
cwt = round((wt/max(wt))*(lc-1)+1); % scale weights to 1..lc
lw = length(wt);
clr = ones(lw,3);
for i=1:lw
  clr(i,:) = clin(cwt(i),:); % get color scaled to weight
end;     

% get point coordinates and frame
[pgr, points, frame] = pgrplot(eig);

% get grid of the density index
expect = eig(1:1,1:3); % set expected, used for options 2 and 3 
ng = 150; % number of grid nodes
err = 0.005; % adjust if contour lines do not end at frame 
opt = 0; % options are 0 to 3
[x,y,z] = pgrdensity(ng, err, opt, expect); 

% set up figure
figure;
hold on;
axis([-1.0 1.0 -1.0 1.0]);
axis('equal');
axis('off');

% contour grid at D = 0.1
contour(x,y,z,9,'Color',[.4 .4 .4]);

% plot points
px = points(:,1); 
py = points(:,2); 
scatter(px, py, 72, clr, 'filled', 'MarkerEdgeColor', 'k');

% plot triangular frame, array of [x,y]
for i = 1:3
  j = i+1;
  if j == 4 
    j = 1;
  end
  lx = [frame(i,1), frame(j,1)];
  ly = [frame(i,2), frame(j,2)]; 
  line('XData', lx, 'YData', ly, 'Color', 'k', 'LineWidth', 2);
end

% plot labels
text(-0.98,0.54,'P','FontSize',24,'HorizontalAlignment', 'Center')
text(0.98,0.54,'G','FontSize',24,'HorizontalAlignment', 'Center')
text(-0.00,-1.14,'R','FontSize',24,'HorizontalAlignment', 'Center')

hold off;
