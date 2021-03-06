% File    : trifabtest.m
% System  : MATLAB
% Purpose : Test program for trifabplot.m and trifabdensity.m
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
%   Programs. Vol 52, No. 6, doi: 10.1130/abs/2020AM-358862.
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
% where [e1,e2,e3] are the normalized eigenvalues of the fabric 
% orientation tensor, and weight is a value used to interpolate the symbol 
% colors. The example data files use shear strain (olbary-Hansen2014.csv) 
% and a domain (0 to 4) based on distance into a shear zone 
% (qtz-Hunter2019.csv) weights.
% 
% Options (opts=0..3) include contours of density, intensity, expected 
% (protolith) density, and expected (portolith) intensity.
%-------------------------------------------------------------------------

% get comma delimited test file of e1,e2,e3,weight
% get comma delimited test file
[filename, pathname] = uigetfile( {'*.csv'});
eig = csvread([pathname, filename]);

% select one of two example files
%eig = csvread(['olbary-Hansen2014.csv']); % olivine barycenters
%eig = csvread(['qtz-Hunter2019.csv']); % quartz c-axes

% create a linear color map from white to red
lc = 256;
c1 = [1,1,1];
c2 = [1,0,0];
clin = [linspace(c1(1),c2(1),lc)', linspace(c1(2),c2(2),lc)', ... 
        linspace(c1(3),c2(3),lc)'];
    
% interpolate weights to colors 
wt = eig(:,4); % weights 
cwt = round((wt/max(wt))*(lc-1)+1); % scale weights to 1..lc
lwt = length(wt);
c = ones(lwt,3);
c(:,:) = clin(cwt(:),:); % get color scaled to weight

% get point coordinates and frame
[pgr,points,frame] = trifabplot(eig);

% get grid of the density index
expect = eig(1:1,1:3); % set expected for options 2 and 3 
ng = 150; % number of grid nodes
err = 0.005; % adjust if lines do not end at frame 
opt = 0; % options are 0 to 3
[x,y,z] = trifabdensity(ng, err, opt, expect); 

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
scatter(px, py, 72, c, 'filled', 'MarkerEdgeColor', 'k');

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
