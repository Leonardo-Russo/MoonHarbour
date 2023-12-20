close all
clear
clc


obj = renderSTL('ISS_2016.stl');

N = size(obj.Vertices, 1);

r = [1000, 0, 0]';

input('Press Enter to start moving.\n')
for i = 1 : N

    obj.Vertices(i, :) = obj.Vertices(i, :) + r';

end



function obj = renderSTL(fileName)

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
fv = stlread(fileName);

%All code below is for rendering the 3D object
obj = patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
view([-135 35]);

end