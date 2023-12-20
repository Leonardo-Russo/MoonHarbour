% Load the STL file
[vertices, ~] = stlread('Torus.stl');
gm = geometryFromMesh('triangulation', delaunayTriangulation(vertices));

% Plot original geometry
figure;
pdegplot(gm, 'FaceLabels','on');
title('Original Geometry');

% Calculate the bounding box and its centroid
minValues = min(vertices);
maxValues = max(vertices);
centroid = (minValues + maxValues) / 2;

% Define your desired center position
desiredCenter = [10, 10, 0];

% Calculate the translation vector
translationVector = desiredCenter - centroid;

% Translate the vertices
translatedVertices = vertices + translationVector;

% Create new geometry from translated vertices
newGm = geometryFromMesh('triangulation', delaunayTriangulation(translatedVertices));

% Plot translated geometry
figure;
pdegplot(newGm, 'FaceLabels','on');
title('Translated Geometry');
