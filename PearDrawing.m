

model = createpde;
%importGeometry(model,'BracketTwoHoles.stl');
geometryFromEdges(model,@lshapeg);
mesh = generateMesh(model);
pdeplot(model)