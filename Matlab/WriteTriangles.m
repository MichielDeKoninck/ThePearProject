
%% Punten gegenereerd uit PDE-modeler inladen
load("PointsRefinementLevel0.mat") % p
load("PointRefinementLevel1.mat") %p1

%% Test: Willekeurige Punten uitlezen
p(1,180)
p1(1,185)

%% Lijst van driehoek gebruiken om bestand aan te maken met als info: Driehoeknummer, Coordinaten van hoekpunten

%triangles, edges, points: krijg je als data in je workspace door mesh te
%genereren in PDE modeler en vervolgens te exporteren. 

trianglePointIndices = triangles(1:3,:); %Nodig welke 3 punten verbonden worden
edgesPointIndices = edges(1:2,:); %Enkel nodig welke 2 punten ze verbinden


AmountofPoints = length(points(1,:));
ListOfIndexesPoints = 1:1:AmountofPoints;
A_points = [ListOfIndexesPoints; points(1,:); points(2,:)];

AmountofEdges = length(edges(1,:)); %Dit zijn enkel de edges die op de rand van de figuur liggen. Niet alle edges van alle driehoeken.
ListOfIndexesEdges = 1:1:AmountofEdges;
A_edges = [ListOfIndexesEdges; edges(1,:); edges(2,:)];

AmountofTriangles = length(triangles(1,:)); 
ListOfIndexesTriangles = 1:1:AmountofTriangles;
A_triangles = [ListOfIndexesTriangles; triangles(1,:); triangles(2,:);triangles(3,:)];


fileID = fopen('pointsinfo.txt','w');
fprintf(fileID,'%6s %10s %10s\n','PointNumber','x','y');
fprintf(fileID,'%d %20f %13f\n',A_points );
fclose(fileID);


fileID = fopen('edgesinfo.txt','w');
fprintf(fileID,'%6s %10s %11s\n','Edgenumber','Point1','Point2');
fprintf(fileID,'%d %17d %11d\n',A_edges);
fclose(fileID);


fileID = fopen('trianglesinfo.txt','w');
fprintf(fileID,'%6s %6s %6s %6s\n','TriangleNumber','Punt1','Punt2', 'Punt3');
fprintf(fileID,'%12d %6d %6d %6d\n',A_triangles );
fclose(fileID);





