%Code om matrix uit matlab in een bestand te schrijven waar C++ even
%gemakkelijk zal aankunnen. 
%Nog niet volledig zeker dat dit op de juiste manier wegschrijft.

%data = rand(8, 10);
data = p
fid  = fopen('FileRefinementlevel1.txt', 'w');
if fid == - 1
  error('Cannot open file for writing');
end
fwrite(fid, ndims(data), 'uint16'); % TODO nog niet zeker van 
fwrite(fid, size(data), 'uint64');
fwrite(fid, data, 'double');
fclose(fid);