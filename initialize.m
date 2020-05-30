function [voxel,dimensions,minMax,allData] = initialize(dataFile,nuclei,presynaptic,postsynaptic,extra)

%open datafile
data=bfopen(dataFile);

%get dimensions of data for preallocation
omeMeta = data{1, 4};
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices

%Get the voxel size for graphing and grouping later
voxelSizeXdefaultValue = omeMeta.getPixelsPhysicalSizeX(0).value();           % returns value in default unit
voxelSizeXdefaultUnit = omeMeta.getPixelsPhysicalSizeX(0).unit().getSymbol(); % returns the default unit type
voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
voxelSizeXdouble = voxelSizeX.doubleValue();                                  % The numeric value represented by this object after conversion to type double
voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER); % in µm
voxelSizeYdouble = voxelSizeY.doubleValue();                                  % The numeric value represented by this object after conversion to type double
voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER); % in µm
voxelSizeZdouble = voxelSizeZ.doubleValue();   

%save the metadata as array and return
voxel=[voxelSizeXdouble,voxelSizeYdouble,voxelSizeZdouble];
dimensions=[stackSizeX,stackSizeY,stackSizeZ];


[allData] = splitToMatrices(data,nuclei,presynaptic,postsynaptic,extra);




%Save mins and maxes for filtering
minMax=[min(allData(:,:,:,1),[],'all'),max(allData(:,:,:,1),[],'all');min(allData(:,:,:,2),[],'all'),max(allData(:,:,:,2),[],'all');min(allData(:,:,:,3),[],'all'),max(allData(:,:,:,3),[],'all');min(allData(:,:,:,4),[],'all'),max(allData(:,:,:,4),[],'all')];

end
