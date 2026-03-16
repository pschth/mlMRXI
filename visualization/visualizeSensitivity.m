function visualizeSensitivity( sourceGrid, S, figNr )

zslice = unique(sourceGrid.positions(:,3));
figure(figNr);
slice(sourceGrid.originalMesh.X, sourceGrid.originalMesh.Y, sourceGrid.originalMesh.Z, reshape(S, sourceGrid.dimensions), [], [], zslice);
colorbar;
alpha(0.75);
view([10 -9]);

end

