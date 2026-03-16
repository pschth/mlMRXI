function phantomNew = interpolatePhantomToNewResolution(phantom, ResOld, ResNew)
% function name self-explanatory (interp does not preserve voxel contents
% correctly). Works for 3D and for cuboidal ROIs.
% % % 
% INPUT
% phantom ... matrix with resolution ResOld, contains voxel contents
% ResOld ... 1 x 3 vector, old resolution
% ResNew ... 1 x 3 scalar/vector, new resolution

% reshape phantom to original resolution
phantom = reshape(phantom, ResOld);

% store original particle mass
particleMass = sum(phantom(:));

% for each dimension scale up to least common multiple and scale down to
% new resolution
for d = 1:3
    leCoMu = lcm(ResOld(d),ResNew(d));
    upscaleFactor = leCoMu/ResOld(d);
    downscaleFactor = leCoMu/ResNew(d);
    
    if upscaleFactor~=1
        % scale up is least common multiple is different from old
        % resolution
        phantomUpscaled = repmat(phantom, [1 1 1 upscaleFactor]);
        switch d
            case 1
                phantomUpscaled = reshape(permute(phantomUpscaled, [4 1 2 3]), ...
                    leCoMu, ResOld(2), ResOld(3));
            case 2
                phantomUpscaled = reshape(permute(phantomUpscaled, [1 4 2 3]), ...
                    ResNew(1), leCoMu, ResOld(3));
            case 3
                phantomUpscaled = reshape(permute(phantomUpscaled, [1 2 4 3]), ...
                    ResNew(1), ResNew(2), leCoMu);
        end
    else
        % otherwise use original phantom resolution
        phantomUpscaled = phantom;
    end
    
    % piecewise summation of the least common multiple grid to gain new
    % resolution
    if downscaleFactor~=1
        switch d
            case 1
                phantom = squeeze(sum(reshape(phantomUpscaled, ...
                    downscaleFactor, ResNew(1), ResOld(2), ResOld(3)), 1));
            case 2
                phantom = squeeze(sum(reshape(phantomUpscaled, ...
                    ResNew(1), downscaleFactor, ResNew(2), ResOld(3)), 2));
            case 3
                phantom = squeeze(sum(reshape(phantomUpscaled, ...
                    ResNew(1), ResNew(2), downscaleFactor, ResNew(3)), 3));
        end
    else
        phantom = phantomUpscaled;
    end
    
end
% scale new phantom values according to new resolution
phantomNew = phantom./sum(phantom(:)) * particleMass; % scale to new voxel masses
end