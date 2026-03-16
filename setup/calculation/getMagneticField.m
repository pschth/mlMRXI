function H = getMagneticField( mrxisetup, I, verbose )
% calculates the magnetic field vectors in the voxels that result from the
% driving fields of the electromagnetic coils defined in mrxisetup. The
% coils are driven by the currents defined in the coil current matrix I.
% 
% INPUT:
% mrxisetup - MRXI setup configuration object
% 
% I - nCoils x nActivations matrix; contains the coil currents in [A] per
% activation; if left blank, consecutive, single coil activations with unit
% current are calculated (i.e., I = eye(nCoils))
% 
% verbose (optional) - boolean; if false, no text messages are generated.
% Default: true
% 
%
% OUTPUT:
% H - (nActivations*nVoxels) x 3 matrix; magnetic field vector resulting in
% each voxel during every activation. The nActivation different activations
% are concatenated below each other.
    
    %%% exceptions and default values
    % get number of coils
    nCoils = mrxisetup.getNumberOfCoils;
    if nargin < 2 || isempty(I)
        % if no coil current matrix is defined or if it is empty, use
        % consecutive, single coil activations with unit current
        I = eye(nCoils);
    end
    
    % if no 'verbose' argument was provided, set default value TRUE
    if nargin < 3 || isempty(verbose)
        verbose = true;
    end
        
    % number of rows in I must equal nCoils
    if size(I,1) ~= nCoils
        error('Wrong size of coil current matrix I.');
    end
    

    %%% code
    % get number of voxels and number of activations
    nVoxels = mrxisetup.getNumberOfVoxels;
    nActivations = size(I,2);
    
    % check if there are coils that are never activated
    activeCoils = ~all(I==0,2);
    
    % use single coilpattern if no assignment was performed; 
    % INFO FOR pyMRXI: this is to catch errors from old code -> not
    % necessary when translating to Python
    if ~isfield(mrxisetup.coilData, 'coilpatternassignment')
        mrxisetup.coilData.coilpatternassignment = ones(nCoils,1);
    end

    % if coils are given as coil patterns, reshape coil segment positions
    % for matrix calculations
    if ~mrxisetup.coilData.isDipole
        % get different used coil patterns of active coils
        usedCoilPatterns = unique(mrxisetup.coilData.coilpatternassignment(activeCoils));
        % get number of different used coil patterns
        nUsedCoilPatterns = length(usedCoilPatterns);
        % get coil segment data for each applied coil pattern
        coilSegments = cell(nUsedCoilPatterns,1); % for storing the different coil patterns
        nSegments = zeros(nUsedCoilPatterns,1); % for storing the number of coil segments of the different coil patterns
        for cp = 1:nUsedCoilPatterns
            % get number of active coils using the current coil pattern
            nCoilsWithCurrentCoilPattern = sum( mrxisetup.coilData.coilpatternassignment(activeCoils) == ...
                                                usedCoilPatterns(cp));
            % get number of coil segment end points for current coil pattern
            % INFO FOR pyMRXI: the coilpatterns should consistently
            % available in cells, so the following if statement should not
            % be necessary, just the line with the first 'size' function
            if iscell(mrxisetup.coilData.coilpattern)
                nSegPoints = size(mrxisetup.coilData.coilpattern{usedCoilPatterns(cp)}, 2);
            else
                nSegPoints = size(mrxisetup.coilData.coilpattern, 2);              
            end
            nSegments(cp) = nSegPoints-1; % number of coil segments of current coil pattern
            coilSegments{cp} = zeros(3,nSegPoints,nCoilsWithCurrentCoilPattern); % initialize appropriate storage matrix for coilSegments cell cp
        end
        % fill coilSegments cell(s)
        coilCounter = ones(nUsedCoilPatterns,1); % counter for which coil pattern was used how many times
        % for every coil...
        for ci = 1:nCoils
            % skip coil if it is never used
            if ~activeCoils(ci)
                continue
            end
            % get the index of the used coil pattern in usedCoilPatterns
            cpFillIdx = find(usedCoilPatterns == mrxisetup.coilData.coilpatternassignment(ci));
            % get current coilData.segments field name
            structname = ['c' num2str(ci)];
            % address correct coil pattern cpFillIdx to copy the coil
            % segment data from the current coil inside at index
            % coilCounter(cpFillIdx)
            coilSegments{cpFillIdx}(:,:,coilCounter(cpFillIdx)) = mrxisetup.coilData.segments.(structname);
            % increment coilCounter(cpFillIdx)
            coilCounter(cpFillIdx) = coilCounter(cpFillIdx) + 1;
        end
        % decrement counter to actual number of coils used
        coilCounter = coilCounter - 1;
    end
    
    
    % initialize all huge matrices for later to predict RAM usage
    Hvec = zeros([3 nVoxels sum(activeCoils)]); % for storing intermediate results
    H = zeros(nVoxels*nActivations, 3); % for the final result
            
    % if coils are approximated as dipoles...
    if mrxisetup.coilData.isDipole
        % print status message if verbose is TRUE
        if verbose
            fprintf('\n\tCalculating theoretical influences of each coil on each voxel...');
        end
        
        % initialize fill index for Hvec
        iActiveCoils = 0;
        % for all coils...
        for ci = 1:nCoils
            % the magnetic field of the dipoles is calculated according to:
            % https://en.wikipedia.org/wiki/Magnetic_dipole -> EQUATION (2)
            % NOTE: the magnetic moment is calculated here as: magnetic
            % moment = coil current * dipole orientation, thus the coil
            % current defined for dipoles is actually the strength of the
            % magnetic moment
            
            % skip coil if it is never used
            if ~activeCoils(ci)
                continue
            end
            
            % increment Hvec index if coil is used
            iActiveCoils = iActiveCoils + 1;
            
            % get all voxel-dipole distance vectors
            rVec = bsxfun(@minus, mrxisetup.ROIData.voxels, mrxisetup.coilData.centers(:,ci));
            
            % get all voxel-dipole scalar distances
            r = rssq(rVec,1);
            
            % for each coil ci get theoretical influences of the coil ci on
            % EVERY voxel (mrxisetup.coilData.orientations are the unit
            % magnetic moment vectors in this case; multiplied with I, they
            % become the final magnetic moment vectors)
            Hvec(:,:,iActiveCoils) = bsxfun(@times, 3*rVec, (mrxisetup.coilData.orientations(:,ci)'*rVec)./r.^5)-... left term in the Wikipedia equation
                           mrxisetup.coilData.orientations(:,ci)./r.^3; % right term in the Wikipedia equation
            % HINT: division by 4*pi and multiplication with I is performed
            % at the end of this function
        end
        % remove garbage
        clear rVec r
    else
%         the magnetic field of the coils is calculated according to:
%         Liebl, M., et al. "Quantitative imaging of magnetic nanoparticles
%         by magnetorelaxometry with multiple excitation coils." Physics in
%         Medicine & Biology 59.21 (2014): 6607. -> EQUATION (6)
        
        % for each coil pattern...
        for cp = 1:nUsedCoilPatterns
            
            % print status message if verbose is TRUE
            if verbose
                if nUsedCoilPatterns > 1
                    fprintf('\n\tCalculating theoretical influences of each coil on each voxel for coil pattern %.0d (%.0d/%.0d): ', ...
                        usedCoilPatterns(cp), cp, nUsedCoilPatterns); 
                else
                    fprintf('\n\tCalculating theoretical influences of each coil on each voxel: '); 
                end
                % initialize precentage for printing the progress of the
                % calculation
                percentage = 0;
                % print current progress
                delCount = fprintf('%.0f%%', percentage);
                pause(1e-9); % pause a nanosecond for printing (only for safety, should work without)
            end

            % for each voxel...
            for vi = 1:nVoxels
                
                % print progress in percent if verbose is TRUE
                if verbose
                % print current progress
                    if floor(vi/nVoxels*100) >= percentage
                        percentage = ceil(vi/nVoxels*100);
                        fprintf(1, repmat('\b',1,delCount)); % delete old percentage
                        delCount = fprintf('%.0f%%', percentage); % print new percentage
                        pause(1e-9); % pause a nanosecond for printing (only for safety, should work without)
                    end
                end

                % get distance vectors between the center of voxel vi and
                % the start points of the coil segments
                s1Vec = bsxfun(@minus, coilSegments{cp}, mrxisetup.ROIData.voxels(:,vi));
                % get distance vectors between the center of voxel vi and
                % the end points of the coil segments (second start point
                % is first end point and so on...)
                s2Vec = circshift(s1Vec, -1, 2);
                % delete the last entries since these correspond to the
                % first start point in s1Vec and to the first end point in
                % s2Vec (which are already covered by s1Vec(:,1,:) and
                % s2Vec(:,1,:))
                s1Vec(:,end,:) = [];
                s2Vec(:,end,:) = [];

                % get absolute distances of the distance vectors from voxel
                % center vi to the start and the end points of the coil
                % segments
                s1 = rssq(s1Vec,1);
                s2 = circshift(s1, -1, 2);
                % reshape to row vectors for later computations
                s1 = s1(:)';
                s2 = s2(:)';

                % reshape distance vectors to 3 x
                % (nSegments*nCoilsOfCurrentCoilPattern) (for faster cross
                % & dot product computation in MATLAB)
                s1Vec = reshape(s1Vec, [3 nSegments(cp)*coilCounter(cp)]);
                s2Vec = reshape(s2Vec, size(s1Vec));

                % get product of absolute distances of s1 and s2
                s1_times_s2 = s1 .* s2;

                % calculate influences of each coil on voxel vi with unit
                % current (term in the brackets in Liebl et al., equation
                % (6) without I)
                scalarArguments = (s1 + s2)./(s1_times_s2.*(s1_times_s2 + sum(s1Vec.*s2Vec,1))); % calculate everything inside the brackets that results in scalar values
                H_temp = bsxfun(@times, cross(s1Vec,s2Vec,1), scalarArguments); % calculate vectorial (cross) operation and multiply with the scalar values
                % HINT: division by 4*pi and multiplication with I is performed
                % at the end of this function
                
                % fill all entries of Hvec that correspond to voxel vi and
                % to the currently used coil pattern usedCoilPatterns(cp)
                Hvec(:,vi,mrxisetup.coilData.coilpatternassignment(activeCoils) == usedCoilPatterns(cp)) = ...
                    sum( ... sum up the contributions of all segments to get the magnetic field vectors of every used coil with the pattern usedCoilPatterns(cp) in voxel vi 
                    reshape( H_temp, [3 nSegments(cp) coilCounter(cp)] ), 2);
            end
        end
        % remove garbage
        clear s1Vec s2Vec s1 s2 s1_times_s2 scalarArguments H_temp
    end
    
    % print status message if verbose is TRUE
    if verbose
        fprintf('\n\tCalculating magnetic field(s) of sequence(s): '); 
        delCount = 0;
    end
    
    % for every coil activation sequence (=columns of matrix I)...
    for seq = 1:nActivations
        % print current progress if verbose is TRUE
        if verbose
            fprintf(1, repmat('\b',1,delCount)); % delete old sequence count
            delCount = fprintf('%.0f / %.0f', seq, nActivations); % print new sequence count
            pause(1e-9); % pause for printing
        end
        
        % calculate fill indices for H for the current activation sequence
        % seq
        iH = ((seq-1)*nVoxels+1) : (seq*nVoxels);
        
        % initialize fill index for Hvec
        iActiveCoils = 0;
        % for every coil...
        for ci = 1:nCoils
            % skip coil if it is never used
            if ~activeCoils(ci)
                continue
            end
            % increment Hvec index if coil is used
            iActiveCoils = iActiveCoils + 1;
            
            % multiply all coils with their respective coil currents
            % applied in the activation sequence 'seq' and add their
            % individual contributions to get the resulting magnetic field
            % vectors in all voxels
            H(iH, :) = H(iH, :) + (Hvec(:,:,iActiveCoils) .* I(ci,seq))';
        end
    end
    % divide H by 4*pi as required by equation (6) in Liebl et al.
    H = H./(4*pi);
    
    % done
    if verbose
        fprintf('\nDone!\n');
    end
end