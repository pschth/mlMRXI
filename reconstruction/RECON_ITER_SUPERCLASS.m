%%%% DO NOT USE! MUCH SLOWER THAN INDEPENDENT CLASSES!


classdef RECON_ITER_SUPERCLASS < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%% PUBLIC MEMBERS
    properties (Access = public)
        
        req; % required algorithm parameters
        aux; % non-mandatory auxiliary algorithm parameters
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE MEMBERS
    properties (Access = protected)
        
        int; % internal algorithm parameters
        info; % informational algorithm parameters
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% PUBLIC FUNCTIONS
    methods (Access = public)
        function o = RECON_ITER_SUPERCLASS(L,b)
            % constructor
            %
            % INPUT:
            % SETUP - MRXI setup configuration object
            %
            % L - nActivations*nSensors x nVoxels matrix; leadfield matrix
            %
            % b - nActivations*nSensor vector; measurement vector
            %
            % 
            
            if ~size(L,1)==length(b)
                error('Matrix dimension and measurement dimension do not match.');
            end
            
            o.setSystem(L,b);
            o = o.initCommon();
                        
        end
        
        function o = copySettings(o, RECON, forceCopy)
            % copies the required and auxiliary settings of an object of
            % the same RECON class
            %
            % INPUT:
            % RECON - object; RECON class object
            %
            % forceCopy (optional) - bool; forces the copying of settings although the
            % wrong RECON class is provided; default: false
            
            if nargin < 2
                forceCopy = false;
            end
            
            if isa(RECON, o.info.className) || forceCopy
                o.req = RECON.req;
                o.aux = RECON.aux;
            else
                error('Wrong RECON class.');
            end
            
        end
        
        function o = setSystem(o,L,b)
            % sets leadfield matrix and measurement vector. Use e.g. after
            % copying the object for use of different system with same
            % parameter settings.
            %
            % INPUT:
            % L (optional) - nActivations*nSensors x nVoxels matrix; leadfield matrix
            %
            % b (optional) - nActivations*nSensor vector; measurement vector
            
            if ~isempty(L)
                o.int.L = L;
            end
            if ~isempty(b)
                o.int.b = b;
            end
        end
        
        function printExplanation(o)
            
            % get all fieldnames
            namesReq = fieldnames(o.req);
            namesAux = fieldnames(o.aux);
            namesInfo= fieldnames(o.info);
            
            % get the longest name
            nLongestName = 0;
            for i = 1:length(namesInfo)
                if length(namesInfo{i})>nLongestName
                    nLongestName = length(namesInfo{i});
                end                
            end
            
            fprintf(['\nAlgorithm: ' o.info.name]);
            % print explanations of required variables
            fprintf('\n\nRequired variables:');
            for r = 1:length(namesReq)
                for i = 1:length(namesInfo)
                    if strcmp(namesReq{r}, namesInfo{i})
                        nFill = nLongestName - length(namesReq{r});
                        fillStr = repmat(' ', 1, nFill);
                        fprintf(['\n' namesReq{r} fillStr '\t' o.info.(namesReq{r})]);
                        continue;
                    end
                end
            end
            
            fprintf('\n\nAuxiliary variables:');
            % print explanations of auxiliary variables
            for a = 1:length(namesAux)
                for i = 1:length(namesInfo)
                    if strcmp(namesAux{a}, namesInfo{i})
                        nFill = nLongestName - length(namesAux{a});
                        fillStr = repmat(' ', 1, nFill);
                        fprintf(['\n' namesAux{a} fillStr '\t' o.info.(namesAux{a})]);
                        continue;
                    end
                end
            end
            fprintf('\n');
        end
        
        function int = getInternals(o)
            % returns the internal parameters of the object
            int = o.int;
        end
        
        function o = start(o)
            % starts the reconstruction algorithm
            
            o.int.tic = tic;
            o.reconstruct();
        end
        
        function [result, corrstorage] = getResult(o)
            % returns the result
            
            result = o.int.result;
            if nargout > 1
                corrstorage = o.int.corrstorage;
            end
        end
        
 
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTIONS
    methods (Access = protected)
        
        function o = initCommon(o)
            % initializes the common parameters of the reconstruction
            % algorithm.
            
            % set all known internal parameters
            o.initInternals();
            % set common informational parameters
            o.initCommonInfo();

            %%% required
            o.req.nIter                 = 1e4;                                  % Number of iterations
            o.req.deltaMin              = 0;                                    % in %, minimum normalized RMS deviation from one iteration step to the next, otherwise abort
            
            %%% auxiliary
            o.aux.xInit                 = zeros( o.int.nVox, 1);    % Estimation of the result (zeros = cold start, estimation = warm start)
            o.aux.xOrig                 = zeros( o.int.nVox, 1);    % Ground truth for correlation control during computation
            o.aux.OutputInterval        = 100;                      % Defines the verbose interval
            o.aux.verbose               = true;
        end
        
        function o = init(o)
            % initializes the individual parameters of the reconstruction
            % algorithm
        end
        
        function o = initInternals(o)
            % defines all required internal parameters known after calling
            % the constructor
            
            % get setup constants
            o.int.nVox = size(o.int.L,2);
            o.int.rnVox = sqrt(o.int.nVox);
            
            % set algorithm parameters
            o.int.i = 0; % algorithm iteration count
            o.int.vi = 0; % verbose count
            
            % initialize deviation variable
            o.int.dev = inf;

        end  
        
        function o = initCommonInfo(o)
            % defines all common informational parameters known after
            % calling the constructor
            
            o.info.runtime = 0; % algorithm runtime
            
            %%% informational (field names must be identical to required
            %%% field names)
            o.info.nIter = 'Scalar - Number of iterations';
            o.info.deltaMin = 'Scalar - Minimum deviation of one iteration step to the next. Otherwise abort.';
            %%% informational (field names must be identical to auxiliary
            %%% field names)
            o.info.xInit = 'nVoxel vector - Initialisation of the result (zeros = cold start, estimation = warm start)';
            o.info.xOrig = 'nVoxel vector - Ground truth for correlation comparison during computation (zeros = no comparison, otherwise = comparison)';
            o.info.OutputInterval = 'Scalar - Defines the verbose interval';
        end
        
        function o = initInfo(o)
            % defines all individual informational parameters known after
            % calling the constructor
        end
        
        function o = reconstruct(o)
            % this is the actual reconstruction algorithm.
        end
        
        function o = printProgress(o)
            if o.aux.verbose
                if o.int.i >= o.int.vi || o.int.i == o.req.nIter || o.int.dev < o.req.deltaMin
                    % update runtime
                    o.int.toc = toc(o.int.tic);
                    o.info.runtime = o.info.runtime + o.int.toc;
                    o.int.tic = tic;

                    % print informational stats
                    if isempty(o.aux.xOrig) || all(o.aux.xOrig(:)==0)
                        fprintf([o.info.name '\tIteration %.0f/%.0f\tDeviationPrev: %.2e/%.2e\tRuntime: %.2f\n'],...
                            o.int.i, o.req.nIter, o.int.dev, o.req.deltaMin, o.info.runtime);
                    else
                        currcorr = corr(o.aux.xOrig(:), o.int.result(:));
                        if o.int.i == 1
                            o.int.fillIdx = 1;
                        else
                            o.int.fillIdx = o.int.fillIdx + 1;
                        end
                        o.int.corrstorage(2,o.int.fillIdx) = currcorr;
                        o.int.corrstorage(1,o.int.fillIdx) = o.int.i;
                        if isnan(currcorr)
                            o.int.i = o.req.nIter;
                        end
                        fprintf([o.info.name '\tIteration %.0f/%.0f\tDeviationPrev: %.2e/%.2e\tRuntime: %.2f\tCorrelation: %.5f\n'],...
                            o.int.i, o.req.nIter, o.int.dev, o.req.deltaMin, o.info.runtime, currcorr);
                    end

                    % increase verbose counter
                    o.int.vi = o.int.vi + o.aux.OutputInterval;
                end
            end
        end
        
    end
end

