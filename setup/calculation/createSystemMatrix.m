% calculates the MRX forward model system matrix "out". Also returns the
% magnetic field vectors H_mag at each voxel location for each activation
% sequence.
% 
% SETUP: standard MRXI config file from the MRXI toolbox
% 
% I: nCoils x nActivations matrix, contains coil current values per
% activation. Default: Identity-matrix
% 
% thresh_RAM: a RAM-threshold in Gb can be stated. If the
% predicted RAM-usage of the calculation is below, the fast version of the
% algorithm is computed. Otherwise, the slower, RAM-efficient version is
% computed. Default: 4 Gb
% 
% verbose: boolean; if false, no text messages are generated. Default: true
% 


%   Author:     Peter Schier
%   E-Mail:     peter.schier@umit.at
%   Institute:  UMIT Hall in Tirol
%   Date:       13-Oct-2017

function [out, H_mag] = createSystemMatrix( SETUP, I, thresh_RAM, verbose )
    
    if nargin < 2
        I = eye(SETUP.getNumberOfCoils);
    elseif isempty(I)
        I = eye(SETUP.getNumberOfCoils);
    end

    if nargin < 3
        thresh_RAM = 4;
    elseif isempty(thresh_RAM)
        thresh_RAM = 4;
    end
    
    if nargin < 4
        verbose=true;
    elseif isempty(verbose)
        verbose=true;
    end
    
    H_mag = getMagneticField( SETUP, I, verbose );
    out = getLeadFieldMatrix( SETUP, H_mag, thresh_RAM, verbose );
end
