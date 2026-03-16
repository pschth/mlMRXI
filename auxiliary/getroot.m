function rootdir = getroot()
% returns project root directory. Requirement: initFunctions.m has to be
% called first

rootdir = strsplit( which('initFunctions.m'),  'initFunctions.m' ); 
rootdir = rootdir{1};

end

