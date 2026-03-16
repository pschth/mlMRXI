function parsave( location, saveStruct )
% rewritten save function for parfor use
% location: specifies save location
% saveStruct: matlab struct containing all variables that are to be saved.

save(location,'-struct','saveStruct');

end

