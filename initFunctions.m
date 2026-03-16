function root = initFunctions()
% adds all folders to the MATLAB path; has to be called in every new
% project

%%%%%% STANDARD PATHS:
root = pwd;
addpath(root);
addpath([root '/auxiliary/']);
addpath([root '/auxiliary/external/phantom3D/']);
addpath([root '/auxiliary/external/subaxis/']);
addpath([root '/auxiliary/external/regu/']);
addpath([root '/auxiliary/external/triangle_ray_intersection/']);
addpath([root '/auxiliary/optimization/']);
addpath([root '/reconstruction/']);
addpath([root '/setup/']);
addpath([root '/setup/generation/']);
addpath([root '/setup/calculation/']);
addpath([root '/visualization/']);
%%%%%%

%%%%%% CUSTOM PATHS:
% custom paths are to be added here if demanded
addpath([root '/experimental_data/']);
addpath([root '/experimental_data/201604_PTB_P_Phantom/']);
addpath([root '/experimental_data/201901_PTB_PTB_Phantom/']);
addpath([root '/experimental_data/202006_6KanalSystem/']);
%%%%%%
end

