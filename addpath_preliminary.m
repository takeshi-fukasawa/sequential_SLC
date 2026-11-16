code_path='C:/Users/fukas/Dropbox/git/';
save_path="C:/Users/fukas/Dropbox/sequential/output/";


sequential_code_path=append(code_path,'sequential_SLC_temp/');

addpath(append(sequential_code_path,'sequential_algorithm/'))

sequential_code_path2=append(sequential_code_path,'sequential_algorithm/');
addpath(append(sequential_code_path2,'SLC'))
addpath(append(sequential_code_path2,'EPL'))
addpath(append(sequential_code_path2,'numerical_derivative'))
addpath(append(sequential_code_path2,'Newton'))
addpath(append(sequential_code_path2,'NFXP'))    
addpath(append(sequential_code_path2,'linear_eq'))
addpath(append(sequential_code_path2,'util'))


%% Fixed Point Iteration Acceleration Methods
FPI_acceleration_methods_path=append(code_path,'spectral/');
addpath(FPI_acceleration_methods_path)
