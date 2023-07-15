arInit;
arLoadModel('model_template');
arLoadData('Tumor_CD_data');
arLoadData('Tumor_HFD_data');
arLoadData('Fat_CD_data');
arLoadData('Fat_HFD_data');
arCompileAll;

%% Constraint parameters
ar.lb = -7*ones(size(ar.lb)); % lower parameter bounds
ar.ub = 4*ones(size(ar.ub)); % upper parameter bounds

% set parameters for parameter estimation and optimization
ar.config.atol = 1e-8;
ar.config.rtol = 1e-8;
ar.config.optim.TolFun = 1e-8;
ar.config.optim.TolX = 1e-8;
ar.config.maxsteps = 1e5;

%% Numerical settings

% Optimizer settings
ar.config.optim.PrecondBandWidth = inf;
ar.config.optim.Display          = 'iter';
ar.config.optim.MaxIter          = 1e4;
