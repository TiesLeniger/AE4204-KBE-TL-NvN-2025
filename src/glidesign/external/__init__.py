# import, initialise MATLAB engine
import matlab.engine
MATLAB_Q3D_ENGINE = matlab.engine.start_matlab()
MATLAB_Q3D_ENGINE.cd(r'Q3D')