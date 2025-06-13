# import, initialise MATLAB engine
import matlab.engine
MATLAB_Q3D_ENGINE = matlab.engine.start_matlab()
MATLAB_Q3D_ENGINE.cd(r'src/glidesign/external/Q3D')

from .q3d_data import Q3DData