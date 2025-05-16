# Python native imports
import os

# Python third party package imports

# Parapy imports
from parapy.geom import FittedCurve, Point
from parapy.core import Attribute, Part, Input
from parapy.core.validate import OneOf, Range
from kbeutils.data import airfoils

class Airfoil(FittedCurve):

