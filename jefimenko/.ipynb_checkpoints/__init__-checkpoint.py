"""jefimenko - A EM simulator based on the Jefimenko equations"""

__version__ = '.1'
__author__ = 'Nehemiah Null <Bombadil224@gmail.com>'
__all__ = []

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import sys
from tqdm import tqdm

from .main import *
from .classes import *
from .graphing import *
from .simulation import *


print('loaded')