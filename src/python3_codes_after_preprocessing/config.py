import os

abspp = os.getcwd()
parent = os.path.abspath(os.path.join(abspp, os.pardir))
parent2 = os.path.abspath(os.path.join(parent, os.pardir))


PROJECT_ROOT = parent2
