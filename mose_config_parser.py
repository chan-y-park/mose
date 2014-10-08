import ConfigParser
import sympy as sym
from math import pi


class MoseConfigParser(ConfigParser.SafeConfigParser):
    """
    This simple class is a wrapper of SafeConfigParser, and returns
    a value of an option as a Python expression,
    whereas RawConfigParser returns a value as a string.
    """
    def get(self, section, option):
        u = sym.Symbol('u')
        value = ConfigParser.SafeConfigParser.get(self, section, option)
        return eval(value)
