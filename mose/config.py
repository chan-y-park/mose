import ConfigParser
import logging

from math import pi


class MoseConfigParser(ConfigParser.SafeConfigParser):
    """
    This simple class is a wrapper of SafeConfigParser.
    """
    def get(self, section, option):
        """
        Returns a value of an option as a Python expression,
        whereas ConfigParser.get() returns a value as a string.
        """
        value = ConfigParser.SafeConfigParser.get(self, section, option)
        return eval(value)


    def getstr(self, section, option):
        """
        Returns a value as a string.
        """
        return ConfigParser.SafeConfigParser.get(self, section, option)


class MoseConfig:
    """
    A container class of the configuration data.
    Saves the configuration data as a Python dictionary.
    """
    def __init__(self):
        self.data = {}
        self.parser = None 


    def __setitem__(self, option, value):
        self.data[option] = value


    def __getitem__(self, option):
        try:
            value = self.data[option]
        except KeyError as e:
            logging.warning('Option {} not specified; use None.'.format(e))
            self.data[option] = None
            value = None
        return value


    def iteritems(self):
        return self.data.iteritems()


    def read(self, config_file):
        """
        Read an .ini file and load the configuration data.
        """
        logging.info('config file: %s', config_file)
        config_parser = MoseConfigParser()

        with open(config_file, 'r') as fp:
            config_parser.readfp(fp)

        for section in config_parser.sections():
            self[section] = {}
            for option in config_parser.options(section):
                if (section == 'fibration'):
                    self[section][option] = (
                        config_parser.getstr(section, option)
                    )
                else:
                    self[section][option] = (
                        config_parser.get(section, option)
                    )

        self.parser = config_parser
