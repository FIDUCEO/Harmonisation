"""
HIRS sensor data
"""

'''___Python Modules____'''

'''___Third Party Modules____'''

'''___NPL Modules___'''
from HIRS_config_functions import return_sensor_data


'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "17/01/2018"
__credits__ = ["Gerrit Holl"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

'''___Constant___'''
channel_num = 9

sensor_data = return_sensor_data(channel_num)

if __name__ == "__main__":
    pass
