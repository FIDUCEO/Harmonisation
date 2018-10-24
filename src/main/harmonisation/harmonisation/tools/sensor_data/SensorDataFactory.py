"""
Factory for sensor data
"""

'''___Built-In Modules___'''
from collections import OrderedDict

'''___Third Party Modules___'''

'''___harmonisation Modules___'''
from AVHRR.AVHRR_11um_3_rescaled_sensor_data import AVHRR11um3Rescaled
from AVHRR.AVHRR_12um_3_rescaled_sensor_data import AVHRR12um3Rescaled
from AVHRR.AVHRR_37um_2_rescaled_sensor_data import AVHRR37um2Rescaled
from AVHRR.AVHRR_11um_3_sensor_data import AVHRR11um3
from AVHRR.AVHRR_12um_3_sensor_data import AVHRR12um3
from AVHRR.AVHRR_37um_2_sensor_data import AVHRR37um2
from AVHRR.AVHRR_11um_4_rescaled_sensor_data import AVHRR11um4Rescaled
from AVHRR.AVHRR_12um_4_rescaled_sensor_data import AVHRR12um4Rescaled
from AVHRR.AVHRR_37um_3_rescaled_sensor_data import AVHRR37um3Rescaled
from AVHRR.AVHRR_11um_4_sensor_data import AVHRR11um4
from AVHRR.AVHRR_12um_4_sensor_data import AVHRR12um4
from AVHRR.AVHRR_37um_3_sensor_data import AVHRR37um3
from AVHRR.AVHRR_sim_3_sensor_data import AVHRRSim3
from AVHRR.AVHRR_sim_4_sensor_data import AVHRRSim4
from HIRS.HIRS__Ch01_3_sensor_data import HIRSCh01_3
from HIRS.HIRS__Ch02_3_sensor_data import HIRSCh02_3
from HIRS.HIRS__Ch03_3_sensor_data import HIRSCh03_3
from HIRS.HIRS__Ch04_3_sensor_data import HIRSCh04_3
from HIRS.HIRS__Ch05_3_sensor_data import HIRSCh05_3
from HIRS.HIRS__Ch06_3_sensor_data import HIRSCh06_3
from HIRS.HIRS__Ch07_3_sensor_data import HIRSCh07_3
from HIRS.HIRS__Ch08_3_sensor_data import HIRSCh08_3
from HIRS.HIRS__Ch09_3_sensor_data import HIRSCh09_3
from HIRS.HIRS__Ch10_3_sensor_data import HIRSCh10_3
from HIRS.HIRS__Ch11_3_sensor_data import HIRSCh11_3
from HIRS.HIRS__Ch12_3_sensor_data import HIRSCh12_3
from HIRS.HIRS__Ch13_3_sensor_data import HIRSCh13_3
from HIRS.HIRS__Ch14_3_sensor_data import HIRSCh14_3
from HIRS.HIRS__Ch15_3_sensor_data import HIRSCh15_3
from HIRS.HIRS__Ch16_3_sensor_data import HIRSCh16_3
from HIRS.HIRS__Ch17_3_sensor_data import HIRSCh17_3
from HIRS.HIRS__Ch18_3_sensor_data import HIRSCh18_3
from HIRS.HIRS__Ch19_3_sensor_data import HIRSCh19_3
from MW.MW____Ch03_3_sensor_data import MWCh03_3


'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "22/10/2018"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class SensorDataFactory:
    def __init__(self):
        self.sensors_data = OrderedDict([(AVHRR11um3Rescaled().name, AVHRR11um3Rescaled()),
                                         (AVHRR12um3Rescaled().name, AVHRR12um3Rescaled()),
                                         (AVHRR37um2Rescaled().name, AVHRR37um2Rescaled()),
                                         (AVHRR11um3().name, AVHRR11um3()),
                                         (AVHRR12um3().name, AVHRR12um3()),
                                         (AVHRR37um2().name, AVHRR37um2()),
                                         (AVHRR11um4Rescaled().name, AVHRR11um4Rescaled()),
                                         (AVHRR12um4Rescaled().name, AVHRR12um4Rescaled()),
                                         (AVHRR37um3Rescaled().name, AVHRR37um3Rescaled()),
                                         (AVHRR11um4().name, AVHRR11um4()),
                                         (AVHRR12um4().name, AVHRR12um4()),
                                         (AVHRR37um3().name, AVHRR37um3()),
                                         (AVHRRSim3().name, AVHRRSim3()),
                                         (AVHRRSim4().name, AVHRRSim4()),
                                         (HIRSCh01_3().name, HIRSCh01_3()),
                                         (HIRSCh02_3().name, HIRSCh02_3()),
                                         (HIRSCh03_3().name, HIRSCh03_3()),
                                         (HIRSCh04_3().name, HIRSCh04_3()),
                                         (HIRSCh05_3().name, HIRSCh05_3()),
                                         (HIRSCh06_3().name, HIRSCh06_3()),
                                         (HIRSCh07_3().name, HIRSCh07_3()),
                                         (HIRSCh08_3().name, HIRSCh08_3()),
                                         (HIRSCh09_3().name, HIRSCh09_3()),
                                         (HIRSCh10_3().name, HIRSCh10_3()),
                                         (HIRSCh11_3().name, HIRSCh11_3()),
                                         (HIRSCh12_3().name, HIRSCh12_3()),
                                         (HIRSCh13_3().name, HIRSCh13_3()),
                                         (HIRSCh14_3().name, HIRSCh14_3()),
                                         (HIRSCh15_3().name, HIRSCh15_3()),
                                         (HIRSCh16_3().name, HIRSCh16_3()),
                                         (HIRSCh17_3().name, HIRSCh17_3()),
                                         (HIRSCh18_3().name, HIRSCh18_3()),
                                         (HIRSCh19_3().name, HIRSCh19_3()),
                                         (MWCh03_3().name, MWCh03_3())])

    def get_names(self):
        return self.sensors_data.keys()

    def get_sensor_data(self, name):
        return self.sensors_data[name].get_sensor_data()


if __name__ == "__main__":
    pass
