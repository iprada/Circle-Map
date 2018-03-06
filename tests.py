import pybedtools as bt
import os

os.chdir("/isdata/kroghgrp/xsh723/projects/circle_map/test_data/realigner/test_command_interfance")

b = bt.BedTool('test_circles.bed')

a = bt.BedTool('test.bed')

a.cat(b).saveas('test.bed')


