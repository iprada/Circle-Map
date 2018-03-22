from coverage import coverage

object = coverage('B06.bam','/isdata/kroghgrp/xsh723/projects/circle_map/test_data/circle_map_coverage/',mismatch=2)
object.find_circles()