"""
Calc Wellpath

Author: Scott K
Created: 2021-05-23

"""

import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# well data inputs
file_name = "Final Field_McClintic A 9H Surveys trimmed.csv"
tgt_azi = 345.44
shl_x, shl_y = 1507302.1, 816452.02
kb = 2709.5

#import survey data
with open(file_name) as file:
	remove_headers = file.readline()[1:]
	
	print("\n*********************")
	print(f"Opened file: {file_name}")
	
	md_list = []
	inc_list = []
	azi_list = []
	
	for i in file:
		temp_list = i.split(",")
		temp_md = float(temp_list[0])
		temp_inc = float(temp_list[1])
		temp_azi = float(temp_list[2])
		
		md_list.append(temp_md)
		inc_list.append(temp_inc)
		azi_list.append(temp_azi)

print("Read file contents")

file.close()
print(f"Closed file {file_name}")
print("*********************\n")

surveys = [md_list, inc_list, azi_list]

# func that takes md, inc, azi inputs and returns a full survey sheet using minimum curvature
def min_curv(surveys, shl_x, shl_y, tgt_azi):
	#unpack survey data
	md_list, inc_list, azi_list = surveys
	#define additional lists
	dls_list = [0]
	rf_list = [1]
	tvd_list = [md_list[0]]
	easting_list = [0]
	northing_list = [0]
	x_list = [shl_x]
	y_list = [shl_y]
	vs_list = [0]
	cdist_list = [0]
	cdir_list = [0]
	
	#initialize running sum variables prior to starting calcs in for-loop
	runsum_tvd = tvd_list[0]
	runsum_easting = 0
	runsum_northing = 0
	runsum_vs = 0
	
	# for-loop calcs all values in standard survey sheet
	for svy in range(1, len(md_list)):
		# calc svy to svy course length
		course_length = md_list[svy] - md_list[svy-1]
		
		# convert inc and azi values to radians
		inc2_rad = math.radians(inc_list[svy])
		azi2_rad = math.radians(azi_list[svy])
		inc1_rad = math.radians(inc_list[svy-1])
		azi1_rad = math.radians(azi_list[svy-1])
		
		# calc DLS and append to list
		dog_leg = math.acos(((math.sin(inc1_rad) * math.sin(inc2_rad) * math.cos(azi2_rad - azi1_rad)) + (math.cos(inc1_rad) * math.cos(inc2_rad))))
		dog_leg_degrees = math.degrees(dog_leg)
		dls = (dog_leg_degrees * 100) / course_length
		dls_list.append(dls)
		
		# calc RF and append to list
		if dog_leg == 0:
			rf = 1
		else:
			rf = math.tan(dog_leg / 2) * (2 / dog_leg)
		rf_list.append(rf)
		
		# calc TVD and append to list
		dtvd = (course_length / 2) * (math.cos(inc1_rad) + math.cos(inc2_rad)) * rf
		runsum_tvd += dtvd
		tvd_list.append(runsum_tvd)
		
		# calc Easting and Northing and append to lists
		deasting = (course_length / 2) * ((math.sin(inc1_rad) * math.sin(azi1_rad)) + (math.sin(inc2_rad) * math.sin(azi2_rad))) * rf
		dnorthing = (course_length / 2) * ((math.sin(inc1_rad) * math.cos(azi1_rad)) + (math.sin(inc2_rad) * math.cos(azi2_rad))) * rf
		runsum_easting += deasting
		runsum_northing += dnorthing
		easting_list.append(runsum_easting)
		northing_list.append(runsum_northing)
		x = shl_x + runsum_easting
		y = shl_y + runsum_northing
		
		# calc closure distance, closure direction, vertical section and append to lists
		cdist = math.hypot(runsum_northing, runsum_easting)
		cdir = math.atan(runsum_easting / runsum_northing)
		
		dir_diff = math.radians(tgt_azi) - cdir
		if math.degrees(dir_diff) >= 360:
			dir_diff -= math.radians(360)
		vs = cdist * math.cos(dir_diff)
		
		if tgt_azi > 270:
			min_azi = tgt_azi + 90 - 360
			max_azi = tgt_azi - 90
		elif tgt_azi < 90:
			max_azi = tgt_azi - 90 + 360
			min_azi = tgt_azi + 90
		else:
			max_azi = tgt_azi + 90
			min_azi = tgt_azi - 90
		
		vs1 = abs(vs_list[svy-1])
		vs2 = abs(vs)
		dvs = abs(vs2 - vs1)
		if azi_list[svy] > min_azi and azi_list[svy] < max_azi:
			dvs *= -1
		runsum_vs += dvs
		vs_list.append(runsum_vs)
		print(md_list[svy], runsum_tvd)
	return easting_list, northing_list, tvd_list, vs_list

def plot_wellbore_xyz(coordinates):
	# unpack xyz coordinates
	x, y, z = coordinates
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	fig.show()
	return

def plot_wellbore_tvd_vs(coordinates):
	# unpack tvd and vs coordinates
	tvd, vs = coordinates
	
	fig, ax = plt.subplots(1)
	plt.xlabel('Vertical Section')
	plt.ylabel('TVD')
	ax.plot(vs, tvd, lw=3, color='black', label='Wellbore')
	
	max_tvd = max(tvd)
	
	max_vs = max(vs)
	min_vs = min(vs)
	
	ax.set_ylim(max_tvd -50, max_tvd + 40)
	ax.set_xlim(min_vs - 500, max_vs + 500)
	ax.invert_yaxis()
	
	plt.legend()
	plt.show()

easting_list, northing_list, tvd_list, vs_list = min_curv(surveys, shl_x, shl_y, tgt_azi)

xyz = easting_list, northing_list, tvd_list

tvd_vs = [tvd_list, vs_list]

plot_wellbore_tvd_vs(tvd_vs)

plot_wellbore_xyz(xyz)