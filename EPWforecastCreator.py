# No way to get atomospheric pressure forecasts
# need to Add rain and weather flaggs

# Needs to be checked for accuracy



#all these imports are standard on most modern python implementations
import math
#import library to do http requests:
import urllib2
#import xml parser called minidom:
from xml.dom.minidom import parseString


pi = math.pi
sin = math.sin
cos = math.cos
acos = math.acos
tan = math.tan
asin = math.asin
radians=math.radians
degrees=math.degrees
exp = math.exp
log=math.log

# Convert radians to degrees
rtd = 180/pi

# Convert degrees to radians
dtr = pi/180


# User input of zipcode
zipcode = str(raw_input("Enter Zipcode of Location: "))
# IF a different zipcode is needed change it below or uncomment the line above and comment the line below.
#zipcode = '60601'


def hour_angle(apparent_solar_time):
	return (15 * (apparent_solar_time(hour_of_day, longitude, equation_of_time, local_standard_time_meridian, local_standard_time)))-180

def apparent_solar_time(hour_of_day, longitude, equation_of_time, local_standard_time_meridian, local_standard_time):
	return local_standard_time(hour_of_day)+(equation_of_time(day_number)/60)+((local_standard_time_meridian(longitude)-longitude)/15)

def local_standard_time(hour_of_day):
	import time
	#this will only work if the program is run in the U.S. 
	localtime = time.localtime(time.time())
	daylight_savings_is_true=localtime.tm_isdst
	return hour_of_day-daylight_savings_is_true
	
def local_standard_time_meridian(longitude):
	import time
	#this will nly work if the program is run in the U.S. 
	localtime = time.localtime(time.time())
	daylight_savings_is_true=localtime.tm_isdst
	time2=time_stamp[0]
	return (int(time2[19:22])-daylight_savings_is_true)*15
		
def equation_of_time(day_number):
	if 1 <= day_number <=106:
		a=-14.2
		b=7
		c=111
	elif 107 <= day_number <=166:
		a=4
		b=106
		c=59
	elif 167 <= day_number <=246:
		a=-6.5
		b=166
		c=80
	elif 247 <= day_number <=365:
		a=16.4
		b=247
		c=113
	#return a*sin(pi*(day_number-b)/c)
	return 229.18*((0.0418*sin((4*pi*(day_number-4)/365.24)+3.5884)-0.0334*sin(2*pi*(day_number-4)/365.24)))
	
def solar_declination(day_number):
	return 23.45 * sin( 2*pi/365 * (day_number + 284) )

def time_of_sunrise(latitude, solar_declination):
	return 12/pi * acos(tan(radians(latitude)) * tan(radians(solar_declination(day_number))))
	
def time_of_sunset(latitude, solar_declination):
	return 12/pi * (2*pi - acos(tan(radians(latitude)) * tan(radians(solar_declination(day_number)))))

def possible_sunshine_hours(latitude, solar_declination):
	return time_of_sunset(latitude, solar_declination)-time_of_sunrise(latitude, solar_declination)

def solar_altitude(
	solar_declination, 
	latitude, 
	hour_angle):
	sa= rtd*asin((sin(solar_declination(day_number)*dtr) * sin(latitude*dtr)) + (cos(solar_declination(day_number)*dtr) \
	* cos(latitude*dtr) * cos(hour_angle(apparent_solar_time)*dtr)))
	if sa<0:
		return 0
	else:
		return sa
		
# check for noon and sign before and after noon
def solar_azimuth(
	solar_altitude,
	latitude, 
	solar_declination
	):
	if apparent_solar_time(hour_of_day, longitude, equation_of_time, local_standard_time_meridian, local_standard_time) == 12:
		return 0
	elif apparent_solar_time(hour_of_day, longitude, equation_of_time, local_standard_time_meridian, local_standard_time) < 12:
		return -rtd*acos(((sin(solar_altitude(solar_declination, latitude, hour_angle)*dtr) * sin(latitude*dtr)) \
		- sin(solar_declination(day_number)*dtr)) / (cos(solar_altitude(solar_declination, latitude, hour_angle)*dtr) * cos(latitude*dtr)))	
	else:
		return rtd*acos(((sin(solar_altitude(solar_declination, latitude, hour_angle)*dtr) * sin(latitude*dtr)) \
		- sin(solar_declination(day_number)*dtr)) / (cos(solar_altitude(solar_declination, latitude, hour_angle)*dtr) * cos(latitude*dtr)))	

# apparent solar irradiation at air mass equals zero(W/m2)

def apparent_solar_irradiation_air_mass_zero(
	day_number
	):
	return 1148 + 57 * cos(day_number*dtr)
	
#atmospheric extinction coefficient

def atmospheric_extinction_coefficient(
	day_number
	):
	return 0.161 - 0.0225 * cos(day_number*dtr)
	
# Estimated Direct Solar Irradiance (W/m2) 

def estimated_direct_solar_irradiance(
	apparent_solar_irradiation_air_mass_zero, 
	atmospheric_extinction_coefficient, 
	solar_altitude
	):
	asiamz = apparent_solar_irradiation_air_mass_zero(
	day_number
	)
	aec = atmospheric_extinction_coefficient(
	day_number
	)
	sa = solar_altitude(
	solar_declination, 
	latitude, 
	hour_angle
	)
	#Prevent division by zero
	if sin(dtr*sa)< .0001:
		return 0
	else:
		return asiamz/exp(aec/sin(dtr*sa))

def horizontal_IR(temp_now, temp_dew, cloud_cover):
	sigma=float(5.6697e-8)
	tempK=float(temp_now+273)
	temp_dewK=float(temp_dew+273)
	sky_emissivity=(.787+.764*log((temp_dewK/273)))*(1+(.0224*cloud_cover)-(.0035*cloud_cover**2)+(.00028*cloud_cover**3))
	return sky_emissivity*sigma*(tempK**4)


def zhang_huang_Est_hourly_Sol_Rad (solar_altitude, cloud_cover, percent_humidity, temp_now, temp_minus_3h, wind_speed):
	global_solar_const= 1355
	c0=.5598
	c1=.4982
	c2=-.6762
	c3=.02842
	c4=-.00317
	c5=.014
	d=-17.853
	k=.843
	sol_alt=math.sin(math.radians((solar_altitude(solar_declination, latitude, hour_angle))))
	j1=global_solar_const*sol_alt
	j2=c0+(c1*cloud_cover)+(c2*math.pow(cloud_cover,2))+(c3*(temp_now-temp_minus_3h))+(c4*percent_humidity)+(c5*wind_speed)
	return ((j1*j2)+d)/k
			
#This is the Direct Insolation Solar Code (DISC) model developed by Dr. E. Maxwell of the National Renewable Energy Laboratory. The model development, #description, and validation is found in "A Quasi-Physical Model for Converting Hourly Global Insolation to Direct Normal Insolation" NREL TR/215-3087, #August, 1987. Solar Energy Research Institute (National Renewable Energy Laboratory ) 1617 Cole Blvd.
#The model converts HOURLY Global Horizontal Data, inserted by the user in column N,
#into direct normal insolation (column U) based upon a quasi-physical relationship between the global clearness index (Kt) and the direct normal #clearness index (Kn).
#Mean Bias errors in derived results are on the order of -50 W/m2 and RMS errrors on the estimated DNI are on the order of +/- 150 W/m2.
#D.R. Myers: National Renewable Energy Laboratory, daryl_myers@nrel.goc



#A15
station_pressure=1000


#A17 Pressure Correction for Air Mass
pressure_correction_air_mass=station_pressure/1013.25


#H4 ETR ETR: Extraterrestrial Solar Irradiance based on solar constant of 1367 W/m2 and Spencer radius vector algorithm (Search #2, 1971):
def extraterrestrial_solar_irradiance(day_number):
	day_angle=6.283185*(day_number-1)/365
	return 1367*(1.00011+0.034221*cos(day_angle)+0.00128*sin(day_angle)+0.000719*cos(2*day_angle)+0.000077*sin(2*day_angle))


#zenith_angle Zenith Angle Zenith angle fo Sun, complement of solar elevation
def zenith_angle(solar_altitude):
	return 90-solar_altitude(solar_declination, latitude, hour_angle)

#M4 AM Pressure adj. Air Mass = [1/cos(zenith)]*[P/Po]
#the relative path length through the atmosphere. AM 1.0=> sun overhead

def pressure_adj_air_mass(zenith_angle, pressure_correction_air_mass):
	if zenith_angle(solar_altitude)<80:
		return ((1/(cos(radians(zenith_angle(solar_altitude)))+0.15/(93.885-zenith_angle(solar_altitude))**1.253))*pressure_correction_air_mass)
	else:
		return 0

#K1_GCI K1 Global Clearness Index. Meas Glob/[ETR*cos(z)]
# K1_GCI=IF(pressure_adj_air_mass>0,GHI/(cos(radians(zenith_angle))*extraterrestrial_solar_irradiance),0)
def K1_GCI(zenith_angle, extraterrestrial_solar_irradiance):
	if pressure_adj_air_mass(zenith_angle, pressure_correction_air_mass)>0:
		return zhang_huang_Est_hourly_Sol_Rad (solar_altitude, cloud_cover, percent_humidity, temp_now, temp_minus_3h, wind_speed)/(cos(radians(zenith_angle(solar_altitude)))*extraterrestrial_solar_irradiance(day_number))
	else:
		return 0

#P4 A
#A=IF(K1_GCI>0,IF(K1_GCI>0.6, -5.743+21.77*K1_GCI-27.49*K1_GCI**2+11.56*K1_GCI**3, IF(K1_GCI<0.6,0.512-1.56*K1_GCI+2.286*K1_GCI**2-2.222*K1_GCI**3)),0)		
def A(K1_GCI):
	if K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)>0:
		if K1_GCI(zenith_angle, extraterrestrial_solar_irradiance) > 0.6:
			return -5.743+21.77*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)-27.49*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)**2+11.56*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)**3
		else:
			return 0.512-1.56*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)+2.286*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)**2-2.222*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)**3
	else:
		return 0

#Q4 B
#B=IF(K1_GCI>0,IF(K1_GCI>0.6, 41.4-118.5*K1_GCI+66.05*K1_GCI**2+31.9*K1_GCI**3, IF(K1_GCI<0.6,0.37+0.962*K1_GCI,0)),0)
def B(K1_GCI):
	if K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)>0:
		if K1_GCI(zenith_angle, extraterrestrial_solar_irradiance) > 0.6:
			return 41.4-118.5*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)+66.05*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)**2+31.9*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)**3
		else:
			return 0.37+0.962*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)
	else:
		return 0

#R4 C
#C=IF(K1_GCI>0,IF(K1_GCI>0.6, -47.01+184.2*K1_GCI-222*K1_GCI**2+73.81*K1_GCI**3, IF(K1_GCI<0.6,-0.28+0.932*K1_GCI-2.048*K1_GCI**2,0)),0)
def C(K1_GCI):
	if K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)>0:
		if K1_GCI(zenith_angle, extraterrestrial_solar_irradiance) > 0.6:
			return -47.01+184.2*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)-222*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)**2+73.81*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)**3
		else:
			return -0.28+0.932*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)-2.048*K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)**2
	else:
		return 0
#S4 Delta Kn
#Delta_Kn=IF(K1_GCI>0,(A+B*EXP(C*pressure_adj_air_mass)),0)
def Delta_Kn(K1_GCI):
	if K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)>0:
		return A(K1_GCI)+B(K1_GCI)*exp(C(K1_GCI)*pressure_adj_air_mass(zenith_angle, pressure_correction_air_mass))
	else:
		return 0

#T4 Knc Computed Direct Beam Clearness Index DNI*/ETRDNI
#Knc=IF(K1_GCI>0,(0.886-0.122*pressure_adj_air_mass+0.0121*(pressure_adj_air_mass)**2-0.000653*(pressure_adj_air_mass)**3+0.000014*(pressure_adj_air_mass)**4),0)
def Knc(pressure_adj_air_mass):
	if K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)>0:
		return 0.886-0.122*pressure_adj_air_mass(zenith_angle, pressure_correction_air_mass)+0.0121*(pressure_adj_air_mass(zenith_angle, pressure_correction_air_mass))**2-0.000653*(pressure_adj_air_mass(zenith_angle, pressure_correction_air_mass))**3+0.000014*(pressure_adj_air_mass(zenith_angle, pressure_correction_air_mass))**4
	else:
		return 0
#U4 Results, Estimated DNI, Computed Direct Normal: ETR(DNI)*Knc
#estimated_DNI=IF(K1_GCI>0,extraterrestrial_solar_irradiance*(Knc-Delta_Kn),0)
def estimated_DNI_DISK_method(extraterrestrial_solar_irradiance, Knc, Delta_Kn):
	if K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)>0:
		return extraterrestrial_solar_irradiance(day_number)*(Knc(pressure_adj_air_mass)-Delta_Kn(K1_GCI))
	else:
		return 0

def diffuse_horizontal_radiation(zhang_huang_Est_hourly_Sol_Rad, estimated_DNI_DISK_method, zenith_angle):
	if zhang_huang_Est_hourly_Sol_Rad (solar_altitude, cloud_cover, percent_humidity, temp_now, temp_minus_3h, wind_speed)>0:
		return zhang_huang_Est_hourly_Sol_Rad (solar_altitude, cloud_cover, percent_humidity, temp_now, temp_minus_3h, wind_speed)-estimated_DNI_DISK_method(extraterrestrial_solar_irradiance, Knc, Delta_Kn)*(sin(radians(90-zenith_angle(solar_altitude))))
	else:
		return 0
###############################################################################################################################################################


 
# Routine to download and parse the weather forecast from NOAA from a zipcode.

#Get Date
from datetime import *
print "time of forecast"
time_of_forecast=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print time_of_forecast
startdate=str(date.today()+timedelta(days=1))
enddate=str(date.today()+timedelta(days=7))
print startdate
print enddate
print "zipcode"
print zipcode





siteurl='http://graphical.weather.gov/xml/sample_products/browser_interface/ndfdXMLclient.php?whichClient=NDFDgenMultiZipCode&lat=&lon=&listLatLon=&lat1=&lon1=&lat2=&lon2=&resolutionSub=&listLat1=&listLon1=&listLat2=&listLon2=&resolutionList=&endPoint1Lat=&endPoint1Lon=&endPoint2Lat=&endPoint2Lon=&listEndPoint1Lat=&listEndPoint1Lon=&listEndPoint2Lat=&listEndPoint2Lon=&zipCodeList='+zipcode+'&listZipCodeList=&centerPointLat=&centerPointLon=&distanceLat=&distanceLon=&resolutionSquare=&listCenterPointLat=&listCenterPointLon=&listDistanceLat=&listDistanceLon=&listResolutionSquare=&citiesLevel=&listCitiesLevel=&sector=&gmlListLatLon=&featureType=&requestedTime=&startTime=&endTime=&compType=&propertyName=&product=time-series&begin='+startdate+'T00%3A00%3A00&end='+enddate+'T00%3A00%3A00&Unit=m&temp=temp&qpf=qpf&pop12=pop12&dew=dew&wspd=wspd&wdir=wdir&sky=sky&rh=rh&Submit=Submit'
#download the file:
file = urllib2.urlopen(siteurl)



#file="/desktop/data2.xml"
#convert to string:
#print siteurl

#from xml.dom.minidom import parseString
#file= open('/Users/troypeters/Desktop/data2.xml', )
data=file.read()
destinationfile=zipcode+time_of_forecast+'.xml'
raw_xml = open (destinationfile, 'w') ## a will append, w will over-write 
#here we are writing the source file content to destination file
raw_xml.write(data)
#providing information that the task is completed
raw_xml.close()

file.close()
dom = parseString(data)

parameter=dom.getElementsByTagName('point')
for node in parameter:
	latitude_name=node.getAttribute('latitude')
	print "Latitude"
	print latitude_name
	longitude_name=node.getAttribute('longitude')
	print "Longitude"
	print longitude_name

		
time_stamp=[]
time_layout=[]		
parameter=dom.getElementsByTagName('time-layout')
for node in parameter:
	type_name=node.getAttribute('time-coordinate')
	print type_name
	units_name=node.getAttribute('summarization')
	print units_name
	time_layout_name=node.getAttribute('time-layout')
	print time_layout_name
	blist=node.getElementsByTagName('layout-key')
	for b in blist:
		par_name= b.childNodes[0].nodeValue
		time_layout.append(par_name)
		print 'look here'
		print par_name
	alist=node.getElementsByTagName('start-valid-time')
	for a in alist:
		value_text= a.childNodes[0].nodeValue
		time_stamp.append(value_text)
		#print value_text


#know when to sart tabulating data

step_subtraction_temp=time_layout[0]
step_subtraction=int(step_subtraction_temp[8:10])
#step_subtraction_tempA=time_layout[2]
#step_subtractionA=int(step_subtraction_tempA[8:10])
print step_subtraction
#print step_subtractionA


hourly_temperature_list=[]
wind_speed_list=[]
wind_direction_list=[]
cloud_cover_list=[]
percent_humidity_list=[]
		
parameter=dom.getElementsByTagName('temperature')
for node in parameter:
	type_name=node.getAttribute('type')
	print type_name
	units_name=node.getAttribute('units')
	print units_name
	time_layout_name=node.getAttribute('time-layout')
	print time_layout_name
	blist=node.getElementsByTagName('name')
	for b in blist:
		par_name= b.childNodes[0].nodeValue
		#print par_name
	alist=node.getElementsByTagName('value')
	for a in alist:
		value_text= a.childNodes[0].nodeValue
		hourly_temperature_list.append(value_text)
		#print value_text
parameter=dom.getElementsByTagName('wind-speed')
for node in parameter:
	type_name=node.getAttribute('type')
	print type_name
	units_name=node.getAttribute('units')
	print units_name
	time_layout_name=node.getAttribute('time-layout')
	print time_layout_name
	blist=node.getElementsByTagName('name')
	for b in blist:
		par_name= b.childNodes[0].nodeValue
		#print par_name
	alist=node.getElementsByTagName('value')
	for a in alist:
		value_text= a.childNodes[0].nodeValue
		wind_speed_list.append(value_text)
		#print value_text
parameter=dom.getElementsByTagName('direction')
for node in parameter:
	type_name=node.getAttribute('type')
	print type_name
	units_name=node.getAttribute('units')
	print units_name
	time_layout_name=node.getAttribute('time-layout')
	print time_layout_name
	blist=node.getElementsByTagName('name')
	for b in blist:
		par_name= b.childNodes[0].nodeValue
		#print par_name
	alist=node.getElementsByTagName('value')
	for a in alist:
		value_text= a.childNodes[0].nodeValue
		wind_direction_list.append(value_text)
		#print value_text
parameter=dom.getElementsByTagName('cloud-amount')
for node in parameter:
	type_name=node.getAttribute('type')
	print type_name
	units_name=node.getAttribute('units')
	print units_name
	time_layout_name=node.getAttribute('time-layout')
	print time_layout_name
	blist=node.getElementsByTagName('name')
	for b in blist:
		par_name= b.childNodes[0].nodeValue
		#print par_name
	alist=node.getElementsByTagName('value')
	for a in alist:
		value_text= a.childNodes[0].nodeValue
		cloud_cover_list.append(value_text)
		#print value_text
parameter=dom.getElementsByTagName('humidity')
for node in parameter:
	type_name=node.getAttribute('type')
	print type_name
	units_name=node.getAttribute('units')
	print units_name
	time_layout_name=node.getAttribute('time-layout')
	print time_layout_name
	blist=node.getElementsByTagName('name')
	for b in blist:
		par_name= b.childNodes[0].nodeValue
		#print par_name
	alist=node.getElementsByTagName('value')
	for a in alist:
		value_text= a.childNodes[0].nodeValue
		percent_humidity_list.append(value_text)
		#print value_text

list_length=len(wind_speed_list)
date_length=len(time_stamp)
day_of_year=[]
year_stamp=[]
month_of_year=[]
day_of_month=[]
hour_of_day_list=[]
for i in range(0,list_length):
	datter=time_stamp[i+step_subtraction]
	datter1=datter[:19]
	datter2=int(datter[19:22])*15
	print datter
	print datter1
	print datter2
	x= datetime.strptime(datter1,'%Y-%m-%dT%H:%M:%S')
	day_of_year.append(x.timetuple().tm_yday)
	year_stamp.append(x.year)
	month_of_year.append(x.month)
	day_of_month.append(x.day)
	hour_of_day_list.append(x.hour)
	#print x, day_of_year, month_of_year, day_of_month, hour_of_day
				
# User input of latitude
#latitude = float(input("Enter Latitude of Location: "))
latitude=float(latitude_name)
longitude=float(longitude_name)
#print 'test latitude'
#print latitude




list_length=len(wind_speed_list)		
print "temperature, dewpoint, windspeed, winddirection, cloudcover, relativehumidity"
for i in range(0,list_length):
	hour_of_day=hour_of_day_list[i]
	day_number=day_of_year[i]
	print "hour_of_day, hour_angle, apparent_solar_time, local_standard_time, local_standard_time_meridian, equation_of_time, solar_declination"
	print hour_of_day,hour_angle(apparent_solar_time), apparent_solar_time(hour_of_day, longitude, equation_of_time, local_standard_time_meridian, local_standard_time), local_standard_time(hour_of_day), local_standard_time_meridian(longitude), equation_of_time(day_number), solar_declination(day_number)

row_number=[]
for i in range(0,list_length):
	rn= 24*(day_of_year[i]-day_of_year[0])+hour_of_day_list[i]
	row_number.append(rn)
	#print row_number[i]

def interpolateGap(ts0, v0, ts1, v1):
	count = (ts1 - ts0) / timestampDistance
	return [ (ts0 + i * timestampDistance, v0 + (v1 - v0) * i / count)
		for i in range(1, count) ]

#Interpolate Cloud Cover
data_cloud_cover=[]
for i in range(0,list_length):
	data_cloud_cover.append((row_number[i], float(cloud_cover_list[i])))
timestampDistance = 1

def fillGap(data_cloud_cover, pos, ts0, v0, ts1, v1):
	data_cloud_cover[pos+1:pos] = interpolateGap(ts0, v0, ts1, v1)

for i in range(len(data_cloud_cover)-1, 0, -1):
	timestamp, value = data_cloud_cover[i]
	previousTimestamp, previousValue = data_cloud_cover[i-1]
	if previousTimestamp + timestampDistance < timestamp:
		fillGap(data_cloud_cover, i-1, previousTimestamp, previousValue, timestamp, value)

data_cloud_coverf = [x[1] for x in data_cloud_cover] #x[column that you do not want to remove]
print data_cloud_coverf

#Interpolate Percent Humidity
data_percent_humidity=[]
for i in range(0,list_length):
	data_percent_humidity.append((row_number[i], float(percent_humidity_list[i])))
timestampDistance = 1

def fillGap(data_percent_humidity, pos, ts0, v0, ts1, v1):
	data_percent_humidity[pos+1:pos] = interpolateGap(ts0, v0, ts1, v1)

for i in range(len(data_percent_humidity)-1, 0, -1):
	timestamp, value = data_percent_humidity[i]
	previousTimestamp, previousValue = data_percent_humidity[i-1]
	if previousTimestamp + timestampDistance < timestamp:
		fillGap(data_percent_humidity, i-1, previousTimestamp, previousValue, timestamp, value)

data_percent_humidityf = [x[1] for x in data_percent_humidity] #x[column that you do not want to remove]
print data_percent_humidityf

#Interpolate Wind Speed
data_wind_speed=[]
for i in range(0,list_length):
	data_wind_speed.append((row_number[i], float(wind_speed_list[i])))
timestampDistance = 1

def fillGap(data_wind_speed, pos, ts0, v0, ts1, v1):
	data_wind_speed[pos+1:pos] = interpolateGap(ts0, v0, ts1, v1)

for i in range(len(data_wind_speed)-1, 0, -1):
	timestamp, value = data_wind_speed[i]
	previousTimestamp, previousValue = data_wind_speed[i-1]
	if previousTimestamp + timestampDistance < timestamp:
		fillGap(data_wind_speed, i-1, previousTimestamp, previousValue, timestamp, value)

data_wind_speedf = [x[1] for x in data_wind_speed] #x[column that you do not want to remove]
print data_wind_speedf

#Interpolate Wind Direction
data_wind_direction=[]
for i in range(0,list_length):
	data_wind_direction.append((row_number[i], float(wind_direction_list[i])))
timestampDistance = 1

def fillGap(data_wind_direction, pos, ts0, v0, ts1, v1):
	data_wind_direction[pos+1:pos] = interpolateGap(ts0, v0, ts1, v1)

for i in range(len(data_wind_direction)-1, 0, -1):
	timestamp, value = data_wind_direction[i]
	previousTimestamp, previousValue = data_wind_direction[i-1]
	if previousTimestamp + timestampDistance < timestamp:
		fillGap(data_wind_direction, i-1, previousTimestamp, previousValue, timestamp, value)

data_wind_directionf = [x[1] for x in data_wind_direction] #x[column that you do not want to remove]
print data_wind_directionf

#Interpolate Air Temperature
data_hourly_temperature=[]
temp_half=len(hourly_temperature_list)/2
for i in range(0,temp_half):
	data_hourly_temperature.append((row_number[i], float(hourly_temperature_list[i])))
timestampDistance = 1

def fillGap(data_hourly_temperature, pos, ts0, v0, ts1, v1):
	data_hourly_temperature[pos+1:pos] = interpolateGap(ts0, v0, ts1, v1)

for i in range(len(data_hourly_temperature)-1, 0, -1):
	timestamp, value = data_hourly_temperature[i]
	previousTimestamp, previousValue = data_hourly_temperature[i-1]
	if previousTimestamp + timestampDistance < timestamp:
		fillGap(data_hourly_temperature, i-1, previousTimestamp, previousValue, timestamp, value)

data_hourly_temperaturef = [x[1] for x in data_hourly_temperature] #x[column that you do not want to remove]
print data_hourly_temperaturef

#Interpolate Dew Point Temperature
data_dew_point_temperature=[]
for i in range(0,temp_half):
	data_dew_point_temperature.append((row_number[i], float(hourly_temperature_list[i+temp_half])))
timestampDistance = 1

def fillGap(data_dew_point_temperature, pos, ts0, v0, ts1, v1):
	data_dew_point_temperature[pos+1:pos] = interpolateGap(ts0, v0, ts1, v1)

for i in range(len(data_dew_point_temperature)-1, 0, -1):
	timestamp, value = data_dew_point_temperature[i]
	previousTimestamp, previousValue = data_dew_point_temperature[i-1]
	if previousTimestamp + timestampDistance < timestamp:
		fillGap(data_dew_point_temperature, i-1, previousTimestamp, previousValue, timestamp, value)

data_dew_point_temperaturef = [x[1] for x in data_dew_point_temperature] #x[column that you do not want to remove]
print data_dew_point_temperaturef

count=int(hour_of_day_list[0])
day=int(day_of_year[0])
hour_of_dayf=[]
day_of_monthf=[]
yearf=[]
day_of_yearf=[]
month_of_yearf=[]
data_length=len(data_wind_directionf)
hour_of_dayf.append(hour_of_day_list[0])
for i in range(0,data_length):
	count=count +1
	f=datetime.fromordinal(day)
	hour_of_dayf.append(count)
	day_of_yearf.append(day)
	day_of_monthf.append(f.day)
	month_of_yearf.append(f.month)
	if hour_of_dayf[i]==23:
		count=0
	if hour_of_dayf[i]==24:
		day=day+1
	print hour_of_dayf[i], day_of_yearf[i], day_of_monthf[i], month_of_yearf[i]
	
# Writes a comma separated variable file
# of day_number, time_of_sunrise, time_of_sunset and possible_sunshine_hours
# also prints data to screen
file_name=zipcode+time_of_forecast+'.csv'
print file_name
with open (file_name, 'w') as f: 
	for i in range(0,data_length):
		hour_of_day=hour_of_dayf[i]
		day_number=day_of_yearf[i]
		if i<3:
			temp_minus_3h=float(data_hourly_temperaturef[i])
		else:
			temp_minus_3h=float(data_hourly_temperaturef[i-3])
		cloud_cover=float((data_cloud_coverf[i])/100)
		percent_humidity=float(data_percent_humidityf[i])
		temp_now=float(data_hourly_temperaturef[i])
		temp_dew=float(data_dew_point_temperaturef[i])
		wind_speed=float(data_wind_speedf[i])
		# Year
		n1=str(2014) #str(year_stamp[i])
		# Month
		n2=str(month_of_yearf[i])
		# Day
		n3=str(day_of_monthf[i])
		# Hour
		n4=str(hour_of_dayf[i])
		# Minute
		n5=str(0)
		# field data source
		a1= '?9?9?9?9E0?9?9?9?9?9?9?9?9?9?9?9?9?9?9?9*9*9*9*9*9'
		# Dry Bulb Temperature C
		n6=str(data_hourly_temperaturef[i])
		# Dew Point Temperature C
		n7=str(data_dew_point_temperaturef[i])
		# Relative Humidity
		n8=str(data_percent_humidityf[i])
		# Atmospheric Station Pressure missing= 999999.
		n9=str(999999.)
		# Extraterrestrial Horizontal Radiation missing= 9999.
		n10=str(9999)
		# Extraterrestrial Direct Normal Radiation missing= 9999.
		n11=str(9999)
		# Horizontal Infrared Radiation Intensity [Calculate]
		HIR=horizontal_IR(temp_now, temp_dew, cloud_cover)
		n12=str(HIR)
		# Global Horizontal Radiation
		print equation_of_time(day_number), longitude, local_standard_time_meridian(longitude), solar_altitude(solar_declination, latitude, hour_angle), apparent_solar_time(hour_of_day, longitude, equation_of_time, local_standard_time_meridian, local_standard_time), local_standard_time(hour_of_day)
		#print solar_altitude(solar_declination, latitude, hour_angle), cloud_cover, percent_humidity, temp_now, temp_minus_3h, wind_speed
		zhehsr=zhang_huang_Est_hourly_Sol_Rad(solar_altitude, cloud_cover, percent_humidity, temp_now, temp_minus_3h, wind_speed)
		if zhehsr<0:
			zhehsrA=0
		else:
			zhehsrA=zhehsr
		print 'zhehsrA', zhehsrA
		n13=str(zhehsrA)
		# Direct Normal Radiation
		DNI=estimated_DNI_DISK_method(extraterrestrial_solar_irradiance, Knc, Delta_Kn)
		#print 'zenith', zenith_angle(solar_altitude)
		#print 'extraterrestrial_solar_irradiance', extraterrestrial_solar_irradiance(day_number)
		#print "DNI", estimated_DNI_DISK_method(extraterrestrial_solar_irradiance, Knc, Delta_Kn)
		#print 'KNC', Knc(pressure_adj_air_mass)
		#print 'pressure_adj_air_mass', pressure_adj_air_mass(zenith_angle, pressure_correction_air_mass)
		#print 'pressure_correction_air_mass', pressure_correction_air_mass
		#print 'K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)', K1_GCI(zenith_angle, extraterrestrial_solar_irradiance)
		n14=str(DNI)
		# Diffuse Horizontal Radiation
		DIFF=diffuse_horizontal_radiation(zhang_huang_Est_hourly_Sol_Rad, estimated_DNI_DISK_method, zenith_angle)
		n15=str(DIFF)
		# Global Horizontal Illuminance
		n16=str(999999)
		# Direct Normal Illuminance
		n17=str(999999)
		# Diffuse Horizontal Illuminance
		n18=str(999999)
		# Zenith Luminance
		n19=str(9999)
		# Wind Direction
		n20=str(data_wind_directionf[i])
		# Wind Speed m/s
		n21=str(data_wind_speedf[i])
		# Total Sky Cover missing=99
		n22=str(99)
		# Opaque Sky Cover 0 to 10
		c_c=float(data_cloud_coverf[i]/10)
		n23=str(c_c)
		#Visibility
		n24=str(9999)
		# Ceiling Height
		n25=str(99999)
		# Present Weather Conditions
		n26=str(9)
		# Present Weather Codes
		n27=str(0)
		#
		n28=str(999)
		#
		n29=str(.999)
		#
		n30=str(999)
		#
		n31=str(99)
		#
		n32=str(999)
		#
		n33=str(999)
		#
		n34=str(99)


		f.write(n1)
		f.write(", ")
		f.write(n2)
		f.write(", ")
		f.write(n3)
		f.write(", ")
		f.write(n4)
		f.write(", ")
		f.write(n5)
		f.write(", ")
		f.write(a1)
		f.write(", ")
		f.write(n6)
		f.write(", ")
		f.write(n7)
		f.write(", ")
		f.write(n8)
		f.write(", ")
		f.write(n9)
		f.write(", ")
		f.write(n10)
		f.write(", ")
		f.write(n11)
		f.write(", ")
		f.write(n12)
		f.write(", ")
		f.write(n13)
		f.write(", ")
		f.write(n14)
		f.write(", ")
		f.write(n15)
		f.write(", ")
		f.write(n16)
		f.write(", ")
		f.write(n17)
		f.write(", ")
		f.write(n18)
		f.write(", ")
		f.write(n19)
		f.write(", ")
		f.write(n20)
		f.write(", ")
		f.write(n21)
		f.write(", ")
		f.write(n22)
		f.write(", ")
		f.write(n23)
		f.write(", ")
		f.write(n24)
		f.write(", ")
		f.write(n25)
		f.write(", ")
		f.write(n26)
		f.write(", ")
		f.write(n27)
		f.write(", ")
		f.write(n28)
		f.write(", ")
		f.write(n29)
		f.write(", ")
		f.write(n30)
		f.write(", ")
		f.write(n31)
		f.write(", ")
		f.write(n32)
		f.write(", ")
		f.write(n33)
		f.write(", ")
		f.write(n34)
		f.write("\n")