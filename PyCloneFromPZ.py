import os
import warnings
from glob import glob
from json import load
from pathlib import Path

from LatLon import lat_lon as ll
from numpy import loadtxt, pi, zeros_like
from obspy import UTCDateTime as utc
from obspy import read_inventory
from obspy.io.xseed import Parser

warnings.filterwarnings("ignore")

"""
Script for making response file from pole-zero files

Inputs:
- All requires info like poles, zeros, digitizer sensitivity and normalization
  constant should be placed in a correct file name and instrument directory.
- Poles and zeros will be read in Hz.
- Sensor sensitivity will be read in V/m/s.
- Digitizer sensitivity will be read in V/count.
- Normalization constant will be read in Hz.
- File "request.json" will be read for loading information required by each station.

Outputs:
- For each network, station a directory will be created in "output" folder.

LogChanges:
2021-05-05 > Initiated.
2021-05-06 > Added two helper functions for creating "request.json" and "STATION0.HYP".
2022-08-25 > user will be asked for using helper function.
"""


def dlsv2json():
    """ Convert dataless to json file.
    """
    from obspy import read_inventory
    dlsv = input("\n+++ Dataless input file:\n")
    instrumentCode = input("+++ instrumentCode:\n")
    dlsv = read_inventory("BP.dlsv")
    with open("request.json", "w")as f:
        for net in dlsv:
            for sta in net:
                lon = sta.longitude
                lat = sta.latitude
                elv = sta.elevation
                sd = sta.start_date.strftime("%Y,%j")
                try:
                    ed = sta.end_date.strftime("%Y,%j")
                except:
                    ed = "2100,001"
                f.write("""{
"network" : "%s",
"station" : "%s",
"channel" : "BH",
"latitude" : %6.3f,
"longitude" : %6.3f,
"elevation" : %d,
"startTime" : "%s",
"endTime" : "%s",
"instrumentCode": "%s"
},\n""" % (net.code, sta.code, lat, lon, elv, sd, ed, instrumentCode))


def clearOutput():
    """Clear output directory 
    """
    if os.path.exists("output"):
        cmd = "rm -r output"
        os.system(cmd)


def dlsv2dic(dlsv, dic):
    """A helper function for making STATION0.HYP from dataless files

    Arguments:
        dlsv (dataless) file in dataless format
        dic (dict) a dictionary contains stations information

    Returns:
        dic (dict) a dictionary contains stations information
    """
    dlsv = read_inventory(dlsv)
    for net in dlsv:
        for sta in net:
            lat = ll.Latitude(sta.latitude)
            lon = ll.Longitude(sta.longitude)
            dic.update({sta.code: (sta.code, lat, lon, sta.elevation)})
    return dic


def dlsv2sta0():
    """Convert dataless to station0.hyp
    """
    dic = {}
    for net in glob(os.path.join("output", "DataLess", "*")):
        for sta in sorted(glob(os.path.join(net, "*"))):
            for dlsv in glob(os.path.join(sta, "*.dlsv")):
                dic = dlsv2dic(dlsv, dic)
    with open("STATION0.HYP", "w") as f:
        for i in dic:
            f.write("  %-4s%02d%05.2f%1s %02d%05.2f%1s%0004d\n" % (dic[i][0],
                                                                   dic[i][1].degree,
                                                                   dic[i][1].decimal_minute, 
                                                                   dic[i][1].get_hemisphere(
            ),
                dic[i][2].degree, dic[i][2].decimal_minute, dic[i][2].get_hemisphere(
            ),
                dic[i][3]))


def cloneFromPZ():
    """Start cloning from PZ files
    """
    with open("request.json") as f:
        requests = load(f)
    for request in requests["stations"]:
        net, sta, seismometer, lat, lon, elv, st, et, ic = request.values()
        # Results will be saved into
        saveDir = os.path.join("output", "DataLess", net, sta)
        Path(saveDir).mkdir(parents=True, exist_ok=True)
        # Working on
        print("Net=%2s;Sta=%4s;Seismometer=%2s;Lat=%6.3f;Lon=%6.3f;Elv=%0004d;Start=%s;End=%s;Inst=%s" % (
            net, sta, seismometer, lat, lon, elv, st, et, ic))
        # Read sensor sensitivity in V/m/s
        sensor_sensitivity = loadtxt(os.path.join(
            "instruments", ic, "sensor_sensitivity.dat"))
        # Read digitizer sensitivity in V/count which will be reversed to get count/V
        digitizer_sensitivity = 1 / \
            loadtxt(os.path.join("instruments", ic, "digitizer_sensitivity.dat"))
        # Read Normalization constant at 1Hz
        A0 = loadtxt(os.path.join("instruments", ic,
                     "normalization_constant.dat"))
        # Read Poles in Hz
        poles = loadtxt(os.path.join("instruments", ic, "poles.dat"))
        # Read Zeros in Hz
        zeros = loadtxt(os.path.join("instruments", ic, "zeros.dat"))
        # Convert A0, Poles and Zeros from Hz to Rad/s
        A0_rad_s = A0*(2*pi)**(poles.size/2 - zeros.size/2)
        poles_rad_s = poles*(2*pi)
        zeros_rad_s = zeros*(2*pi)
        # Read RESP template file for Z,N,E channels
        p = Parser(os.path.join("template", "generic.dlsv"))
        blk = p.blockettes
        mult = int(len(blk[58])/3)
        for i, cha in enumerate(["E", "N", "Z"]):
            blk[50][0].network_code = net
            blk[50][0].station_call_letters = sta
            blk[50][0].site_name = sta
            blk[50][0].latitude = lat
            blk[50][0].longitude = lon
            blk[50][0].elevation = elv
            blk[50][0].start_effective_date = utc(st)
            blk[50][0].end_effective_date = utc(et)
            blk[52][i].channel_identifier = "%s%s" % (seismometer, cha)
            blk[52][i].location_identifier = ""
            blk[52][i].latitude = lat
            blk[52][i].longitude = lon
            blk[52][i].elevation = elv
            blk[52][i].start_date = utc(st)
            blk[52][i].end_date = utc(et)
            blk[53][i].number_of_complex_poles = int(poles_rad_s.size / 2)
            blk[53][i].real_pole = poles_rad_s[:, 0].tolist()
            blk[53][i].imaginary_pole = poles_rad_s[:, 1].tolist()
            blk[53][i].real_pole_error = zeros_like(poles_rad_s[:, 0]).tolist()
            blk[53][i].imaginary_pole_error = zeros_like(
                poles_rad_s[:, 1]).tolist()
            blk[53][i].number_of_complex_zeros = int(zeros_rad_s.size / 2)
            blk[53][i].real_zero = zeros_rad_s[:, 0].tolist()
            blk[53][i].imaginary_zero = zeros_rad_s[:, 1].tolist()
            blk[53][i].real_zero_error = zeros_like(zeros_rad_s[:, 0]).tolist()
            blk[53][i].imaginary_zero_error = zeros_like(
                zeros_rad_s[:, 1]).tolist()
            blk[53][i].A0_normalization_factor = A0_rad_s
            blk[53][i].normalization_frequency = 1.0
            # stage sequence number 1, seismometer gain
            blk[58][i*mult].sensitivity_gain = sensor_sensitivity[i]
            # stage sequence number 2, digitizer gain
            blk[58][i*mult+1].sensitivity_gain = digitizer_sensitivity[i]
            # stage sequence number 0, overall sensitivity
            blk[58][(i+1)*mult-1].sensitivity_gain = sensor_sensitivity[i] * \
                digitizer_sensitivity[i]
        # Save in dlsv format
        outFile = os.path.join(saveDir, "%s.%s_%s.dlsv" %
                               (net, sta, st.replace(",", "_")))
        p.write_seed(outFile)
        # Uncomment the following line if you want to get RESP files
        RESPDir = os.path.join("output", "RESP")
        Path(RESPDir).mkdir(parents=True, exist_ok=True)
        os.system("obspy-dataless2resp %s" % (outFile))
        os.system("mv RESP* %s" % (RESPDir))
        # Uncomment the following line if you want to get SCML files
        SCMLDir = os.path.join("output", "SCML")
        Path(SCMLDir).mkdir(parents=True, exist_ok=True)
        os.system("dlsv2inv %s %s" % (outFile, "%s.%s_%s.xml" %
                                      (net, sta, st.replace(",", "_"))))
        os.system("mv *xml %s" % (SCMLDir))


# GO!
mode = input("\n+++ Mode?\n1-Start cloning from PZ files,\n2-Convert DataLess to STATION0.HYP,\n3-Convert DataLess to json file.\n")
if mode == "1":
    cloneFromPZ()
elif mode == "2":
    dlsv2sta0()
elif mode == "3":
    dlsv2json()
else:
    print("Wrong choice, Bye!")
