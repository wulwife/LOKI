#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 16:23:38 2021

@author: Peidong SHI
@email: speedshi@hotmail.com
"""


import os
import obspy
from obspy import UTCDateTime


def vector2trace(datainfo, data, dir_output='./'):
    """
    Write a data vector to an obspy trace.
    
    Parameters
    ----------
    datainfo : dictionary
        contains information about the station and data, includes:
            datainfo['station_name']: str, the name of the station, required;
            datainfo['channel_name']: str, the channel name of the trace, required;
                                      NOTE len(channel_name) MUST <= 3;
            datainfo['dt']: time sampling interval of the data in second, required;
            datainfo['starttime']: datetime, the starting time of the trace, required;
            datainfo['network']: str, network name of the trace, optional;
            datainfo['location']: str,  location name of the trace, optional;
            
    data : numpy vector
        the data vector to be written, shape: npts*1.
    dir_output : str
        the directory for output file.

    Returns
    -------
    None.

    """
    
    
    trace = obspy.Trace()  # initilize an empty obspy trace
    
    # set the trace header information
    trace.stats.station = datainfo['station_name']
    trace.stats.channel = datainfo['channel_name']
    trace.stats.delta = datainfo['dt']
    trace.stats.starttime = UTCDateTime(datainfo['starttime'])
    if 'network' in datainfo:
        trace.stats.network = datainfo['network']
    if 'location' in datainfo:
        trace.stats.location = datainfo['location']
    
    # assign data to the trace
    trace.data = data
    
    # set the displayed datetime format in the output filename
    # NOTE here output until second
    timeformat = "%Y%m%dT%H%M%SZ"  
    
    # make sure output directory exist
    if not os.path.exists(dir_output):
        os.makedirs(dir_output)
    
    # output to file
    if trace.id[0] == '.':
        nametag = trace.id[1:] + '.' + trace.stats.starttime.datetime.strftime(timeformat) + '.mseed'
    else:
        nametag = trace.id + '.' + trace.stats.starttime.datetime.strftime(timeformat) + '.mseed'
    fname = os.path.join(dir_output, nametag)
    trace.write(fname, format="MSEED")
    
    del trace
    
    return


