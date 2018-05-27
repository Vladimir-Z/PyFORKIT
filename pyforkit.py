#!/usr/bin/env python3
# coding: utf8

"""
This script was created for reading MigroMag files and plotting FORC
(first-order reversal curves) diagram as well as related figures. 
This script was inspired by FORKIT program which was written 
Gary Acton (http://paleomag.ucdavis.edu/software-forcit.html) 
and contains some quotes from it.
"""

from __future__ import print_function  # this line only for compability with Python 2

import argparse
import collections
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import UnivariateSpline

__author__ = "Vladimir Zverev (vladimir.zverev@urfu.ru), Alla Dobroserdova"


def create_parser():
    cmd_parser = argparse.ArgumentParser(description=("script reads file .frc and creates FORC"))
    cmd_parser.add_argument(
        '-i', '--filename', type=str, help='relative path to the *.frc file')
    cmd_parser.add_argument(
        '-f',
        '--only-forc',
        action='store_true',
        help='Plot only FORC diagram',
        default=False)
    cmd_parser.add_argument(
        '-sf',
        '--smoothing-factor',
        type=bandwidth_type,
        help='this factor is used for smoothing FORC diagram.',
        default=3)
    cmd_parser.add_argument(
        '-r',
        '--plot-raw-data',
        action='store_true',
        help='plot all curves directrly from MicroMag-file',
        default=False)
    return cmd_parser


def bandwidth_type(x):
    x = int(x)
    if 1 <= x <= 5:
        raise argparse.ArgumentTypeError(
            "Smoothing factor should lies between 1 and 5")
    return x


MicroMag_data = collections.namedtuple(
    "Rawdata", "params calibration_measurement source_data NData")


def convert_lists_to_array(dict_with_data):
    for k in dict_with_data:
        dict_with_data[k] = np.array(dict_with_data[k])
    return dict_with_data


def add_data_to_dict(data, index, dict_name):
    if index in dict_name:
        dict_name[index].append(data)
    else:
        dict_name[index] = [data]


def smooth_array(in_array):
    """
    from forkit:
    "Smoothing is just a 5 point running average, which is done twice."
    """
    array = in_array.copy()
    N = len(array)
    if N > 4:
        for times_count in range(2):
            array[0] = np.sum(array[0:3]) / 3.0
            array[1] = np.sum(array[0:4]) * 0.25
            array[N - 1] = np.sum(array[N - 3:N]) / 3.0
            array[N - 2] = np.sum(array[N - 4:N]) * 0.25
            weight = np.ones(5)
            tmp = np.convolve(weight / weight.sum(), array, mode='valid')
            array = np.concatenate([array[0:2], tmp, array[N - 2:N]])
    return array


def rridge(forc_params, HaHrM_data):
    """
    from forkit
    "get the quasi-reversible ridge, which is the FORC distribution along
     the Bb or Hb axis for Bc or Hc ~ 0."
    """

    h = HaHrM_data[:, 0]
    hr = HaHrM_data[:, 1]
    moment = HaHrM_data[:, 2]

    # from forkit: "Hb=Ha=Hr along the ridge:
    hb = hr
    dm_dh = np.diff(moment) / np.diff(h)
    dm_dh = smooth_array(dm_dh)

    # from forkit:
    # "Take the second derivative (See Equation 1 of Pike, Physical Rev. B,
    # 68, 2003) Note, the equation below assumes that the FORCs are extended with
    # constant moment for Ba < Br using the moment value measured at Ba=Br.
    # This would give dm/dh = 0 for all grid points where Ba < Br. Thus, the
    # second derivative is just equal to dmdh(k)/(hb(k)-hb(k-1)) because
    # [dmdh(k)-dmdh(k-1)] = dmdh(k) - O = dmdh(k)"

    dm2_dh_dhb = -0.5 * dm_dh / np.diff(hb)
    dm2_dh_dhb = smooth_array(dm2_dh_dhb)

    # Another method is using the spline
    # degree_spline = 2
    # spline = UnivariateSpline(h, moment, k=degree_spline)
    # dm_dh = spline.derivative(n=1)
    # dm2_dh2 = -0.5*spline.derivative(n=2)
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # # ax.plot(h, moment, 'ko', ms=2, label='data')
    # ax.plot(rev_h, dm_dh(rev_h), 'k', label='2th deg spline')
    # ax.plot(h[1:], np.diff(moment)/np.diff(h), 'g', label='np diff')
    # plt.show()

    # 1.02 was introduced for plotting purpose
    dnormalizer = np.max(dm2_dh_dhb) * 1.02
    ddnorm = dm2_dh_dhb / dnormalizer

    with open("rridge.out", "w") as rridge_file:
        for i in range(dm_dh.shape[0]):
            print(
                "{:e} {:e} {:e} {:e}".format(hb[i + 1], dm_dh[i],
                                             dm2_dh_dhb[i], ddnorm[i]),
                file=rridge_file)


def plot_B_vs_M(path_for_saving, source_data, smoment_list):
    cm_subsection = np.linspace(0, 1, len(source_data.keys()))
    colors = [matplotlib.cm.YlGnBu(x) for x in cm_subsection]
    fig = plt.figure(figsize=(6.4, 6))
    gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[1, 3])
    ax_BvsM_text = plt.subplot(gs[0])
    ax = plt.subplot(gs[1])

    fig_txt = ("Sample ID and related information\n\n"
               "Sample ID: {} \n"
               "Maximum Moment: {:.1e} ($A m^2$)").format(
                   params['sampid'], Mnorm)
    ax_BvsM_text.set_ylim(0, 0.5)
    ax_BvsM_text.text(0.0, 0.35, fig_txt, fontsize=12, multialignment='left')
    ax_BvsM_text.set_axis_off()

    ax.set_title("FORC Paths in Hysteresis Space")
    ax.set_xlabel("B ($mT$)")
    ax.set_ylabel("M ($A m^2$)")
    ax.grid(color='whitesmoke')
    for key, data in source_data.items():
        if key % 4 == 0:
            ax.plot(data[:, 0], smoment_list[key], color=colors[key], lw=2)

    print(path_for_saving + "B_vs_M.pdf")
    fig.savefig(path_for_saving + "B_vs_M.pdf")


def plot_Ba_vs_Br(path_for_saving, source_data, hr):
    fig = plt.figure(figsize=(6.4, 6))
    gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[1, 3])
    ax_text = plt.subplot(gs[0])
    ax = plt.subplot(gs[1])

    # fig_BvsM, (ax_BvsM_text, ax_BvsM) = plt.subplots(2, 1, figsize=(6.4, 9))
    fig_txt = "Sample ID: {}".format(params['sampid'])
    ax_text.set_ylim(0, 0.5)
    ax_text.text(0.0, 0.35, fig_txt, fontsize=12, multialignment='left')
    ax_text.set_axis_off()

    ax.set_title("FORC Paths in Magnetic Field Space")
    ax.set_xlabel("Ba ($mT$)")
    ax.set_ylabel("Br ($mT$)")
    for key, data in source_data.items():
        y = np.ones(data.shape[0]) * hr[key]
        ax.plot(data[:, 0], y, color='grey', lw=1)

    fig.savefig(path_for_saving + "Ba_vs_Br.pdf")


def plot_drift(path_for_saving, hdrift, mdrift, smdrift):

    fig = plt.figure(figsize=(6.4, 6))
    gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[1, 3])
    ax_text = plt.subplot(gs[0])
    ax = plt.subplot(gs[1])

    # fig_BvsM, (ax_BvsM_text, ax_BvsM) = plt.subplots(2, 1, figsize=(6.4, 9))
    fig_txt = "Sample ID: {}".format(params['sampid'])
    ax_text.set_ylim(0, 0.5)
    ax_text.text(0.0, 0.35, fig_txt, fontsize=12, multialignment='left')
    ax_text.set_axis_off()

    ax.set_title("Drift of Moment (black) & Variation of Field (blue)")
    ax.set_xlabel("Ba ($mT$)")
    ax.set_ylabel("Relative Drift")
    x = range(len(hdrift))
    ax.plot(x, hdrift, color='b')
    ax.plot(x, mdrift, color='k', label='Drift of moment')
    ax.plot(x, smdrift, color='g', label='Smoothed drift of moment')
    ax.legend(loc=0)
    fig.subplots_adjust(left=0.18)
    fig.savefig(path_for_saving + "drift.pdf")


def plot_Ba_Br_moment(path_for_saving, HaHrM_data):

    h = HaHrM_data[:, 0]
    hr = HaHrM_data[:, 1]
    moment = HaHrM_data[:, 2]

    import matplotlib.tri as mtri
    fig = plt.figure(figsize=(6.4, 6))
    gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[1, 3])
    ax_text = plt.subplot(gs[0])
    ax = plt.subplot(gs[1])

    fig_txt = ("Sample ID and related information\n\n"
               "Sample ID: {} \n"
               "Maximum Moment: {:.1e} ($A m^2$)").format(
                   params['sampid'], Mnorm)
    ax_text.set_ylim(0, 0.5)
    ax_text.text(0.0, 0.35, fig_txt, fontsize=12, multialignment='left')
    ax_text.set_axis_off()

    triang = mtri.Triangulation(h, hr)
    cf = ax.tripcolor(
        triang, moment, shading='flat', cmap=plt.get_cmap('nipy_spectral'))
    fig.colorbar(cf)
    ax.set_title("Moment Variation")
    ax.set_xlabel("Ba ($mT$)")
    ax.set_ylabel("Br ($mT$)")
    fig.savefig(path_for_saving + "Ba_vs_Br_vs_moment.pdf")


def quadric_func(x, a1, a2, a3, a4, a5, a6):
    res = a1 + a2 * x[0] + a3 * x[0] * x[0] + a4 * x[1] + a5 * x[1] * x[1] + a6 * x[0] * x[1]
    return res.ravel()


def find_approximation_of_rho(x_grid, y_grid, grid_values):
    """
    function finds value of rho from article 
    Roberts, Pike and Verosub, Journal of Geophysical Research, VOL. 105, N. B12
    """
    # todo: replace on polyvander and leastsqr
    from scipy.optimize import curve_fit

    x_1d = x_grid.flatten()  # or ravel - that is question
    y_1d = y_grid.flatten()
    xdata = np.vstack((x_1d, y_1d))

    popt, pcov = curve_fit(quadric_func, xdata, grid_values.ravel())

    # fig, ax = plt.subplots()
    # ax.contourf(x_grid, y_grid, grid_values, 50)
    # ax.plot(x_grid.flatten(), y_grid.flatten(), "ok", label = "Data")
    # fitted = quadric_func(xdata, *popt)
    # fig2, ax2 = plt.subplots()
    # ax2.contourf(x_grid, y_grid, fitted.reshape(x_grid.shape), 50)
    # plt.show()

    return -popt[-1]


def plot_forc2D(path_for_saving, HaHrM_data, smoothing_factor=3):

    print("Plotting forc2d has started. It takes a while...")

    h = HaHrM_data[:, 0]
    hr = HaHrM_data[:, 1]
    moment = HaHrM_data[:, 2]

    # experimental data has a little variation but
    # data on regular grid is needed (for a simplification of 2D interpolation) therefore mean distance is calculated.
    # cKDTree is simple and fast way to find neighbourhoods and distances
    import scipy
    XY = np.c_[h, hr]  # it is needed to convert data for use with KDTree
    spatial_tree = scipy.spatial.cKDTree(XY)
    distances, points = spatial_tree.query(XY, k=2, p=1)
    # print (distances)
    step = np.mean(distances[:, 1:])

    minh, maxh = np.min(h), np.max(h)
    uni_h = np.linspace(
        minh, maxh, num=np.ceil(abs(maxh - minh) / step), endpoint=True)

    minhr, maxhr = np.min(hr), np.max(hr)
    uni_hr = np.linspace(
        minhr, maxhr, num=np.ceil(abs(maxhr - minhr) / step), endpoint=True)

    X, Y = np.meshgrid(uni_h, uni_hr)
    M = scipy.interpolate.griddata((h, hr), moment, (X, Y), method="nearest")
    SF = smoothing_factor
    rho = np.zeros(M.shape)
    win_withd = uni_hr.shape[0] // 2
    for i in range(M.shape[0])[SF:-SF]:
        for j in range(i, M.shape[1])[SF:-SF]:
            if win_withd - i + j < M.shape[1] - SF and i < 2 * M.shape[0] - j:
                i_slice = slice(i - SF, i + SF + 1)
                j_slice = slice(j - SF, j + SF + 1)
                rho[i, j] = find_approximation_of_rho(X[i_slice, j_slice],
                                                      Y[i_slice, j_slice],
                                                      M[i_slice, j_slice])

    # plt.contourf(X, Y, rho, 50)

    rho_new = np.zeros((win_withd, win_withd))
    newX = np.zeros((win_withd, win_withd))
    newY = np.zeros((win_withd, win_withd))
    for i in range(rho_new.shape[0]):
        for j in range(rho_new.shape[1]):
            a = win_withd - i + j
            b = win_withd + i + j
            rho_new[i, j] = rho[a, b]
            newX[i, j] = 0.5 * (X[a, b] - Y[a, b])
            newY[i, j] = 0.5 * (X[a, b] + Y[a, b])

    fig = plt.figure(figsize=(6.4, 6))
    gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[1, 3])
    ax_text = plt.subplot(gs[0])
    ax = plt.subplot(gs[1])

    fig_txt = ("Sample ID and related information\n\n"
               "Sample ID: {} \n"
               "Smoothing factor: {}").format(params['sampid'], SF)
    ax_text.set_ylim(0, 0.5)
    ax_text.text(0.0, 0.35, fig_txt, fontsize=12, multialignment='left')
    ax_text.set_axis_off()
    from matplotlib.ticker import MaxNLocator
    from matplotlib.colors import BoundaryNorm
    levels = MaxNLocator(nbins=50).tick_values(rho_new.min(), rho_new.max())
    cmap = plt.get_cmap('cool')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cf = ax.contourf(newX, newY, rho_new, levels=levels, cmap=cmap, norm=norm)
    cf2 = ax.contour(
        cf, levels=cf.levels[::5], colors='k', origin='lower', linewidths=0.25)
    fig.colorbar(cf)
    ax.plot(newX.flatten(), newX.flatten() * 0, color='k', linewidth=0.25)
    ax.set_title("FORC Diagram")
    ax.set_xlabel("Bu ($mT$)")
    ax.set_ylabel("Bc ($mT$)")
    fig.savefig(path_for_saving + "forc2D.pdf")


def write_out_files(path_for_saving, MicroMag_data, HaHrM_data_raw, HaHrM_data,
                    hdrift, smdrift, smoment_list):
    path_for_saving = path_for_saving + "out_files/"
    try:
        os.makedirs(path_for_saving)
    except OSError:
        if not os.path.isdir(path_for_saving):
            raise

    calibration_measurement = MicroMag_data.calibration_measurement
    source_data = MicroMag_data.source_data

    with open("drift.raw", "w") as drift_file:
        for key, data in calibration_measurement.items():
            print(
                "{} {:e} {:e} {:e} {:e}".format(
                    key + 1, hdrift[key], mdrift[key], data[0, 0], data[0, 1]),
                file=drift_file)

    with open(path_for_saving + "Ha-vs-Hr-vs-M.raw", "w") as Ha_Hr_M_file:
        for elem in HaHrM_data_raw:
            print(
                "{:e} {:e} {:e} {} {} {}".format(elem[0], elem[1], elem[2],
                                                 elem[3], elem[4], elem[5]),
                file=Ha_Hr_M_file)

    with open(path_for_saving + "drift.out", "w") as drift_file:
        for key, data in calibration_measurement.items():
            print(
                "{} {:e} {:e} {:e} {:e}".format(key + 1, hdrift[key],
                                                smdrift[key], data[0, 0],
                                                data[0, 1]),
                file=drift_file)

    with open(path_for_saving + "Ha-vs-Hr-vs-M.out", "w") as Ha_Hr_M_file:
        for elem in HaHrM_data:
            print(
                "{:e} {:e} {:e} {} {} {}".format(elem[0], elem[1], elem[2],
                                                 elem[3], elem[4], elem[5]),
                file=Ha_Hr_M_file)

    with open(path_for_saving + "H-vs-M.out", "w") as H_M_file:
        for key, data in source_data.items():
            title = "> New FORC"
            if int(key) == 0:
                title = "> First FORC (a single point)"
            print(title, file=H_M_file)
            smoment = smoment_list[key]
            for i in range(data.shape[0]):
                print(
                    "{:e} {:e}".format(data[i, 0], smoment[i]), file=H_M_file)

    with open(path_for_saving + "Ha-vs-Hr.out", "w") as Ha_Hr_file:
        for key, data in source_data.items():
            title = "> New FORC"
            if int(key) == 0:
                title = "> First FORC (a single point)"
            print(title, file=Ha_Hr_file)
            for i in range(data.shape[0]):
                print("{:e} {:e}".format(data[i, 0], hr[key]), file=Ha_Hr_file)

    rridge(MicroMag_data.params, np.array(HaHrM_data))


def read_micromag_file(micromag_filename):
    # from forkit:
    # "Read Micromag data, starting with the header lines
    # The first test determines if the data is in the old
    # format, which has a single header line with 14 parameters.
    # In the new format, the first line begins with the word MicroMag"
    micromag_file = open(micromag_filename, "r")
    first_line = micromag_file.readline().strip()
    if not "MicroMag 2900/3900 Data File" == first_line:
        #http://www.conrad-observatory.at/cmsjoomla/phocadownload/variforc_usermanual_chapter03.pdf
        raise NotImplementedError(
            "You should rewrite this file to take into account new format")

    # Read in header information for data files with new format
    params = {}

    is_data_begin = False

    while not is_data_begin:
        line = micromag_file.readline().strip()
        if line:
            if line[0] == "\"":
                params["sampid"] = line
            elif ":" in line:
                parameter, value = [x.strip() for x in line.split(":")]
                params[parameter] = value
            elif "=" in line:
                parameter, value = [x.strip() for x in line.split("=")]
                if not "N/A" in value:
                    value = float(value)
                params[parameter] = value
            elif line[0] == "+" or line[0] == "-":
                is_data_begin = True


    MicroMag_data.params = params

    try:
        MicroMag_data.NData = int(params["NData"])
    except KeyError:
        print(("Header of you file does not contain 'NData'."
               "May be it is time to rewrite parser"))
        import sys
        sys.exit()

    # calibration_measurement is dictionary for calibration measurement data
    # source_data is dictionary for calibration FORC data
    calibration_measurement = {
        0: [[float(x.strip()) for x in line.split(",")]]
    }
    source_data = {}
    nonempty_count = 1
    empty_count = 0
    is_calibration_data = True

    for line in micromag_file:
        tmp = line.strip()
        if tmp and not "Data File ends" in tmp:
            nonempty_count += 1
            data = [float(x.strip()) for x in line.split(",")]
            idx = int(empty_count / 2)
            if is_calibration_data:
                add_data_to_dict(data, idx, calibration_measurement)
            else:
                add_data_to_dict(data, idx, source_data)
        else:
            empty_count += 1
            is_calibration_data = not is_calibration_data
    micromag_file.close()

    if nonempty_count - MicroMag_data.NData > 0:
        raise ValueError(
            ("Number of lines with data does not correspond",
             " to NData field in *.frc file. Please, check your file")
        )

    MicroMag_data.source_data = convert_lists_to_array(source_data)
    MicroMag_data.calibration_measurement = convert_lists_to_array(
        calibration_measurement)
    
    print ("frc-file reading has finished.")
    return MicroMag_data


if __name__ == "__main__":

    cmd_parser = create_parser()
    cmd_arguments = cmd_parser.parse_args()
    micromag_filename = cmd_arguments.filename
    dir_name = os.path.dirname(micromag_filename) + "/"
    if dir_name == "/":
        dir_name = ""

    MicroMag_data = read_micromag_file(micromag_filename)
    calibration_measurement = MicroMag_data.calibration_measurement
    source_data = MicroMag_data.source_data
    params = MicroMag_data.params

    if cmd_arguments.plot_raw_data:
        cm_subsection = np.linspace(0, 1, len(source_data.keys()))
        # http://matplotlib.org/users/colormaps.html
        colors = [matplotlib.cm.YlGnBu(x) for x in cm_subsection]
        fig, ax = plt.subplots()
        for k in source_data.keys():
            data = source_data[k]
            ax.plot(data[:, 0], data[:, 1], color=colors[k])
        fig.savefig(dir_name + "raw_reversals.pdf")

    units_of_measure = params.get("Units of measure",
                                  "cgs")  # cgs is default value

    # from forcit:
    # "Note the conversions are from SI Teslas to SI mT with moment remaining in A m^2
    # or are from cgs (Oe and emu) to SI (mT and A m^2)."
    params_for_change = ["Hb1", "Hb2", "Hc1", "Hc2", "HCal", "HNcr", "HSat"]
    if "SI" in units_of_measure:
        for p in params_for_change:
            params[p] *= 1000.0
        params["Slope corr."] *= 1000.0

        for key, data in source_data.items():
            data[:, 0] = data[:, 0] * 1000.0
        for key, data in calibration_measurement.items():
            data[:, 0] *= 1000.0
    else:
        for p in params_for_change:
            params[p] *= 0.1
        params["Moment range"] *= 0.001
        params["Slope corr."] *= 0.01
        for key, data in source_data.items():
            data[:, 0] *= 0.1
            data[:, 1] *= 0.001
        for key, data in calibration_measurement.items():
            data[:, 0] *= 0.1
            data[:, 1] *= 0.001

    slope = params["Slope corr."]

    # it is time create data structures for plotting

    h = np.zeros(MicroMag_data.NData)
    FORC_number = len(source_data)
    hr = np.zeros(FORC_number)
    hdrift = np.zeros(FORC_number)
    mdrift = np.zeros(FORC_number)

    # from forcit:
    # "This sets the maximum moment to 2% larger than the first drift point
    # All other moments are normalized by this value."
    denom = calibration_measurement[0]
    Mnorm = (denom[0, 1] + slope * denom[0, 0]) * 1.02

    for key, data in calibration_measurement.items():
        hdrift[key] = data[:, 0] / denom[0, 0]
        mdrift[key] = data[:, 1] / denom[0, 1]


    moments_list = []
    total_rec_counter = 0

    HaHrM_data_raw = []
    for key, data in source_data.items():
        moment = (data[:, 1] + slope * data[:, 0]) / Mnorm
        moments_list.append(moment)
        hr[key] = data[0, 0]
        for i in range(data.shape[0]):
            h[total_rec_counter] = data[i, 0]
            total_rec_counter += 1
            HaHrM_data_raw.append([
                data[i, 0], hr[key], moment[i], key + 1, i + 1,
                total_rec_counter
            ])

    smdrift = smooth_array(mdrift)
    smoment_list = []

    total_rec_counter = 0
    HaHrM_data = []
    for key, data in source_data.items():
        smoment = moments_list[key] / smdrift[key]
        smoment_list.append(smoment)
        for i in range(data.shape[0]):
            total_rec_counter += 1
            HaHrM_data.append([
                data[i, 0], hr[key], smoment[i], key + 1, i + 1,
                total_rec_counter
            ])

    # todo: write tests and uncoment next line
    # write_out_files(dir_name, MicroMag_data, HaHrM_data_raw, HaHrM_data, hdrift, smdrift, smoment_list)

    if not cmd_arguments.only_forc:
        HaHrM_data = np.array(HaHrM_data)
        plot_B_vs_M(dir_name, MicroMag_data.source_data, smoment_list)
        plot_Ba_vs_Br(dir_name, MicroMag_data.source_data, hr)
        plot_drift(dir_name, hdrift, mdrift, smdrift)
        plot_Ba_Br_moment(dir_name, HaHrM_data)

    plot_forc2D(dir_name, HaHrM_data, cmd_arguments.smoothing_factor)

    if not dir_name:
        dir_name = "current folder"
    print(
        "Images have been saved near frc-file, in {}".format(dir_name))
