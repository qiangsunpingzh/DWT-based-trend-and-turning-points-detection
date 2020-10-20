##########################################DWT-TTD###############################################################
"""
# -*- coding: utf-8 -*-
This code was used to figure out the processes of DWT-TTD for a given time series of signals from EXCEL
Author: qiangsun@cau.edu.cn  
College of Land Science and Technology, China Agricultural University, Beijing, China
"""
##########################################DWT-TTD###############################################################
# packages
import pywt
import numpy as np
import xlrd
import pywt.data
import matplotlib.pyplot as plt
from scipy import optimize

def DWT_A(mw,ts,lt):
    """
    :param mw: mother Wavelet, i.e., dmey
    :param ts: a time series signal
    :param lt: level of WT based on length of signal
    :return: approximation series of DWT
    """
    mode = pywt.Modes.sym
    w = pywt.Wavelet(mw)
    a = ts
    ca = []
    for i in range(lt):
        (a, d) = pywt.dwt(a, w, mode)
        ca.append(a)
    rec_a = []
    for i, coeff in enumerate(ca):
        coeff_list = [coeff, None] + [None] * i
        rec_a.append(pywt.waverec(coeff_list, w))
    return rec_a

def DWT_D(mw,ts,lt):
    """
    :param mw: mother Wavelet, i.e., dmey
    :param ts: a time series signal
    :param lt: level of WT based on length of signal
    :return: detail series of DWT
    """
    mode = pywt.Modes.sym
    w = pywt.Wavelet(mw)
    a = ts
    cd = []
    for i in range(lt):
        (a, d) = pywt.dwt(a, w, mode)
        cd.append(d)
    rec_d = []
    for i, coeff in enumerate(cd):
        coeff_list = [None, coeff] + [None] * i
        rec_d.append(pywt.waverec(coeff_list, w))
    return rec_d

def selected_A(mw,ts,lt, selected_level):
    for i, y in enumerate(DWT_A(mw, ts, lt)):
        if i == selected_level-1:
            return y

def selected_D(mw,ts,lt, selected_level):
    for i, y in enumerate(DWT_D(mw, ts, lt)):
        if i == selected_level-1:
            return y

def TTD_trend(mw,ts,lt, level):
    """
    :param mw: mother Wavelet, i.e., dmey
    :param ts: a time series signal
    :param lt: level of WT based on length of signal
    :param level: level of approximation series selected for trend analysis
    :return: [k,b, zc] are slope, intercept, and M-K test results, respectively
    """
    import math
    X = np.linspace(1,len(ts),len(ts))
    for i, y in enumerate(DWT_A(mw,ts,lt)):
        if i == level-1:
            Y = y[0:len(ts)]
            s = 0
            length = len(Y)
            # M-K
            for m in range(0, length - 1):
                for n in range(m + 1, length):
                    if Y[n] > Y[m]:
                        s = s + 1
                    elif Y[n] == Y[m]:
                        s = s + 0
                    else:
                        s = s - 1
            vars = length * (length - 1) * (2 * length + 5) / 18
            if s > 0:
                zc = (s - 1) / math.sqrt(vars)
            elif s == 0:
                zc = 0
            else:
                zc = (s + 1) / math.sqrt(vars)
            # Slope
            def residuls(p):
                k, b = p
                return Y - (k * X + b)

            r = optimize.leastsq(residuls, [1, 0])
            k, b = r[0]
            return [k,b, zc]

def cluster(items, key_func):
    items = sorted(items)
    clusters = [[items[0]]]
    for item in items[1:]:
        cluster = clusters[-1]
        last_item = cluster[-1]
        if key_func(item, last_item):
            cluster.append(item)
        else:
            clusters.append([item])
    return clusters

def TTD_timing(mw,ts,lt, level):
    """
    :param mw: mother Wavelet, i.e., dmey
    :param ts: a time series signal
    :param lt: level of WT based on length of signal
    :param level: level of approximation series selected for trend analysis
    :return: timing indexs of turning points
    """
    Timing = []
    for j, d in enumerate(DWT_D(mw,ts,lt)):
        if j == level-1:
            for m in range(0, len(data) + 1):
                if m == 0 or m == len(data):
                    Timing.append(m)
                else:
                    if np.sign(d[m] - d[m - 1]) == -np.sign(d[m + 1] - d[m]):
                        Timing.append(m)

    # for selected valley/peak points, using cluster approach eliminate errors due to peak and valley
    # fluctuations. from each cluster, we select median as final turning points
    Timing_cluster = cluster(Timing, lambda curr, prev: curr - prev < 2)

    Timing_median = []
    for i in Timing_cluster:
        Timing_median.append(int(np.median(i)))

    return Timing_median

def TTD_type(mw,ts,lt, level):
    FX0 = TTD_timing(mw, ts, lt, level)
    for i, a in enumerate(DWT_A(mw,ts,lt)):
        if i == level-2:
            ts = []
            for n in FX0:
                ts.append(a[n])

            # using successive turning points to identify types using discriminative function. 1,2,3,4 is peak, valley,
            # interrupted decrease, interrupted increase. the 0 and end is set as "0"
            Type = []
            for n in range(len(FX0)):
                if n < len(FX0) - 1 and n > 0:
                    if np.sign(ts[n + 1] - ts[n]) == -np.sign(ts[n] - ts[n - 1]) == -1:
                        Type.append(1)
                    elif np.sign(ts[n + 1] - ts[n]) == -np.sign(ts[n] - ts[n - 1]) == 1:
                        Type.append(2)
                    elif np.sign(ts[n + 1] - ts[n]) == np.sign(ts[n] - ts[n - 1]) == -1:
                        Type.append(3)
                    elif np.sign(ts[n + 1] - ts[n]) == np.sign(ts[n] - ts[n - 1]) == 1:
                        Type.append(4)


            return Type


def TTD_magnitude(mw,ts,lt, level):
    FX0 = TTD_timing(mw, ts, lt, level)
    for i, a in enumerate(DWT_A(mw,ts,lt)):
        if i == level-2:
            ts = []
            for n in FX0:
                ts.append(a[n])

            mag = []
            for i3 in range(len(ts)):
                if 0 < i3 < len(ts) - 1:
                    magi = ts[i3 + 1] - ts[i3]
                    mag.append(magi)
            return mag


def TTD_RMSE(mw, ts,lt, level):
    for i, a in enumerate(DWT_A(mw,ts,lt)):
        if i == level-2:
            y = a[0:391]
            x = np.linspace(1, len(ts), len(ts))
            FX0 = TTD_timing(mw, ts, lt, level)
            RMSE = []
            for i in FX0[1:][:-1]:
                def piecewise_linear1(x, k1, k2):
                    return np.piecewise(x, [x < i, x >= i], [lambda x: k1 * x + y[i] - k1 * i,
                                                         lambda x: k2 * x + y[i] - k2 * i])

                p, e = optimize.curve_fit(piecewise_linear1, x, y)
                yt = piecewise_linear1(x, *p)
                RSS = np.sum([(yt[i] - y[i]) * (yt[i] - y[i]) for i in range(0, 390, 1)])
                rmse = np.sqrt(RSS / len(y))
                RMSE.append(rmse)

            return RMSE

def BIC(ts,predict,k):
    RSS = np.sum([(ts[i] - predict[i]) * (ts[i] - predict[i]) for i in range(0, 390, 1)])
    mean_predict = np.mean(predict)
    var = np.sum([(predict[i] - mean_predict) * (predict[i] - mean_predict) for i in range(0, 390, 1)])/len(predict)
    bic = k*np.log(len(ts))+2*len(ts)*np.log(2.506*np.sqrt(var))+RSS/var
    return bic

def TTD_BIC(mw, ts,lt, level):
    mag = TTD_magnitude(mw, ts, lt, level)
    RMSE = TTD_RMSE(mw, ts, lt, level)
    Type = TTD_type(mw, ts, lt, level)
    FX0 = TTD_timing(mw, ts, lt, level)[1:][:-1]
    abs_mag = np.abs(mag)
    zipped = zip(FX0, abs_mag)
    sort_zipped = sorted(zipped, key=lambda x: (x[1], x[0]), reverse=True)
    result = zip(*sort_zipped)
    timing_sorted, abs_mag_sorted = [list(x) for x in result]
    for i, a in enumerate(DWT_A(mw,ts,lt)):
        if i == level-2:
            # model selection based on BIC of each case
            # case1: linear regression
            # case2: one turning points with max change magnitude
            # case2: two turning points with max change magnitude
            # case3: ..................
            BIC_CASE = []
            case = np.linspace(0,len(timing_sorted),len(timing_sorted)+1)
            for case0 in case:
                x0 = np.linspace(1, len(ts), len(ts))
                y0 = a[0:391]

                if case0 == 0:
                    # case0
                    def residuls(p):
                        k, b = p
                        return y0 - (k * x0 + b)

                    r = optimize.leastsq(residuls, [1, 0])
                    k, b = r[0]
                    p0 = [k * i + b for i in x0]
                    BIC_case0 = BIC(y0, p0, 2)
                    BIC_CASE.append(BIC_case0)

                if case0 > 0:
                    tp = sorted(timing_sorted[0:int(case0)])
                    # fit with known breakpoint locations
                    y_tp = []
                    for j in tp:
                        y_tp.append(y0[j])
                    import pwlf
                    tp.append(390)
                    tp.insert(0,0)
                    xtp = np.array(tp)
                    # initialize piecewise linear fit with your x and y data
                    my_pwlf = pwlf.PiecewiseLinFit(x0, y0)
                    # will terminate)
                    # my_pwlf.fit_with_breaks_force_points(xtp, tp, y_tp)
                    my_pwlf.fit_with_breaks(xtp)
                    # predict for the determined points
                    yHat = my_pwlf.predict(x0)
                    paranumber = case0+2
                    bic_case = BIC(y0, yHat, paranumber)
                    BIC_CASE.append(bic_case)

            return BIC_CASE

def TTD_valid(mw, ts,lt, level):
    mag = TTD_magnitude(mw, ts, lt, level)
    RMSE = TTD_RMSE(mw, ts, lt, level)
    Type = TTD_type(mw, ts, lt, level)
    FX0 = TTD_timing(mw, ts, lt, level)[1:][:-1]
    abs_mag = np.abs(mag)
    if len(FX0)>0:
        zipped = zip(FX0, abs_mag)
        sort_zipped = sorted(zipped, key=lambda x: (x[1], x[0]), reverse=True)
        result = zip(*sort_zipped)
        timing_sorted, abs_mag_sorted = [list(x) for x in result]

        BIC_CASE = TTD_BIC(mw, ts,lt, level)
        index = BIC_CASE.index(min(BIC_CASE))
        valid_tp = timing_sorted[0:index]


        valid_rmse =[]
        valid_TYPES = []
        valid_mag = []
        # all of valid turning points
        for i in valid_tp:
            id = FX0.index(i)
            valid_rmse.append(RMSE[id])
            valid_TYPES.append(Type[id])
            valid_mag.append(mag[id])
    else:
        valid_tp = []
        valid_rmse = []
        valid_TYPES = []
        valid_mag = []
    return [valid_tp, valid_rmse, valid_TYPES, valid_mag]

def TTD_one(mw, ts, lt, level):
    valid_tp = TTD_valid(mw, ts,lt, level)[0]
    valid_rmse = TTD_valid(mw, ts,lt, level)[1]
    valid_TYPES = TTD_valid(mw, ts,lt, level)[2]
    valid_mag = TTD_valid(mw, ts,lt, level)[3]
    if len(valid_tp) ==0:
        one_tp = 0
        one_rmse = 0
        one_TYPES = 0
        one_mag = 0
    else:
        # select one using RMSE from all of valid turning points
        one_id = valid_rmse.index(min(valid_rmse))

        one_tp = valid_tp[one_id]
        one_rmse = valid_rmse[one_id]
        one_TYPES = valid_TYPES[one_id]
        one_mag = valid_mag[one_id]

    return [one_tp, one_rmse, one_TYPES, one_mag]


if __name__ == '__main__':
    Ex_data = xlrd.open_workbook("C:/Users/dell/Desktop/excel/modis.xlsx")
    Sh_data = Ex_data.sheets()[1]
    data = Sh_data.row_values(0)
    # mother wavelet and mode
    mode = pywt.Modes.sym
    w = 'dmey'
    """
    #------------------------------------------------print-------------------------------------------------------
    # print results
    print TTD_trend(w, data, 8, 8)
    print TTD_timing(w, data, 8, 7)
    print TTD_type(w, data, 8, 7)
    print TTD_magnitude(w, data, 8, 7)
    print TTD_BIC(w, data, 8, 7)
    print TTD_valid(w, data, 8, 7)
    print TTD_one(w, data, 8, 7)
    #------------------------------------------------print-------------------------------------------------------
    """
    """
    ------------------------------------------------plot--------------------------------------------------------
    # plot
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    fig = plt.figure(figsize=(3, 7))
    plt.rc('font', family='Times New Roman', size=10)
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    yformate = FormatStrFormatter('%1.3f')
    labels = ["2001", "2006", "2011", "2016"]
    # plot of original signal
    x = np.linspace(1, 391, 391)
    ax1 = plt.subplot(8, 1, 1)
    ax1.plot(x, data)
    ax1.set_xticks([n for n in range(12, 392, 115)])
    ax1.set_xlim(0, len(data) - 1)
    miny = np.min(data)
    maxy = np.max(data)
    meany = (miny + maxy) / 2
    ax1.set_ylim(miny - 0.005, maxy + 0.005)
    ax1.set_yticks([miny - 0.005, meany, maxy + 0.005])
    ax1.set_yticklabels([miny - 0.005, meany, maxy + 0.005])
    ax1.yaxis.set_major_formatter(yformate)
    ax1.set_xticklabels([])

    ax2 = plt.subplot(8, 1, 2)
    ax2.plot(x,selected_A(w, data, 8, 8)[0:391])
    ypt = [TTD_trend(w, data, 8, 8)[0]*i + TTD_trend(w, data, 8, 8)[1] for i in x]
    ax2.plot(x,ypt)
    ax2.set_xticks([n for n in range(12, 392, 115)])
    ax2.set_xlim(0, len(data) - 1)
    miny = np.min(ypt)
    maxy = np.max(ypt)
    meany = (miny + maxy) / 2
    ax2.set_ylim(miny - 0.005, maxy + 0.005)
    ax2.set_yticks([miny - 0.005, meany, maxy + 0.005])
    ax2.set_yticklabels([miny - 0.005, meany, maxy + 0.005])
    ax2.yaxis.set_major_formatter(yformate)
    ax2.set_xticklabels([])


    ax3 = plt.subplot(8, 1, 3)
    ax3.plot(x, selected_D(w, data, 8, 7)[0:391])
    ax3.set_xticks([n for n in range(12, 392, 115)])
    ax3.set_xlim(0, len(data) - 1)
    miny = np.min(selected_D(w, data, 8, 7)[0:391])
    maxy = np.max(selected_D(w, data, 8, 7)[0:391])
    meany = (miny + maxy) / 2
    ax3.set_ylim(miny - 0.005, maxy + 0.005)
    ax3.set_yticks([miny - 0.005, meany, maxy + 0.005])
    ax3.set_yticklabels([miny - 0.005, meany, maxy + 0.005])
    ax3.yaxis.set_major_formatter(yformate)
    ax3.set_xticklabels([])
    FX_D = [selected_D(w, data, 8, 7)[i] for i in TTD_timing(w, data, 8, 7)]
    ax3.scatter(TTD_timing(w, data, 8, 7),FX_D, s=5, c="black")

    ax4 = plt.subplot(8, 1, 4)
    ax4.plot(x,selected_A(w, data, 8, 7)[0:391])
    ts = [selected_A(w, data, 8, 7)[i] for i in TTD_timing(w, data, 8, 7)]
    ax4.scatter(TTD_timing(w, data, 8, 7), ts, s=5, c="black")
    ax4.plot(TTD_timing(w, data, 8, 7), ts)
    ax4.set_xticks([n for n in range(12, 392, 115)])
    ax4.set_xlim(0, len(data) - 1)
    miny = np.min(selected_A(w, data, 8, 7)[0:391])
    maxy = np.max(selected_A(w, data, 8, 7)[0:391])
    meany = (miny + maxy) / 2
    ax4.set_ylim(miny - 0.005, maxy + 0.005)
    ax4.set_yticks([miny - 0.005, meany, maxy + 0.005])
    ax4.set_yticklabels([miny - 0.005, meany, maxy + 0.005])
    ax4.yaxis.set_major_formatter(yformate)
    ax4.set_xticklabels([])


    ax5 = plt.subplot(8, 1, 5)
    ax5.scatter(TTD_timing(w, data, 8, 7)[1:][:-1], TTD_type(w, data, 8, 7), s=5, c="black")
    ax5.set_xticks([n for n in range(12, 392, 115)])
    ax5.set_xlim(0, len(data) - 1)
    miny = 1
    maxy = 4
    meany = 2.50
    ax5.set_ylim(miny - 0.5, maxy + 0.5)
    ax5.set_yticks([miny - 0.5, meany, maxy + 0.5])
    ax5.set_yticklabels([miny - 0.5, meany, maxy + 0.5])
    ax5.yaxis.set_major_formatter(yformate)
    ax5.set_xticklabels([])

    ax6 = plt.subplot(8, 1, 6)
    ax6.plot(x, x * 0, c="black")
    ax6.scatter(TTD_timing(w, data, 8, 7)[1:][:-1], TTD_magnitude(w, data, 8, 7), s=5, c="black")
    ax6.set_xticks([n for n in range(12, 392, 115)])
    ax6.set_xlim(0, len(data) - 1)
    miny = np.min(TTD_magnitude(w, data, 8, 7))
    maxy = np.max(TTD_magnitude(w, data, 8, 7))
    meany = (miny + maxy) / 2
    ax6.set_ylim(miny - 0.005, maxy + 0.005)
    ax6.set_yticks([miny - 0.005, meany, maxy + 0.005])
    ax6.set_yticklabels([miny - 0.005, meany, maxy + 0.005])
    ax6.yaxis.set_major_formatter(yformate)
    ax6.set_xticklabels([])

    ax7 = plt.subplot(8, 1, 7)
    RMSE = TTD_RMSE(w, data, 8, 7)
    ax7.scatter(TTD_timing(w, data, 8, 7)[1:][:-1], RMSE, s=5, c="black")
    ax7.set_xticks([n for n in range(12, 392, 115)])
    ax7.set_xlim(0, len(data) - 1)
    miny = np.min(RMSE)
    maxy = np.max(RMSE)
    meany = (miny + maxy) / 2
    ax7.set_ylim(miny - 0.005, maxy + 0.005)
    ax7.set_yticks([miny - 0.005, meany, maxy + 0.005])
    ax7.set_yticklabels([miny - 0.005, meany, maxy + 0.005])
    ax7.yaxis.set_major_formatter(yformate)
    ax7.set_xticklabels(labels)

    ax8 = plt.subplot(8, 1, 8)
    BIC = TTD_BIC(w, data, 8, 7)
    case = ['case%s'%i for i in range(0,len(BIC))]
    ax8.bar(case, BIC)
    miny = np.min(BIC)
    maxy = np.max(BIC)
    meany = (miny + maxy) / 2
    ax8.set_ylim(miny - 20, maxy + 20)
    ax8.set_yticks([miny - 20, meany, maxy + 20])
    ax8.set_yticklabels([miny - 20, meany, maxy + 20])
    ax8.yaxis.set_major_formatter(yformate)
    ax8.set_xticklabels(case)
    
    plt.show()
    ------------------------------------------------plot--------------------------------------------------------
    """

    # apply for time series of images
    import arcpy
    import os
    def raster(raster, output_folder, filename, blocksize):
        """
        :param raster: the time series of images inputted
        :param output_folder: 
        :param filename: 
        :param blocksize: 
        :return:
        """
        arcpy.env.overwriteOutput = True
        blocksize = blocksize
        P_documents_out = os.path.join(output_folder, r"block.tif")
        myRaster = arcpy.Raster(raster)
        arcpy.env.outputCoordinateSystem = raster
        TP_list = []
        blockno = 0
        for x in range(0, myRaster.width, blocksize):
            for y in range(0, myRaster.height, blocksize):
                # Lower left coordinate of block (in map units)
                mx = myRaster.extent.XMin + x * myRaster.meanCellWidth
                my = myRaster.extent.YMin + y * myRaster.meanCellHeight
                # Upper right coordinate of block (in cells)
                nx = min([x + blocksize, myRaster.width])
                ny = min([y + blocksize, myRaster.height])
                # Extract data block
                myData = arcpy.RasterToNumPyArray(myRaster, arcpy.Point(mx, my), nx - x, ny - y)
                TP_block = []
                for i in range(myData.shape[1]):
                    for j in range(myData.shape[2]):
                        ts = np.array(myData[:, i, j])
                        one = TTD_one(w, ts, 8, 7)
                        one.append(TTD_trend(w, ts, 8, 8)[0])
                        one.append(TTD_trend(w, ts, 8, 8)[2])
                        print one
                        TP_block.append(one)
                TP = np.array(TP_block)
                TP_data = TP.ravel(order = 'F').reshape(6,myData.shape[1],myData.shape[2])
                K0_Block = arcpy.NumPyArrayToRaster(TP_data, arcpy.Point(mx, my), myRaster.meanCellWidth,
                                                    myRaster.meanCellHeight)


                K0_temp = ('_%i.' % blockno).join(P_documents_out.rsplit('.', 1))
                K0_Block.save(K0_temp)
                TP_list.append(K0_temp)
                del myData
                del TP
                del TP_data
                del K0_Block
                blockno += 1

        arcpy.MosaicToNewRaster_management(";".join(TP_list[0:]), output_folder, filename, "", "64_BIT",
                                           "", "1", "FIRST", "FIRST")
        for K0_item in TP_list:
            if arcpy.Exists(K0_item):
                arcpy.Delete_management(K0_item)

        return filename


    print raster("G:\DWTTTD\sg_gv.tif", "G:\DWTTTD", "ttd.tif", 20)
