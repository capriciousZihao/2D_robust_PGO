# -*- coding: utf-8 -*-
import numpy as np
import time
import math
import pylab
#np.ones((j2.shape[0], 1))*   np.ones((j2.shape[0], 1)) *
def calErr42Trajectory(j1Truth, j2, result):
    j2[:, 1] = j2[:, 1]+(j1Truth[0, 1]-j2[0, 1])
    j2[:, 2] = j2[:, 2] + (j1Truth[0, 2] - j2[0, 2])

    print(j1Truth[0, 3], j2[0, 3])
    print(math.sin(j1Truth[0, 3] - j2[0, 3]))

    e00 = math.cos(np.mean(j1Truth[1:10, 3]) - np.mean(j2[1:10, 3]))
    e01 = math.sin(np.mean(j1Truth[1:10, 3]) - np.mean(j2[1:10, 3]))
    # e00 = math.cos((j1Truth[0, 3]) - (j2[0, 3]))
    # e01 = math.sin((j1Truth[0, 3]) - (j2[0, 3]))

    #rotate = np.mat('math.cos(j1Truth[0, 3]-j2[0, 3])  math.cos(j1Truth[0, 3]-j2[0, 3]); -math.sin(j1Truth[0, 3]-j2[0, 3])  math.cos(j1Truth[0, 3]-j2[0, 3])')
    for i in range(j2.shape[0]):
        print("original j2[%d] %d %d" % (i, j2[i, 1],j2[i,2]))
        print("cos:%f   sin:%f" % (e00, e01))
        x = j2[i, 1]
        y = j2[i, 2]
        j2[i, 1] = x * e00 - y * e01

        j2[i, 2] = x * e01 + y * e00
        print("transformed j2[%d] %d %d" % (i, j2[i, 1], j2[i, 2]))




    if j1Truth.shape != j2.shape :
        multiplyCoeffi = round(j2.shape[0]/float(j1Truth.shape[0])*100-1)
        print ("基准矩阵和待求矩阵的维度不一致")
        time.sleep(3)
        for i in range(0, j1Truth.shape[0]):
            a, b, c = j1Truth[i, [1, 2, 3]] - j2[round(i*multiplyCoeffi/100), [1, 2, 3]]  # 顶点上的X轴误差
            #b = j1Truth[i, 2] - j2[i, 2]  # 顶点上的Y轴误差
            #c = j1Truth[i, 3] - j2[i, 3]  # 顶点上的方向误差
            d = np.sqrt(a * a + b * b)  # 每一个顶点上的位置误差
            result[i, :] = [a, b, c, d]
    else :
        for i in range(0, j1Truth.shape[0]):
            if j1Truth[i, 0] - j2[i, 0] >= 0.001:
                print("计算误差时两个Vertex序号不一致")
                exit(0)
            a = j1Truth[i, 1] - j2[i, 1]  # 顶点上的X轴误差
            b = j1Truth[i, 2] - j2[i, 2]  # 顶点上的Y轴误差
            c = j1Truth[i, 3] - j2[i, 3]  # 顶点上的方向误差
            d = np.sqrt(a * a + b * b)  # 每一个顶点上的位置误差
            result[i, :] = [a, b, c, d]



def setTitle(fileName):
    if fileName.find(MM) != -1:
        algorithm = Max-Mixture
        if (fileName.find(DSC) != -1) or (fileName.find(RRR) != -1):
            print("file name contain more than one algorithm, can not determine the right name")
            time.sleep(3)
    elif fileName.find(DSC) != -1:
        algorithm = DSC
        if fileName.find(RRR) != -1:
            print("file name contain more than one algorithm, can not determine the right name")
            time.sleep(3)
    elif fileName.find(RRR) != -1:
        algorithm = RRR


    title = algorithm+parameter+dataset
    return title


def loadData(flieName):
    inFile = open(flieName, 'r')  # 以只读方式打开某fileName文件
    count = len(inFile.readlines())  # 计算行数

    # 定义两个空list，用来存放文件中的数据
    X = []
    y = []

    for line in inFile:
        trainingSet = line.split(' ')  # 对于每一行，按','把数据分开，这里是分成两部分
        X.append(trainingSet[2])  # 第一部分，即文件中的第一列数据逐一添加到list X 中
        y.append(trainingSet[3])  # 第二部分，即文件中的第二列数据逐一添加到list y 中

    return (X, y)  # X,y组成一个元组，这样可以通过函数一次性返回


# !/usr/bin/python
# Requirements:
# sudo apt-get install python-argparse

"""
This script computes the absolute trajectory error from the ground truth
trajectory and the estimated trajectory.
"""


def align(mt, ms):
    """Align two trajectories using the method of Horn (closed-form).

    Input:
    model -- first trajectory (3xn)
    data -- second trajectory (3xn)

    Output:
    rot -- rotation matrix (3x3)
    trans -- translation vector (3x1)
    trans_error -- translational error per point (1xn)

    """
    model = mt.transpose()
    data  = ms.transpose()
    if model.shape[0] == 2:
        print("input data is se2")
        addM = np.zeros([1, model.shape[1]])

        print("model shape: ")
        print(model.shape)

        print("addM shape: ")
        print(addM.shape)

        model = np.r_[model, addM]
        data = np.r_[data, addM]

    dimension = model.shape[0]
    np.set_printoptions(precision=3, suppress=True)
    # modelTanspose = model.transpose()
    # print(modelTanspose.shape)
    # dataTranspose = data.transpose()
    model_zerocentered = model - model.mean(1).reshape(3, 1)
    data_zerocentered =  data - data.mean(1).reshape(3, 1)
    if dimension == 3:
        W = np.zeros((3, 3))
        for column in range(model.shape[1]):
            W += np.outer(model_zerocentered[:, column], data_zerocentered[ :, column])
        U, d, Vh = np.linalg.linalg.svd(W.transpose())
        S = np.matrix(np.identity(3))
        if (np.linalg.det(U) * np.linalg.det(Vh) < 0):
            S[2, 2] = -1
        rot = U * S * Vh
        trans = data.mean(1).reshape(3, 1) - rot * (model.mean(1).reshape(3, 1))

        # model_aligned = np.dot(rot, model) #+ trans
        model_aligned = np.dot(rot, model)  + trans
        alignment_error = model_aligned - data

        trans_error = np.sqrt(np.sum(np.multiply(alignment_error, alignment_error), 0)).A[0]

        # pylab.figure(10)
        # pylab.plot(model_aligned.transpose()[:, 0], model_aligned.transpose()[:, 1], 'r', label="aligned test date")
        # pylab.plot(model[0, :], model[1, :], 'g', label="GT")
        # pylab.plot(data[0, :], data[1, :], 'b--', label="original test date")
        # pylab.legend(loc="lower right", shadow=True, fancybox=True)
        #
        # pylab.show()
        #
        # exit(0)
        return rot, trans, trans_error, model_aligned
    else:
        W = np.zeros((2, 2))
        for column in range(model.shape[0]):
            W += np.outer(model_zerocentered[:, column], data_zerocentered[:, column])
        U, d, Vh = np.linalg.linalg.svd(W.transpose())
        S = np.matrix(np.identity(2))
        if (np.linalg.det(U) * np.linalg.det(Vh) < 0):
            S[1, 1] = -1
        rot = U * S * Vh
        print(modelTanspose.shape)
        # print(modelTanspose.mean(1))
        sodfa = rot * (modelTanspose.mean(1).reshape(2,1))
        print("sodfa shape: ")
        print(sodfa.shape)

        trans = dataTranspose.mean(1).reshape(2, 1) - sodfa
        print("trans shape: ")
        print(trans.shape)
         # model_aligned = rot * modelTanspose + trans
        model_aligned = np.dot(rot, modelTanspose)+ trans
        print("model_aligned shape: ")
        print(model_aligned.shape)

        alignment_error = model_aligned - dataTranspose

        trans_error = np.sqrt(np.sum(np.multiply(alignment_error, alignment_error), 0)).A[0]



        return rot, trans, trans_error, model_aligned

#
# def plot_traj(ax, stamps, traj, style, color, label):
#     """
#     Plot a trajectory using matplotlib.
#
#     Input:
#     ax -- the plot
#     stamps -- time stamps (1xn)
#     traj -- trajectory (3xn)
#     style -- line style
#     color -- line color
#     label -- plot legend
#
#     """
#     stamps.sort()
#     interval = numpy.median([s - t for s, t in zip(stamps[1:], stamps[:-1])])
#     x = []
#     y = []
#     last = stamps[0]
#     for i in range(len(stamps)):
#         if stamps[i] - last < 2 * interval:
#             x.append(traj[i][0])
#             y.append(traj[i][1])
#         elif len(x) > 0:
#             ax.plot(x, y, style, color=color, label=label)
#             label = ""
#             x = []
#             y = []
#         last = stamps[i]
#     if len(x) > 0:
#         ax.plot(x, y, style, color=color, label=label)
#
#
# if __name__ == "__main__":
#     # parse command line
#     parser = argparse.ArgumentParser(description='''
#     This script computes the absolute trajectory error from the ground truth trajectory and the estimated trajectory.
#     ''')
#     parser.add_argument('first_file', help='ground truth trajectory (format: timestamp tx ty tz qx qy qz qw)')
#     parser.add_argument('second_file', help='estimated trajectory (format: timestamp tx ty tz qx qy qz qw)')
#     parser.add_argument('--offset', help='time offset added to the timestamps of the second file (default: 0.0)',
#                         default=0.0)
#     parser.add_argument('--scale', help='scaling factor for the second trajectory (default: 1.0)', default=1.0)
#     parser.add_argument('--max_difference',
#                         help='maximally allowed time difference for matching entries (default: 0.02)', default=0.02)
#     parser.add_argument('--save', help='save aligned second trajectory to disk (format: stamp2 x2 y2 z2)')
#     parser.add_argument('--save_associations',
#                         help='save associated first and aligned second trajectory to disk (format: stamp1 x1 y1 z1 stamp2 x2 y2 z2)')
#     parser.add_argument('--plot', help='plot the first and the aligned second trajectory to an image (format: png)')
#     parser.add_argument('--verbose',
#                         help='print all evaluation data (otherwise, only the RMSE absolute translational error in meters after alignment will be printed)',
#                         action='store_true')
#     args = parser.parse_args()
#
#     first_list = associate.read_file_list(args.first_file)
#     second_list = associate.read_file_list(args.second_file)
#
#     matches = associate.associate(first_list, second_list, float(args.offset), float(args.max_difference))
#     if len(matches) < 2:
#         sys.exit(
#             "Couldn't find matching timestamp pairs between groundtruth and estimated trajectory! Did you choose the correct sequence?")
#
#     first_xyz = numpy.matrix([[float(value) for value in first_list[a][0:3]] for a, b in matches]).transpose()
#     second_xyz = numpy.matrix(
#         [[float(value) * float(args.scale) for value in second_list[b][0:3]] for a, b in matches]).transpose()
#     rot, trans, trans_error = align(second_xyz, first_xyz)
#
#     second_xyz_aligned = rot * second_xyz + trans
#
#     first_stamps = first_list.keys()
#     first_stamps.sort()
#     first_xyz_full = numpy.matrix([[float(value) for value in first_list[b][0:3]] for b in first_stamps]).transpose()
#
#     second_stamps = second_list.keys()
#     second_stamps.sort()
#     second_xyz_full = numpy.matrix(
#         [[float(value) * float(args.scale) for value in second_list[b][0:3]] for b in second_stamps]).transpose()
#     second_xyz_full_aligned = rot * second_xyz_full + trans
#
#     if args.verbose:
#         print
#         "compared_pose_pairs %d pairs" % (len(trans_error))
#
#         print
#         "absolute_translational_error.rmse %f m" % numpy.sqrt(numpy.dot(trans_error, trans_error) / len(trans_error))
#         print
#         "absolute_translational_error.mean %f m" % numpy.mean(trans_error)
#         print
#         "absolute_translational_error.median %f m" % numpy.median(trans_error)
#         print
#         "absolute_translational_error.std %f m" % numpy.std(trans_error)
#         print
#         "absolute_translational_error.min %f m" % numpy.min(trans_error)
#         print
#         "absolute_translational_error.max %f m" % numpy.max(trans_error)
#     else:
#         print
#         "%f" % numpy.sqrt(numpy.dot(trans_error, trans_error) / len(trans_error))
#
#     if args.save_associations:
#         file = open(args.save_associations, "w")
#         file.write("\n".join(
#             ["%f %f %f %f %f %f %f %f" % (a, x1, y1, z1, b, x2, y2, z2) for (a, b), (x1, y1, z1), (x2, y2, z2) in
#              zip(matches, first_xyz.transpose().A, second_xyz_aligned.transpose().A)]))
#         file.close()
#
#     if args.save:
#         file = open(args.save, "w")
#         file.write("\n".join(["%f " % stamp + " ".join(["%f" % d for d in line]) for stamp, line in
#                               zip(second_stamps, second_xyz_full_aligned.transpose().A)]))
#         file.close()
#
#     if args.plot:
#         import matplotlib
#
#         matplotlib.use('Agg')
#         import matplotlib.pyplot as plt
#         import matplotlib.pylab as pylab
#         from matplotlib.patches import Ellipse
#
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#         plot_traj(ax, first_stamps, first_xyz_full.transpose().A, '-', "black", "ground truth")
#         plot_traj(ax, second_stamps, second_xyz_full_aligned.transpose().A, '-', "blue", "estimated")
#
#         label = "difference"
#         for (a, b), (x1, y1, z1), (x2, y2, z2) in zip(matches, first_xyz.transpose().A,
#                                                       second_xyz_aligned.transpose().A):
#             ax.plot([x1, x2], [y1, y2], '-', color="red", label=label)
#             label = ""
#
#         ax.legend()
#
#         ax.set_xlabel('x [m]')
#         ax.set_ylabel('y [m]')
#         plt.savefig(args.plot, dpi=90)
#
