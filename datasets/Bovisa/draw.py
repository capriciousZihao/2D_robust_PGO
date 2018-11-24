# -*- coding: utf-8 -*-

import time
import pylab
import numpy as np
import sys
import math
sys.path.append('/zihao/.PyCharmCE2018.1/config/scratches/')#引入的模块路经

import scratch_1#应用计算误差模块
from scratch_1 import calErr42Trajectory, align #引用模块中的计算误差函数

#
# VERTEX_SE2 2045 76.6810000 35.4100000 -2.6720000
# VERTEX_SE2 2046 76.6480000 35.3940000 -2.7140000
# VERTEX_SE2 2047 76.6140000 35.3790000 -2.7560000
# VERTEX_SE2 2048 76.5800000 35.3660000 -2.7980000
# VERTEX_SE2 2049 76.5450000 35.3550000 -2.8410000
# VERTEX_SE2 2050 76.5090000 35.3450000 -2.8800000
#
# EDGE_SE2 2045 2046 0.0366682 -0.0006652 -0.0420000 4000 0 0 4000 0 40000
# EDGE_SE2 2046 2047 0.0371591 -0.0004497 -0.0420000 4000 0 0 4000 0 40000
# EDGE_SE2 2047 2048 0.0363930 -0.0007422 -0.0420000 4000 0 0 4000 0 40000
# EDGE_SE2 2048 2049 0.0366599 -0.0014335 -0.0430000 4000 0 0 4000 0 40000
# EDGE_SE2 2049 2050 0.0373467 -0.0011075 -0.0390000 4000 0 0 4000 0 40000



# theta = -0.0420
# xa = 0.0366682
# ya = -0.0006652
#
# thetab = -2.6720000
# xb = 76.6810000
# yb = 35.4100000
#
# a = np.zeros((3, 3))
# b = np.zeros((3, 1))
# a = [[np.cos(theta), -np.sin(theta), xa],
#      [np.sin(theta), np.cos(theta), ya],
#      [0, 0, 1]]
# print(a)
#
# b = [[xb], [yb], [thetab]]
# print(b)
# c = np.dot(a, b)
# print(c)
# exit(0)

# v = np.zeros((3, 3))
# e = np.zeros((3, 3))
# print (v)
# vertex = [76.681, 35.41, -2.672]
# v[0, 0] = np.cos(vertex[2])
# v[0, 1] = -np.sin(vertex[2])
# v[0, 2] = vertex[0]
# v[1, 0] = np.sin(vertex[2])
# v[1, 1] = np.cos(vertex[2])
# v[1, 2] = vertex[1]
# v[2, 0] = 0
# v[2, 1] = 0
# v[2, 2] = 1
# edge = [0.0366682, -0.0006652, -0.042]
# e[0, 0] = np.cos(edge[2])
# e[0, 1] = -np.sin(edge[2])
# e[0, 2] = edge[0]
# e[1, 0] = np.sin(edge[2])
# e[1, 1] = np.cos(edge[2])
# e[1, 2] = edge[1]
# e[2, 0] = 0
# e[2, 1] = 0
# e[2, 2] = 1
#
# update = np.dot(v, e)
# print(update)
# print(update[0, 2], update[1, 2])
#
# vertex = [76.648, 35.394, -2.714]
# v[0, 0] = np.cos(vertex[2])
# v[0, 1] = -np.sin(vertex[2])
# v[0, 2] = vertex[0]
# v[1, 0] = np.cos(vertex[2])
# v[1, 1] = np.cos(vertex[2])
# v[1, 2] = vertex[1]
# v[2, 0] = 0
# v[2, 1] = 0
# v[2, 2] = 1
# print(v)
#
# exit(0)





def plotData(X, y):
    length = len(y)

    pylab.figure(1)

    pylab.plot(X, y, 'rx')
    pylab.xlabel('meters')
    pylab.ylabel('meters')

    pylab.show()  # 让绘制的图像在屏幕上显示出来

    # (X,y) = loadData('Bovisa04_GT_GPS.g2o')  #Bovisa04_GT_GPS.g2o

##########打开 基准 文件############################
pathGT="/home/zihao/cprogram/data/sorted/bovisa/"
#nameGT="city10000_groundtruth_with_initialGuess.g2o"
#nameGT="Bovisa06_GT_GPS.g2o"
nameGT="Bovisa04_GT_GPS.g2o"
fullNameGT=pathGT+nameGT
GTfile= open(fullNameGT, 'r')      # 以只读方式打开基准值文件
###############################################




#########打开 测试 文件##############################
path="/home/zihao/cprogram/code/robustSLAM/robuts_slam_trans_dis_2/build/examples/"
#path="/home/zihao/cprogram/rrr/build/examples/result/ringCity/"
#path = "/home/zihao/cprogram/data/sorted/Bovisa/result/DSC/"
name="flexible-solved.g2o"

fullName=path+name
#inFile = open("/home/zihao/cprogram/data/sorted/city10000/result/DSC/city10000_DSC_0_5.g2o", 'r')  # 以只读方式打开某fileName文件    Bovisa04_GT_GPS.g2o
inFile = open(fullName, 'r')
################################################

#########求出待测数据文件中的顶点，边，回环的个数########
countVertex = 0
countEdge=0
countLoopClo=0
for line in inFile:
    trainingSet = line.split(' ')         # 对于每一行，按','把数据分开，这里是分成两部分
    if trainingSet[0] == "VERTEX_SE2":    # 求顶点的个数
        countVertex = countVertex + 1
    elif trainingSet[0] == "EDGE_SE2":    # 求边的个数
        if abs(eval(trainingSet[1])-eval((trainingSet[2])))!=1: # 求回环的个数
            countLoopClo = countLoopClo + 1
        else:
            countEdge = countEdge + 1     # 普通边的个数
edgeMatrix = np.zeros((countEdge, 5))        # 定义存贮 边   的矩阵
vertexMatrix = np.zeros((countVertex, 4))    # 定义存定 点   的矩阵
loopCloMatrix = np.zeros((countLoopClo, 5))  # 定义存定 回环  的矩阵
###################################################

# ######求出基准数据中的顶点个数 and the number of true loop closures，然后将信息放到它的矩阵中############
countVertexGT = 0
countLcEdgeGT = 0
for line in GTfile:
    trainingSet = line.split(' ')         # 对于每一行，按','把数据分开，这里是分成两部分
    if trainingSet[0] == "VERTEX_SE2" or trainingSet[0] == "VERTEX_SE3":    # 求顶点的个数
        countVertexGT = countVertexGT + 1
    elif trainingSet[0].find("EDGE_SE2") != -1 or trainingSet[0].find("EDGE_SE3") != -1:
        if abs(eval(trainingSet[1])-eval((trainingSet[2]))) != 1: # 求回环的个数
            countLcEdgeGT = countLcEdgeGT + 1
# edgeMatrix = np.zeros((countEdge, 5))        # 定义存贮 边   的矩阵
GTvertexMatrix = np.zeros((countVertexGT, 4))    # 定义存定 GT点   的矩阵
GT_Lc_Matrix = np.zeros((countLcEdgeGT, 2))   #store real loop closure

# loopCloMatrix = np.zeros((countLoopClo, 5))  # 定义存定 回环  的矩阵
GTfile.seek(0)
countVertex = 0
countLcEdgeGT = 0
for line in GTfile:
    trainingSet = line.split(' ')
    if trainingSet[0] == "VERTEX_SE2" or trainingSet[0] == "VERTEX_SE3":    # 存储顶点信息
        GTvertexMatrix[countVertex, :] = [eval(trainingSet[1]), eval(trainingSet[2]), eval(trainingSet[3]), eval(trainingSet[4])]
        countVertex = countVertex + 1
    elif trainingSet[0].find("EDGE") != -1:  # check if this line is edge or not
        if abs(eval(trainingSet[1])-eval(trainingSet[2])) > 1.1:              # 存储loop closure信息
            GT_Lc_Matrix[countLcEdgeGT, :] = [trainingSet[1], trainingSet[2]]
            countLcEdgeGT = countLcEdgeGT+1
if countLcEdgeGT != GT_Lc_Matrix.shape[0]:
    print("number of stored real loop closure is smaller than it shoud be")
    time.sleep(4)
############################################################


############将测试文件中的点、边、回环分别存储到矩阵中###############
countVertex = 0
countEdge = 0
countLoopClo = 0
inFile.seek(0)
for line in inFile:
    trainingSet = line.split(' ')         # 对于每一行，按','把数据分开，这里是分成两部分
    if trainingSet[0] == "VERTEX_SE2":    # 存储顶点信息
        #print(trainingSet[0])
        #print(eval(trainingSet[1]))
        #print(type(line))
        vertexMatrix[countVertex, :] = [eval(trainingSet[1]), eval(trainingSet[2]), eval(trainingSet[3]), eval(trainingSet[4])]
        countVertex = countVertex + 1
    elif trainingSet[0] == "EDGE_SE2":                          # 如果是边
        if abs(eval(trainingSet[1])-eval((trainingSet[2])))>1.1: # 如果是回环边
            loopCloMatrix[countLoopClo, :] = [eval(trainingSet[1]), eval(trainingSet[2]), eval(trainingSet[3]),
                                            eval(trainingSet[4]), eval(trainingSet[5])]
            countLoopClo = countLoopClo + 1
#        else:                                                   # 普通里程计边
#            vertexMatrix[countVertex, :] = [eval(trainingSet[1]), eval(trainingSet[2]), eval(trainingSet[3]),
#                                            eval(trainingSet[4])]
    else:
        if trainingSet[0] == "EDGE_SE2_MIXTURE":   # 暂时无法处理混合边
            continue
        print("%s 不是EDGE_SE2也不是VERTEX_SE2" % trainingSet )
        #input("请输入：")
        time.sleep(3)

print (countLoopClo)
time.sleep(3)
#####################对于基准值和测试值数据个数不一致的情况，提取测试数据#########################################
timesTT = vertexMatrix.shape[0] / float(GTvertexMatrix.shape[0])
# print("i: %f" % i )
# timesTT = math.floor(i)
# print("timesTT: %d" % timesTT )
if timesTT < 1:
    print("truthdata中的数据更多，无法处理")
    exit(0)
else:
    testCollectData = np.zeros((GTvertexMatrix.shape[0], 2))
    for l in range(GTvertexMatrix.shape[0]):
        testCollectData[l, :] = vertexMatrix[math.floor(l*timesTT), 1:3]
print(testCollectData[0:20, 0])
pylab.figure(4)
pylab.plot(testCollectData[:, 0], testCollectData[:, 1], 'r*',label="test date")
# pylab.show()
# exit(0)
########比较基准文件中的点数目与测试文件中的点的数目###################

if countVertex != countVertexGT:
    print("countVertex: %d   countVertexGT: %d" % (countVertex, countVertexGT))
    print("%d vertexes in GT while %d vertexes in test file" % (countVertexGT, countVertex))
    time.sleep(3)
    sameVertex = 0
    if countVertex >= countVertexGT:
        ErrMatrix = np.zeros((countVertexGT, 4))  # store position error of vertexes
    elif countVertex < countVertexGT:
        print("GT(%d) file contains more vertex than test file(%d)" % (countVertexGT, countVertex))
        time.sleep(3)
else:
    ErrMatrix = np.zeros((countVertexGT, 4))

    #pylab.figure(1)
    #pylab.plot(vertexMatrix[:, 1], vertexMatrix[:, 2], 'r',label="test date")
    #pylab.plot(GTvertexMatrix[:, 1], GTvertexMatrix[:, 2], 'g',label="GT data")
    #pylab.legend(loc="lower right",shadow=True, fancybox=True)
    #pylab.show()
    #str = input("请输入：")
    #exit(0)
    #time.sleep(3)
###############################################################

################计算点的ATE######################################
backVertexMatrix = np.zeros(vertexMatrix.shape)
backVertexMatrix = vertexMatrix
calErr42Trajectory(GTvertexMatrix, vertexMatrix, ErrMatrix)
mn1, mn2, mn3, mn4 = align(GTvertexMatrix[:, 1:3], testCollectData)

###############################################################


fig1 = pylab.figure(1)
ax1 = fig1.add_subplot(1, 2, 1)
ax1.plot(vertexMatrix[:, 1], vertexMatrix[:, 2], 'b')   # 绘制所有的Vertex

LC_Valid = np.zeros((loopCloMatrix.shape[0], 1))
if loopCloMatrix.shape[0]:
    print(loopCloMatrix[0, :])
else:
    print("no loop closure")
    time.sleep(3)

#### 判断回环的真假
for i in range(loopCloMatrix.shape[0]):
    if (vertexMatrix[loopCloMatrix[i, 0], 0] != loopCloMatrix[i, 0]) or (vertexMatrix[loopCloMatrix[i, 1], 0] != loopCloMatrix[i, 1]):
        print("读取回环对应顶点信息有误！！！")
        sys.exit(0)
    for j in range(GT_Lc_Matrix.shape[0]):  # check the loopClosure is True or not, true LC is green, wrong LC is red
        if GT_Lc_Matrix[j, 0] == loopCloMatrix[i, 0] and GT_Lc_Matrix[j, 1] == loopCloMatrix[i, 1]:
            LC_Valid[i] = 1
            break
        elif GT_Lc_Matrix[j, 0] == loopCloMatrix[i, 1] and GT_Lc_Matrix[j, 1] == loopCloMatrix[i, 0]:
            LC_Valid[i] = 1
            break
        else:
            LC_Valid[i] = 0
    if i%100==0:
        print("i=%d" % i)

#### 画回环，真回环用绿色，假的用红色
for  i in range(loopCloMatrix.shape[0]):
    x = [vertexMatrix[loopCloMatrix[i, 0], 1], vertexMatrix[loopCloMatrix[i, 1], 1]]
    y = [vertexMatrix[loopCloMatrix[i, 0], 2], vertexMatrix[loopCloMatrix[i, 1], 2]]
    if LC_Valid[i]:
        ax1.plot(x, y, 'g--')
    else:
        ax1.plot(x, y, 'r--')
ax1.set_title(name[:-4])


# pylab.show()#让绘制的图像在屏幕上显示出来
ax = fig1.add_subplot(1, 2, 2)

#ax.plot(range(ErrMatrix.shape[0]), ErrMatrix[:, 3], vertexMatrix[:, 1], 'g')
ax.plot(range(ErrMatrix.shape[0]), ErrMatrix[:, 3], 'g')
meanErr = np.mean(ErrMatrix[:, 3])

#pylab.plot(vertexMatrix[:, 1], vertexMatrix[:, 2], 'r', label="test date")
#pylab.plot(GTvertexMatrix[:, 1], GTvertexMatrix[:, 2], 'g', label="GT data")
#pylab.legend(loc="lower right", shadow=True, fancybox=True)

tit = 'error of each vertex, mean error: '
ax.set_title(tit+str(int(meanErr*100)/100.0))
pylab.savefig(path+"%s.png" % name)

pylab.figure(2)
pylab.plot(vertexMatrix[:, 1], vertexMatrix[:, 2], 'r--', label="transformed test date")
pylab.plot(GTvertexMatrix[:, 1], GTvertexMatrix[:, 2], 'g--', label="GT")
pylab.plot(testCollectData[:, 0], testCollectData[:, 1], 'b*', label="original test date")

pylab.legend(loc="lower right", shadow=True, fancybox=True)

pylab.savefig(path+"%s_.png" % name)

# mn4
pylab.figure(3)
print(mn4.shape)
print(mn4[1, :])
pylab.plot(mn4.transpose()[:, 0], mn4.transpose()[:, 1], 'r', label="aligned test date")
pylab.plot(GTvertexMatrix[:, 1], GTvertexMatrix[:, 2], 'g', label="GT")
pylab.plot(testCollectData[:, 0], testCollectData[:, 1], 'b--', label="original test date")
pylab.legend(loc="lower right", shadow=True, fancybox=True)
print(mn1)
print(mn2)
print("误差均值：%f" % mn3.mean())
pylab.show()  # 让绘制的图像在屏幕上显示出来
# plotData(X,y)
