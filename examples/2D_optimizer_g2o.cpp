/*
 * test.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: yasir
 */

#include "include/RRR.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <sys/time.h>
#include <assert.h>

#include "g2o/types/slam2d/edge_se2.h"
#include "g2o/types/slam2d/vertex_se2.h"
#include "g2o/types/slam3d/edge_se3.h"
#include "g2o/types/slam3d/vertex_se3.h"
#include "g2o/core/eigen_types.h"
#include "g2o/solvers/cholmod/linear_solver_cholmod.h"
#include "g2o/core/optimization_algorithm_levenberg.h"

#include "g2o/core/sparse_optimizer.h"

#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"


#include <string>


typedef RRR < G2O_Interface
				<
				g2o::VertexSE2, g2o::EdgeSE2
				>
			>
			RRR_2D_G2O;

typedef RRR < G2O_Interface
				<
				g2o::VertexSE3, g2o::EdgeSE3
				>
			>
			RRR_3D_G2O;

// #include "g2o/core/optimization_algorithm_levenberg.h"
// // #include "g2o/core/optimization_algorithm_gauss_newton.h"

// // #include "g2o/solvers/csparse/linear_solver_csparse.h"
// #include "g2o/solvers/cholmod/linear_solver_cholmod.h"

using namespace g2o;

/*!
 * This example show how to apply RRR to already existing instance of g2o optimizer.
 * To simulate this, a graph is read into the optimizer and then a pointer of this
 * optimizer is passed to an instance of RRR.
 *
 *
 */
using namespace Eigen;

int main(int argc, char** argv)
{

	if(argc < 2)
	{
		std::cerr<<"Please specify a graph file to read " <<std::endl;
		std::cerr<<argv[0]<<" graph_file clusteringThreshold[default=50]"<<std::endl;
		return -1;
	}

	// std::cin.get();
	// assert(5 == 7);

	//以写入和在文件末尾添加的方式打开.txt文件，没有的话就创建该文件。
    ofstream time_out;
    time_out.open("/home/zihao/Documents/robust_slam/timeConsumingEachPart/timeComputation.txt", std::ios::out | std::ios::app);  
       if (!time_out.is_open())
        return 0;

	double clusteringThreshold = 1;
	int nIter = 4;
	int num_all_vertex = 0;
	string  fullName;

	struct timeval t1, t2, clusterStart, clusterEnd;

	if(argc == 3)
	{
	  clusteringThreshold = atof(argv[2]);
	  // cout<<"clusteringThreshold: "<<clusteringThreshold<<endl;
	}


	std::cout<<std::endl;
	std::cout<<"---------------------------------------------"<<std::endl;
	std::cout<<"  Robust SLAM based on Consistent clustering "<<std::endl;
	std::cout<<"   Zihao Xu, Cesar Cadena , Roland Siegwart  "<<std::endl;
	std::cout<<"                ASL, ETH 2018                "<<std::endl;
	std::cout<<"---------------------------------------------"<<std::endl;
	std::cout<<std::endl;
	
	
	std::cout<<"Reading from file :"<<argv[1]<<std::endl;
	if(argc == 3)
	{
		std::cout<<"Clustering threshold (t_g) :" << clusteringThreshold<<std::endl;
	}

	fullName = argv[1];
	int posfile = fullName.find_last_of('/');
	string path = fullName.substr(0, posfile+1); 
	cout<<"path: "<<path<<endl;
	// exit(0);
	string  fileName=fullName.substr(posfile+1);
		cout<<"fileName: "<<fileName<<endl;

    fileName = fileName.substr(0, fileName.length()-4);
	cout<<"fileName"<<fileName <<endl;

	string g2ofile,resultfile;

	if(argc == 3)
	{
		g2ofile = path+"resolved_"+fileName+'_'+to_string(int(clusteringThreshold*100))+".g2o";
		resultfile = path+"result_"+fileName+'_'+to_string(int(clusteringThreshold*100))+".g2o";
	}
	else
	{
		g2ofile = path+"resolved_"+fileName+'_'+".g2o";
		resultfile = path+"result_"+fileName+'_'+".g2o";
	}

    cout<<"g2ofile:"<<g2ofile<<endl;
    cout<<"path+fileName: "<<resultfile<<endl;

	
    const char *g2f=g2ofile.data();
    const char *resultFile=resultfile.data();

	// exit(0);
	/* Allocate a g2o optmizer :
	 * It is important that the solver be define in case we want to
	 * reason on a graph already in memory.
	 * */
    
// typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  		SlamBlockSolver;
// typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> 	SlamLinearSolver;
// 	g2o::SparseOptimizer optimizer;
// 	SlamLinearSolver* linearSolver = new SlamLinearSolver();
// 	linearSolver->setBlockOrdering(false);
// 	SlamBlockSolver* blockSolver = new SlamBlockSolver(std::unique_ptr<SlamLinearSolver>(linearSolver));
// 	g2o::OptimizationAlgorithmGaussNewton* solver   = 
// 		new g2o::OptimizationAlgorithmGaussNewton(std::unique_ptr<SlamBlockSolver>(blockSolver));





  typedef BlockSolver< BlockSolverTraits<-1, -1> >  SlamBlockSolver;
  typedef LinearSolverCholmod<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;
  // typedef LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;
   // allocating the optimizer
  SparseOptimizer optimizer;
  auto linearSolver = g2o::make_unique<SlamLinearSolver>();
  linearSolver->setBlockOrdering(false);
  OptimizationAlgorithmGaussNewton* solver = new OptimizationAlgorithmGaussNewton(
    g2o::make_unique<SlamBlockSolver>(std::move(linearSolver)));




 //    g2o::SparseOptimizer optimizer;
	// std::unique_ptr<g2o::BlockSolverTraits<-1, -1>::LinearSolverType> linearSolver;
	// linearSolver = g2o::make_unique<g2o::LinearSolverCholmod<SlamBlockSolver::PoseMatrixType>>();
	// g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton(
	//     g2o::make_unique<g2o::BlockSolver_6_3>(std::move(linearSolver)));	


	// g2o::SparseOptimizer optimizer;
	// std::unique_ptr<g2o::BlockSolver_6_3::LinearSolverType> linearSolver;
	// linearSolver = g2o::make_unique<g2o::LinearSolverCholmod<g2o::BlockSolver_6_3::PoseMatrixType>>();
	// g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(
	//     g2o::make_unique<g2o::BlockSolver_6_3>(std::move(linearSolver)));	


	// OptimizationAlgorithmLevenberg * solver;
 //    solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<g2o::BlockSolverX>(
	// 			g2o::make_unique<g2o::LinearSolverCholmod<g2o::BlockSolverX::PoseMatrixType>>()));	

  	// g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(
   //  	g2o::make_unique<BlockSolverX>(g2o::make_unique<LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>>()));



	optimizer.setAlgorithm(solver);
	// optimizer.setLevenberg(true);


	// SparseOptimizer optimizer;
 //    OptimizationAlgorithmLevenberg * solver;
 //    solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<g2o::BlockSolverX>(
 //    				g2o::make_unique<g2o::LinearSolverCholmod<g2o::BlockSolverX::PoseMatrixType>>()));
 //    optimizer.setAlgorithm(solver);	

	/* load the graph file in the optimizer */
	optimizer.load(argv[1]);


	/* Initialized RRR with the parameters defined */
	RRR_2D_G2O rrr(clusteringThreshold, nIter, fullName);
	// RRR_2D_G2O rrr( nIter, fullName);

	/* Pass the current optimizer pointer to rrr */
// struct timeval t1, t2, clusterStart, clusterEnd;
	gettimeofday(&clusterStart, NULL);

	rrr.setOptimizer(&optimizer,argv[1]);
	// sleep(1);

	gettimeofday(&clusterEnd, NULL);


	/* Find loop closure that are consistent
	 * If the function is passed a bool variable with value true,
	 * it will automatically elimiate all the wrong loops from the
	 * original optimizer. Otherwise, the function removeIncorrectLoops()
	 * can be called to do the same.
	 */

// struct timeval t1, t2, clusterStart, clusterEnd;
	// gettimeofday(&clusterStart, NULL);

	// rrr.robustify();
	rrr.robustify_simplify();  
	// rrr.robustify_simplify_large_first();

	gettimeofday(&t2, NULL);
	//那么函数f运行所花的时间为
	double deltaT = (t2.tv_sec-clusterStart.tv_sec)  + (t2.tv_usec-clusterStart.tv_usec)/1000000.0;// 微秒
	cout<<" "<<endl;
	cout<<"total time consumed: "<<deltaT<<endl;
	cout<<" "<<endl;

	double deltaT_clustering = (clusterEnd.tv_sec-clusterStart.tv_sec) + (clusterEnd.tv_usec-clusterStart.tv_usec) / 1000000.0;// 微秒

	cout<<"clustering time consumed: "<<deltaT_clustering<<endl;

	time_out << fileName << endl; 
    // time_out << "clustering_time: "<< deltaT_clustering <<" intra_time: "<< rrr.deltaT <<" total_time: " <<deltaT <<"\n";  　
    time_out << "clustering_time: "<< deltaT_clustering <<" intra_time: "<< rrr.deltaT/ 1000000.0 <<
    	" inter_time: "<< deltaT - deltaT_clustering - rrr.deltaT/ 1000000.0<<" total_time: " <<deltaT <<"\n";

    time_out << "all_test_closure "<< rrr.allTest <<" true_closure "<< rrr.allTrue<<
    	" accepted_false " <<rrr.allfalseAccepted <<" lost_good "<< rrr.lostGood<<"\n";
    time_out.close();
	/**
	 * If we didn't remove the wrong loop closures earlier, remove them now
	 */
	rrr.removeIncorrectLoops();

	/*
	 * This will write a graph with only the correct loop closures, even when we have not
	 * elimiated incorrect ones using one of the methods above.
	 * */
	std::cout<<std::endl;
	std::cout<<"Output written to flexible-solved.g2o"<<std::endl;

	 // rrr.write(g2f);
	rrr.write("flexible-solved.g2o");

	// rrr.write_resolved_result(g2f, resultFile);

	// rrr.saveOptimizedResult(OPre);

	/**
	 * Since we have removed the incorrect ones, this file would be the same
	 * as the one above.
	 */
	
	//optimizer.save("g2o-saved.g2o");


	return 0;
}





