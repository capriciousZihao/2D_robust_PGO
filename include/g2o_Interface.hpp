// RRR - Robust Loop Closing over Time
// Copyright (C) 2014 Y.Latif, C.Cadena, J.Neira
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef BACKEND_G2O_HPP_
#define BACKEND_G2O_HPP_

#include "backEndInterface.hpp"

#include "g2o/types/slam2d/vertex_se2.h"

#include "g2o/types/slam3d/vertex_se3.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/robust_kernel.h"
#include <g2o/core/robust_kernel_impl.h>
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include <Eigen/Geometry>
#include <numeric>
#include <stdio.h>

#include<iostream>

using namespace std;

typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  		SlamBlockSolver;
typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> 	SlamLinearSolver;

// typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  		SlamBlockSolver;
// typedef g2o::LinearSolverCholmod<SlamBlockSolver::PoseMatrixType> 	SlamLinearSolver;

template<typename VertexType, typename EdgeType>
class G2O_Interface : public BaseBackend
{
	typedef
			std::map< IntPair, g2o::HyperGraph::Edge* >
			IntPairToEdgePtrMap;

	typedef
			std::vector< g2o::HyperGraph::Edge* >
			EdgePtrVector;

	typedef g2o::SparseOptimizer OptimizerType;

	EdgePtrVector 					odomteryEdges;
	IntPairToEdgePtrMap				loopclosureEdges;

	
	g2o::OptimizableGraph::EdgeSet 	activeEdges, odoEdgeSegment;
	g2o::OptimizableGraph::VertexSet 	activeVertexes;
	// Transform<double, 3, 1>         try_to_get_vertex_position;
	Eigen::Matrix<double, 3, 3>     M33d;


	bool initialized;

public:
	g2o::SparseOptimizer*			optimizer;
	int odoSize;
	G2O_Interface()
	{
		optimizer = NULL;
		initialized = false;
	}

	bool setOptimizer(void* opt)
	{
		if(optimizer==NULL)
		{
			optimizer = (OptimizerType*)opt;
			initialized = true;
			store();
		}
		else
		{
			std::cerr<<"Already existing optimizer?"<<std::endl;
			return false;
		}
		return true;
	}


	bool getLoopClosures(IntPairSet& loops)
	{
		if(optimizer == NULL)
		{
			std::cerr<<"Please read in a g2o file or pass the pointer to an existing optimizer before calling getLoopClosures()"<<std::endl;
			return false;
		}


		g2o::OptimizableGraph::EdgeSet::iterator
		eIt = optimizer->edges().begin(),
		eEnd = optimizer->edges().end();

		for( ; eIt!=eEnd ; eIt++)
		{
			int e1 = (*eIt)->vertices()[0]->id();
			int e2 = (*eIt)->vertices()[1]->id();
			if(std::abs(e1-e2) > 1)
			{
				// if(e1 == 2480)
				// {
				// 	cout<<"find the desired node!"<<endl;
				// 	exit(0);
				// }
				loops.insert(IntPair(e1,e2));
				loopclosureEdges[IntPair(e1,e2)] = *eIt;

			}
			else
			{
				odomteryEdges.push_back(*eIt);
			}
		}
		odoSize = odomteryEdges.size();
		// cout<<"loop size: "<<loops.size()<<endl;
		// printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
		// exit(0);

		for(int i=0; i<odomteryEdges.size(); i++)
		{
			if(odomteryEdges[i]->vertex(0)->id() != i or odomteryEdges[i]->vertex(1)->id() != (i+1))
			{
				cout<<"check if it is odometry edge fail"<<endl;
				cout<<"i is "<<i<<" but first vertex id is "<< 
				 	odomteryEdges[i]->vertex(0)->id()<<" and second vertex id is "<<odomteryEdges[i]->vertex(1)->id()<<endl;
				exit(0);
			}
		}

		return true;
	}

	virtual ~G2O_Interface(){}

	int vertexCount() { return optimizer->vertices().size(); };
	int edgeCount() { return optimizer->edges().size(); };
	


	bool store() // store the current graph state
	{
		g2o::OptimizableGraph::VertexIDMap::iterator
			vIt = optimizer->vertices().begin(),
			vEnd = optimizer->vertices().end();

		for(; vIt!=vEnd ; vIt++ )	static_cast< VertexType* >(vIt->second)->push();

		return true;
	}

	bool restore() // restore from backup the previous estimate of the graph
	{
		g2o::OptimizableGraph::VertexIDMap::iterator
			vIt = optimizer->vertices().begin(),
			vEnd = optimizer->vertices().end();

		for(; vIt!=vEnd ; vIt++ )	static_cast< VertexType* >(vIt->second)->pop();

		// HACK : store again!
		store();
		return true;

	}

	bool read(
			const char* filename
			)
	{
		if(optimizer != NULL)
		{
			std::cerr<<"An allocated optimizer already exists "<<std::endl;
			return false;
		}
		optimizer = new g2o::SparseOptimizer;
		//SlamLinearSolver* linearSolver = new SlamLinearSolver();
		std::unique_ptr<SlamLinearSolver> linearSolver (new SlamLinearSolver());
		linearSolver->setBlockOrdering(false);
		//SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
		std::unique_ptr<SlamBlockSolver> blockSolver ( new SlamBlockSolver (std::move(linearSolver)));
		//g2o::OptimizationAlgorithmGaussNewton* solverGauss   = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);
		g2o::OptimizationAlgorithmGaussNewton* solverGauss   = new g2o::OptimizationAlgorithmGaussNewton(std::move(blockSolver));
		optimizer->setAlgorithm(solverGauss);


		if(!optimizer->load(filename)){
			std::cerr<<"Can't find file to read : "<<filename<<std::endl;;
			return false;
		}
		store();
		initialized = true;
		return true;
	}

	// Optimize on a given set of loop closure links
	// @return: The error of each link + overall error of the active part of the graph as the last element of the vector

	// IntPairDoubleMap optimize(const IntPairSet& activeLoops, const int nIterations,
	// 	 std::vector<std::vector<std::pair<std::pair<int, int>, double> > > & odoEdgeRelate2LC_Error)
	IntPairDoubleMap optimize(const IntPairSet& activeLoops, const int nIterations,
		std::vector<std::pair<std::pair<int, int>, double> > & odoEdgeRelateLC_Error,
		std::vector<std::array<double, 5>> & odoEdgeError)
	{

	// EdgePtrVector 					odomteryEdges;
	// IntPairToEdgePtrMap				loopclosureEdges;
	
	// g2o::OptimizableGraph::EdgeSet 	activeEdges;
	// g2o::OptimizableGraph::VertexSet 	activeVertexes;

		std::pair<std::pair<int, int>, double> odoErrorElement;
		int node, endNode;
		std::array<double, 5> element_error;
		double information_Angle;
		
		odoEdgeRelateLC_Error.clear();
		odoEdgeError.clear();

		activeEdges.clear();
		activeEdges.insert(odomteryEdges.begin(),odomteryEdges.end());

		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			activeEdges.insert( loopclosureEdges[*it]);
		}

		IntPairDoubleMap loopClosureLinkError;
		restore();
		//
		optimizer->setVerbose(false);
	    //optimizer->setVerbose(true);
		optimizer->findGauge()->setFixed(true);
		///////////////////
		//std::cerr<<"optimizer->findGauge() output: "<<(*(optimizer->findGauge()))<<std::endl;;
		//optimizer->vertex(0)->setFixed(true);
		optimizer->initializeOptimization(activeEdges);
		optimizer->optimize(nIterations,false);
		optimizer->computeActiveErrors();

		// for(int i = 0; i < odoSize; i++)
		// {
		// 	// odoEdgeError.push_back(dynamic_cast< EdgeType* >(odomteryEdges[i])->chi2());
		// 	element_error[0]  = (dynamic_cast< EdgeType* >(odomteryEdges[i])->chi2());
		// 	element_error[1]  = (dynamic_cast< EdgeType* >(odomteryEdges[i])->error())(0);
		// 	element_error[2]  = (dynamic_cast< EdgeType* >(odomteryEdges[i])->error())(1);
		// 	element_error[3]  = (dynamic_cast< EdgeType* >(odomteryEdges[i])->error())(2);
		// 	information_Angle = (dynamic_cast< EdgeType* >(odomteryEdges[i])->information())(2,2);
		// 	element_error[4]  = element_error[3]*element_error[3]*information_Angle;
		// 	odoEdgeError.push_back(element_error);
		// }

		double sumLoopChieErr = 0;
		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			//dynamic_cast< EdgeType* >(loopclosureEdges[*it])->computeError();
			loopClosureLinkError[*it] = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();
			sumLoopChieErr = sumLoopChieErr + loopClosureLinkError[*it];
		}

		// NOTE : The number of edges involved is the last element
		loopClosureLinkError[IntPair(-1,0)] = optimizer->activeChi2(); // There can be no links with negative IDS
		loopClosureLinkError[IntPair(-1,-1)] = optimizer->activeEdges().size();

		loopClosureLinkError[IntPair(-2,0)] = sumLoopChieErr; // There can be no links with negative IDS
		loopClosureLinkError[IntPair(-2,-1)] = activeLoops.size();		
		//
		return loopClosureLinkError;

	}

	// IntPairDoubleMap optimize_oneMember(const IntPairSet& activeLoops, const int nIterations,
	// 	EdgePtrVector & manmadeEdges, const std::vector<std::pair<int,int > > & artificialLoops)
	// {
	// 	if(activeLoops.size() == 1)
	// 	{
	// 		std::pair<std::pair<int, int>, double> odoErrorElement;
	// 		int node, endNode;
	// 		std::array<double, 5> element_error;
	// 		double information_Angle;

	// 		activeEdges.clear();
	// 		activeEdges.insert(odomteryEdges.begin(),odomteryEdges.end());

	// 		for(
	// 			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
	// 			it!=end;
	// 			it++)
	// 		{
	// 			activeEdges.insert( loopclosureEdges[*it]);
	// 		}
	// 		for(int i =0; i < manmadeEdges.size(); i++)
	// 		{
	// 			activeEdges.insert( manmadeEdges[i]);
	// 		}

	// 		IntPairDoubleMap loopClosureLinkError;
	// 		restore();
	// 		//
	// 		optimizer->setVerbose(false);
	// 	    //optimizer->setVerbose(true);
	// 		optimizer->findGauge()->setFixed(true);
	// 		///////////////////
	// 		//std::cerr<<"optimizer->findGauge() output: "<<(*(optimizer->findGauge()))<<std::endl;;
	// 		//optimizer->vertex(0)->setFixed(true);
	// 		optimizer->initializeOptimization(activeEdges);

	// 		// optimizer->addEdge(manmadeEdges[0]);
	// 		// optimizer->addEdge(manmadeEdges[1]);
	// 		// optimizer->addEdge(manmadeEdges[2]);
	// 		// optimizer->addEdge(manmadeEdges[3]);

	// 		optimizer->optimize(nIterations,false);
	// 		optimizer->computeActiveErrors();

	// 		double sumLoopChieErr = 0;
	// 		for(
	// 			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
	// 			it!=end;
	// 			it++)
	// 		{
	// 			//dynamic_cast< EdgeType* >(loopclosureEdges[*it])->computeError();
	// 			loopClosureLinkError[*it] = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();
	// 			sumLoopChieErr = sumLoopChieErr + loopClosureLinkError[*it];
	// 			cout<<"loop "<<(*it).first<<" "<<(*it).second<<" has error "<<loopClosureLinkError[*it]<<endl;
	// 		}

	// 		cout<<"artificialLoops size is "<<artificialLoops.size()<<endl;
	// 		for(int i = 0 ; i < manmadeEdges.size(); i++)
	// 		{
	// 			//dynamic_cast< EdgeType* >(loopclosureEdges[*it])->computeError();
	// 			// cout<<"debug 1"<<endl;
	// 			loopClosureLinkError[artificialLoops[i]] = dynamic_cast< EdgeType* >(manmadeEdges[i])->chi2();
	// 			// cout<<"debug 2"<<endl;
	// 			sumLoopChieErr = sumLoopChieErr + loopClosureLinkError[artificialLoops[i]];
	// 			// cout<<"debug 3"<<endl;
	// 			cout<<"loop "<<(artificialLoops[i]).first<<" "<<(artificialLoops[i]).second<<
	// 				" has error "<<loopClosureLinkError[artificialLoops[i]]<<endl;
	// 		}
	// 		// NOTE : The number of edges involved is the last element
	// 		loopClosureLinkError[IntPair(-1,0)] = optimizer->activeChi2(); // There can be no links with negative IDS
	// 		loopClosureLinkError[IntPair(-1,-1)] = optimizer->activeEdges().size();

	// 		loopClosureLinkError[IntPair(-2,0)] = sumLoopChieErr; // There can be no links with negative IDS
	// 		loopClosureLinkError[IntPair(-2,-1)] = activeLoops.size();		
	// 		//
	// 		return loopClosureLinkError;
	// 	}
	// 	else
	// 	{
	// 		cout<<"the cluster size is not 1, but this function is for 1 menber clsuter optimize"<<endl;
	// 		printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
	// 		exit(0);
	// 	}
	// }

	void huber_chi2(double &  e, g2o::Vector3 & rho)
	{
		double  _delta = 1;
	  	number_t dsqr = _delta * _delta;
	  	if (e <= dsqr) { // inlier
	    	rho[0] = e;
	    	// rho[1] = 1.;
	    	// rho[2] = 0.;
	  	} 
	  	else { // outlier
	    	number_t sqrte = sqrt(e); // absolut value of the error
	    	rho[0] = 2*sqrte*_delta - dsqr; // rho(e)   = 2 * delta * e^(1/2) - delta^2
	    	// rho[1] = _delta / sqrte;        // rho'(e)  = delta / sqrt(e)
	    	// rho[2] = - 0.5 * rho[1] / e;    // rho''(e) = -1 / (2*e^(3/2)) = -1/2 * (delta/e) / e
	  	}
	}


	void cauchy_chi2(double &  e2, g2o::Vector3 & rho)
	{
		double  _delta = 1;
	   	number_t dsqr = _delta * _delta;
	 	number_t dsqrReci = 1. / dsqr;
	  	number_t aux = dsqrReci * e2 + 1.0;
	  	rho[0] = dsqr * log(aux);
	  	// rho[1] = 1. / aux;
	  	// rho[2] = -dsqrReci * std::pow(rho[1], 2);
	}	
	// IntPairDoubleMap optimize_active_robustKernel(const IntPairSet& activeLoops, const int nIterations,
	// 	int & startRegion, int & endRegion, std::vector<std::array<double, 5>> & odoEdgeError,
	// 	string ba, string aa)
	IntPairDoubleMap optimize_active_robustKernel(const IntPairSet& activeLoops, const int nIterations,
		int startRegion, int endRegion, std::vector<std::array<double, 5>> & odoEdgeError)	
	{

	// EdgePtrVector 					odomteryEdges;
	// IntPairToEdgePtrMap				loopclosureEdges;
	
	// g2o::OptimizableGraph::EdgeSet 	activeEdges;
	// g2o::OptimizableGraph::VertexSet 	activeVertexes;
		cout<<" optimize_active_robustKernel"<<endl;
		// cout<<"nodes range form "<<startRegion<<" to "<<endRegion<<endl;
		std::pair<std::pair<int, int>, double> odoErrorElement;
		int node, endNode;
		std::array<double, 5> element_error;
		double information_Angle;

		odoEdgeError.clear();

		activeEdges.clear();
		activeEdges.insert(odomteryEdges.begin(),odomteryEdges.end());
		// activeEdges.insert(odomteryEdges.begin() + startRegion,odomteryEdges.begin() + endRegion);
		// for(
		// 	EdgePtrVector::const_iterator it = odomteryEdges.begin(), end = odomteryEdges.end();
		// 	it!=end;
		// 	it++)
		// for(
		// 	int i = 0;
		// 	i < odomteryEdges.size();
		// 	i++)			
		// {
  //        	// 核函数
  //        	// dynamic_cast< EdgeType* >(loopclosureEdges[*it])->setRobustKernel( new g2o::RobustKernelHuber() );	
  //        	dynamic_cast< EdgeType* >(odomteryEdges[i])->setRobustKernel( new g2o::RobustKernelCauchy() );
		// 	activeEdges.insert( odomteryEdges[i]);
		// }


		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			// if(*it == (std::pair<int,int> (4190, 4085)))
			// 	continue;

         	// 核函数
         	// g2o::RobustKernelHuber();
         	// RobustKernelCauchy
         	// dynamic_cast< EdgeType* >(loopclosureEdges[*it])->setRobustKernel( new g2o::RobustKernelHuber() );	
         	dynamic_cast< EdgeType* >(loopclosureEdges[*it])->setRobustKernel( new g2o::RobustKernelCauchy() );

        	// loopclosureEdges[*it]->setRobustKernel( new g2o::RobustKernelHuber() );			
			// loopclosureEdges[*it].
			activeEdges.insert( loopclosureEdges[*it]);
		}
		cout<<"nIterations is "<<nIterations<<endl;

		IntPairDoubleMap loopClosureLinkError;
		restore();
		//
		optimizer->setVerbose(false);
	    // optimizer->setVerbose(true);

        endRegion = 0;
		optimizer->vertex(endRegion)->setFixed(true);

		// optimizer->vertices().clear();

		optimizer->initializeOptimization(activeEdges);

		optimizer->optimize(nIterations,false);
		optimizer->computeActiveErrors();
		// optimizer->vertex(endRegion)->setFixed(false);

		double sumLoopChieErr = 0, sum_cauchy =0.0, sum_huber = 0.0, middle_chi;
		g2o::Vector3 rho;
		if(activeLoops.size() >= 1)
		{
			for(
				IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
				it!=end;
				it++)
			{
				middle_chi = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();	
				sumLoopChieErr = sumLoopChieErr + middle_chi;
				// loopClosureLinkError[*it] = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();	

				cout<<"loop "<<(*it).first<<" "<<(*it).second<<"  has error "<<middle_chi<<endl;

				// cauchy 
				// huber_chi2(middle_chi, rho);
				cauchy_chi2(middle_chi, rho);
				// loopClosureLinkError[*it] = rho[0];
				loopClosureLinkError[*it] = middle_chi;
				sum_cauchy = sum_cauchy + loopClosureLinkError[*it];	

	  			cout<<"*   cauchy    has robust error "<<rho[0]<<"  sum_cauchy:"<<sum_cauchy<<" "<<"  sum_chi:"<<sumLoopChieErr<<endl;
			}

			// NOTE : The number of edges involved is the last element
			// loopClosureLinkError[IntPair(-1,0)] = optimizer->activeChi2(); // There can be no links with negative IDS
			loopClosureLinkError[IntPair(-1,0)] = optimizer->activeRobustChi2();
			loopClosureLinkError[IntPair(-1,-1)] = optimizer->activeEdges().size();

			loopClosureLinkError[IntPair(-2,0)] = sum_cauchy; // There can be no links with negative IDS
			loopClosureLinkError[IntPair(-2,-1)] = activeLoops.size();		
			//
			return loopClosureLinkError;
		}
		else
		{
			for(
				IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
				it!=end;
				it++)
			{
				middle_chi = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();	
				sumLoopChieErr = sumLoopChieErr + middle_chi;
				// loopClosureLinkError[*it] = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();	

				cout<<"loop "<<(*it).first<<" "<<(*it).second<<"  has error "<<middle_chi<<endl;

				// cauchy 
				// huber_chi2(middle_chi, rho);
				cauchy_chi2(middle_chi, rho);
				loopClosureLinkError[*it] = rho[0];
				sum_cauchy = sum_cauchy + loopClosureLinkError[*it];	

	  			cout<<"*   cauchy    has robust error "<<rho[0]<<"  sum_cauchy:"<<sum_cauchy<<" "<<"  sum_chi:"<<sumLoopChieErr<<endl;
			}

			// NOTE : The number of edges involved is the last element
			// loopClosureLinkError[IntPair(-1,0)] = optimizer->activeChi2(); // There can be no links with negative IDS
			loopClosureLinkError[IntPair(-1,0)] = optimizer->activeChi2();
			loopClosureLinkError[IntPair(-1,-1)] = optimizer->activeEdges().size();

			loopClosureLinkError[IntPair(-2,0)] = middle_chi; // There can be no links with negative IDS
			loopClosureLinkError[IntPair(-2,-1)] = activeLoops.size();		
			//
			cout<<"end optimize_active_robustKernel"<<endl;
			return loopClosureLinkError;
		}
	}

	IntPairDoublePairMap optimize_active_robustKernel_pair_return(const IntPairSet& activeLoops, const int nIterations,
		int startRegion, int endRegion, std::vector<std::array<double, 5>> & odoEdgeError)	
	{

	// EdgePtrVector 					odomteryEdges;
	// IntPairToEdgePtrMap				loopclosureEdges;
	
	// g2o::OptimizableGraph::EdgeSet 	activeEdges;
	// g2o::OptimizableGraph::VertexSet 	activeVertexes;
		cout<<" optimize_active_robustKernel_pair_return"<<endl;
		// cout<<"nodes range form "<<startRegion<<" to "<<endRegion<<endl;
		std::pair<std::pair<int, int>, double> odoErrorElement;
		int node, endNode;
		std::array<double, 5> element_error;
		double information_Angle;
		DoublePair error_pair; 

		odoEdgeError.clear();

		activeEdges.clear();
		activeEdges.insert(odomteryEdges.begin(),odomteryEdges.end());
		// activeEdges.insert(odomteryEdges.begin() + startRegion,odomteryEdges.begin() + endRegion);
		// for(
		// 	EdgePtrVector::const_iterator it = odomteryEdges.begin(), end = odomteryEdges.end();
		// 	it!=end;
		// 	it++)
		// for(
		// 	int i = 0;
		// 	i < odomteryEdges.size();
		// 	i++)			
		// {
  //        	// 核函数
  //        	// dynamic_cast< EdgeType* >(loopclosureEdges[*it])->setRobustKernel( new g2o::RobustKernelHuber() );	
  //        	dynamic_cast< EdgeType* >(odomteryEdges[i])->setRobustKernel( new g2o::RobustKernelCauchy() );
		// 	activeEdges.insert( odomteryEdges[i]);
		// }


		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			// if(*it == (std::pair<int,int> (4190, 4085)))
			// 	continue;

         	// 核函数
         	// g2o::RobustKernelHuber();
         	// RobustKernelCauchy
         	// dynamic_cast< EdgeType* >(loopclosureEdges[*it])->setRobustKernel( new g2o::RobustKernelHuber() );	
         	dynamic_cast< EdgeType* >(loopclosureEdges[*it])->setRobustKernel( new g2o::RobustKernelCauchy() );

        	// loopclosureEdges[*it]->setRobustKernel( new g2o::RobustKernelHuber() );			
			// loopclosureEdges[*it].
			activeEdges.insert( loopclosureEdges[*it]);
		}
		cout<<"nIterations is "<<nIterations<<endl;

		IntPairDoublePairMap loopClosureLinkError;
		restore();
		//
		optimizer->setVerbose(false);
	    // optimizer->setVerbose(true);

        endRegion = 0;
		optimizer->vertex(endRegion)->setFixed(true);

		// optimizer->vertices().clear();

		optimizer->initializeOptimization(activeEdges);

		optimizer->optimize(nIterations,false);
		optimizer->computeActiveErrors();
		// optimizer->vertex(endRegion)->setFixed(false);

		double sumLoopChieErr = 0, sum_cauchy =0.0, sum_huber = 0.0, middle_chi;
		g2o::Vector3 rho;
		if(activeLoops.size() >= 1)
		{
			for(
				IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
				it!=end;
				it++)
			{
				error_pair.second = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();	
				sumLoopChieErr = sumLoopChieErr + error_pair.second;
				// loopClosureLinkError[*it] = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();	

				cout<<"loop "<<(*it).first<<" "<<(*it).second<<"  has error "<<error_pair.second<<endl;

				// cauchy 
				// huber_chi2(middle_chi, rho);
				cauchy_chi2(error_pair.second, rho);
				error_pair.first = rho[0];
				// loopClosureLinkError[*it] = rho[0];
				loopClosureLinkError[*it] = error_pair;
				sum_cauchy = sum_cauchy + loopClosureLinkError[*it].first;	

	  			cout<<"*   cauchy    has robust error "<<rho[0]<<"  sum_cauchy:"<<sum_cauchy<<" "<<"  sum_chi:"<<sumLoopChieErr<<endl;
			}

			// NOTE : The number of edges involved is the last element
			// loopClosureLinkError[IntPair(-1,0)] = optimizer->activeChi2(); // There can be no links with negative IDS
			loopClosureLinkError[IntPair(-1,0)] = std::pair<double, double>(optimizer->activeRobustChi2(), 0);
			loopClosureLinkError[IntPair(-1,-1)] = std::pair<double, double>(optimizer->activeEdges().size(), 0);

			loopClosureLinkError[IntPair(-2,0)] = std::pair<double, double>(sum_cauchy, 0); // There can be no links with negative IDS
			loopClosureLinkError[IntPair(-2,-1)] = std::pair<double, double>(activeLoops.size(), 0);		
			//
			cout<<"end optimize_active_robustKernel_pair_return"<<endl;
			return loopClosureLinkError;
		}
	}
	IntPairDoubleMap optimize_robustKernel(const IntPairSet& activeLoops, const int nIterations,
		 std::vector<std::array<double, 5>> & odoEdgeError)	
	{

	// EdgePtrVector 					odomteryEdges;
	// IntPairToEdgePtrMap				loopclosureEdges;
	
	// g2o::OptimizableGraph::EdgeSet 	activeEdges;
	// g2o::OptimizableGraph::VertexSet 	activeVertexes;
		cout<<" "<<endl<<"optimize_robustKernel"<<endl;
		cout<<"nIterations is "<<nIterations<<endl;
		// cout<<"nodes range form "<<startRegion<<" to "<<endRegion<<endl;
		std::pair<std::pair<int, int>, double> odoErrorElement;
		int node, endNode;
		std::array<double, 5> element_error;
		double information_Angle;

		activeEdges.clear();
		activeEdges.insert(odomteryEdges.begin(),odomteryEdges.end());

		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			// if(*it == (std::pair<int,int> (4190, 4085)))
			// 	continue;

         	// 核函数
         	// g2o::RobustKernelHuber();
         	dynamic_cast< EdgeType* >(loopclosureEdges[*it])->setRobustKernel( new g2o::RobustKernelCauchy() );	
        	// loopclosureEdges[*it]->setRobustKernel( new g2o::RobustKernelHuber() );			
			// loopclosureEdges[*it].
			activeEdges.insert( loopclosureEdges[*it]);
		}

		IntPairDoubleMap loopClosureLinkError;
		restore();
		//
		optimizer->setVerbose(false);
	    //optimizer->setVerbose(true);


		optimizer->findGauge()->setFixed(true);

  //       int endRegion = 0;
		// optimizer->vertex(endRegion)->setFixed(true);		
// 
		// optimizer->vertices().clear();

		optimizer->initializeOptimization(activeEdges);

		optimizer->optimize(nIterations,false);
		optimizer->computeActiveErrors();
		optimizer->findGauge()->setFixed(false);

		double sumLoopChieErr = 0, sum_cauchy =0.0, sum_huber = 0.0, middle_chi;
		g2o::Vector3 rho;
		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			middle_chi = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();	
			// loopClosureLinkError[*it] = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();	

			cout<<"loop "<<(*it).first<<" "<<(*it).second<<"  has error "<<middle_chi<<endl;

			// cauchy 
			// huber_chi2(middle_chi, rho);
			cauchy_chi2(middle_chi, rho);
			loopClosureLinkError[*it] = rho[0];
			sumLoopChieErr = sumLoopChieErr + loopClosureLinkError[*it];	

  			cout<<"*   cauchy robust error "<<rho[0]<<"    accumulated error "<<sumLoopChieErr<<endl;
		}

		// NOTE : The number of edges involved is the last element
		loopClosureLinkError[IntPair(-1,0)] = optimizer->activeRobustChi2(); // There can be no links with negative IDS
		loopClosureLinkError[IntPair(-1,-1)] = optimizer->activeEdges().size();

		loopClosureLinkError[IntPair(-2,0)] = sumLoopChieErr; // There can be no links with negative IDS
		loopClosureLinkError[IntPair(-2,-1)] = activeLoops.size();		
		//
		cout<<"end optimize_robustKernel"<<endl;
		return loopClosureLinkError;
	}

	IntPairDoubleMap optimize_active(const IntPairSet& activeLoops, const int nIterations,
		int & startRegion, int & endRegion, std::vector<std::array<double, 5>> & odoEdgeError,
		string ba, string aa)
	{

	// EdgePtrVector 					odomteryEdges;
	// IntPairToEdgePtrMap				loopclosureEdges;
	
	// g2o::OptimizableGraph::EdgeSet 	activeEdges;
	// g2o::OptimizableGraph::VertexSet 	activeVertexes;
		cout<<" optimize_active"<<endl;
		cout<<"nodes range form "<<startRegion<<" to "<<endRegion<<endl;
		std::pair<std::pair<int, int>, double> odoErrorElement;
		int node, endNode;
		std::array<double, 5> element_error;
		double information_Angle;

		odoEdgeError.clear();

		activeEdges.clear();
		// activeEdges.insert(odomteryEdges.begin(),odomteryEdges.end());

		activeEdges.insert(odomteryEdges.begin() + startRegion,odomteryEdges.begin() + endRegion);


		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			// if(*it == (std::pair<int,int> (4190, 4085)))
			// 	continue;
			activeEdges.insert( loopclosureEdges[*it]);
		}

		IntPairDoubleMap loopClosureLinkError;
		restore();
		//
		optimizer->setVerbose(false);
	    //optimizer->setVerbose(true);

		optimizer->vertex(endRegion)->setFixed(true);

		// optimizer->vertices().clear();

		optimizer->initializeOptimization(activeEdges);

		optimizer->optimize(nIterations,false);
		optimizer->computeActiveErrors();
		optimizer->vertex(endRegion)->setFixed(false);

		double sumLoopChieErr = 0;
		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			//dynamic_cast< EdgeType* >(loopclosureEdges[*it])->computeError();
			loopClosureLinkError[*it] = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();

			sumLoopChieErr = sumLoopChieErr + loopClosureLinkError[*it];
		}

		// NOTE : The number of edges involved is the last element
		loopClosureLinkError[IntPair(-1,0)] = optimizer->activeChi2(); // There can be no links with negative IDS
		loopClosureLinkError[IntPair(-1,-1)] = optimizer->activeEdges().size();

		loopClosureLinkError[IntPair(-2,0)] = sumLoopChieErr; // There can be no links with negative IDS
		loopClosureLinkError[IntPair(-2,-1)] = activeLoops.size();		
		//
		cout<<"end optimize_active"<<endl;
		return loopClosureLinkError;
	}

	bool write(const char* filename, const IntPairSet& correctLoops)
	{
		restore();
		std::ofstream out(filename);

		activeEdges.clear();
		activeEdges.insert(odomteryEdges.begin(),odomteryEdges.end());

		for ( IntPairSet::const_iterator it = correctLoops.begin(), end = correctLoops.end();
				it!=end;
				it++)
		{
			activeEdges.insert(loopclosureEdges[*it]);
		}
		optimizer->saveSubset(out,activeEdges);
		out.close();

		return true;
	}
	// optimizer.save("../data/sphere_after.g2o");
	bool saveOptimizedResult(const char* filename, IntPairSet & activeLoops)
	{
		restore();

		std::ofstream out(filename);

		activeEdges.clear();
		activeEdges.insert(odomteryEdges.begin(),odomteryEdges.end());

		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			// if(*it == (std::pair<int,int> (4190, 4085)))
			// 	continue;
			activeEdges.insert( loopclosureEdges[*it]);
		}
		
		//
		optimizer->setVerbose(false);
	    //optimizer->setVerbose(true);
		optimizer->findGauge()->setFixed(true);
		///////////////////
		//std::cerr<<"optimizer->findGauge() output: "<<(*(optimizer->findGauge()))<<std::endl;;
		//optimizer->vertex(0)->setFixed(true);
		printf("This is in %s on line %d\n",  __FILE__, __LINE__);
		optimizer->initializeOptimization(activeEdges);
		printf("This is in %s on line %d\n",  __FILE__, __LINE__);

		optimizer->optimize(10,false);
		printf("This is in %s on line %d\n",  __FILE__, __LINE__);

		optimizer->computeActiveErrors();
		printf("This is in %s on line %d\n",  __FILE__, __LINE__);

		// store();
		printf("This is in %s on line %d\n",  __FILE__, __LINE__);
		// optimizer->save(filename);
		optimizer->saveSubset(out,activeEdges);
		printf("This is in %s on line %d\n",  __FILE__, __LINE__);
		cout<<"active chi2 in saveOptimized Result : "<<optimizer->activeChi2()<<endl; 
		printf("This is in %s on line %d\n",  __FILE__, __LINE__);
		out.close();
		return true;
	}

	int edgeDimension()
	{
		return EdgeType::Dimension;
	}

	bool isInitialized()
	{
		return initialized;
	}

	bool removeEdges(const IntPairSet& falseLinks)
	{
		for (IntPairSet::const_iterator it = falseLinks.begin(), end = falseLinks.end(); it!=end; ++it)
		{
			if(optimizer->removeEdge(loopclosureEdges[*it]))
			{
				loopclosureEdges.erase(*it); // Clear up from our records as well
			}
		}
		return true;
	}
};


#endif /* BACKEND_G2O_HPP_ */
