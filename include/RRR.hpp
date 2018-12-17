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


#ifndef RRR_HPP_
#define RRR_HPP_


#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <fstream> 
#include <sys/time.h>

#include "g2o_Interface.hpp"
#include "g2o/types/slam2d/se2.h"
#include "g2o/types/slam2d/edge_se2.h"

#include "cluster.hpp"
#include "utils.hpp"

#include "g2o_interface_se2_zihao.hpp"


template <class BackEndInterface> //!< @param BackEndInterace one of the back-end
class RRR
{
	Clusterizer clusterizer;
	

	double clusteringThreshold;
	int nIterations;
	int ID_IGNORE;
	string cluster_filename;
	string OptimizedResultWithoutRecheck;//file name of the optimized result with     out     doubt loop check
	string OptimizedResultWithRecheck;	//file name of the optimized result with doubt loop check
	string OptimizedResultWithReRecheck, OptimizedResultBeforInter;
	string final_suvived_loops_fileName, OptimizedResultBeforactive, OptimizedResultAfteractive;
	string OptimizedResultIntra ;
	string OptimizedResultIntraInter ;
	string OptimizedResultIntraConflict ;
	string OptimizedResultIntraConflictInter ;	
   // const char *OPre=OptiRefile.data();
	std::vector<std::pair<std::pair<int, int>, double> >  odoEdgeRelateLC_Error;
	std::vector<std::array<double, 5>>  odoEdgeError;
	int edgeDimension;

	IntSet _goodSet;
	// double fiveNineBelief = 28*2;
	const double twoNineBelief = 11.345, fiveNineBelief = 28, nineNineBelief = 44.8413, twelveBelief = 58.9413, 
		upperLimitChi = 77.4, tenUpperLimit = 774;


public:
	double deltaT;
	int allfalseAccepted = 0;
	int allTest = 0, lostGood = 0;
	int allTrue = 0;
	std::vector<std::pair<int, double>> SmallSumChiCluster;
	std::pair<int, double> eleForSamllSCCluster;
	find_element find_ele;
	BackEndInterface* gWrapper;
	std::map<int, double> overallChi2SET4singleElementCluster;
	int small_numIterations = 5;
	int overOptimizeIteration =12;

	RRR(
		double ClusteringThreshold = 1,
		int numIterations = 5,
		string Cluster_filename = "unmane"
			):
	clusteringThreshold(ClusteringThreshold),
				nIterations(numIterations),
				cluster_filename(Cluster_filename)
	{
		ID_IGNORE 			= -2;
		gWrapper 			= new BackEndInterface;
		edgeDimension 		= gWrapper->edgeDimension();

		string fullName = cluster_filename;
		int posfile = fullName.find_last_of('/');
		string path = fullName.substr(0, posfile+1); 
		string  fileName=fullName.substr(posfile+1);
	    fileName = fileName.substr(0, fileName.length()-4);

	    // resultfile = path+"result_"+fileName+'_'+".g2o";

		// clusterizer.nameofclusterfile   =    "zihao_cluster"+fileName+'_'+to_string(int(clusteringThreshold*100))+".g2o";
		clusterizer.nameofclusterfile   =    "zihao_cluster"+fileName+"_.g2o";

		OptimizedResultWithoutRecheck = path+"result_zihao_"+fileName+"_.g2o";

		OptimizedResultBeforInter =  path+"result_zihao_"+fileName+".g2o";

		OptimizedResultIntra =  path+"result_zihao_"+fileName+"_1.g2o";
		OptimizedResultIntraInter = path+"result_zihao_"+fileName+"_2.g2o";

	}



	bool read(const char* filename)
	{
		// cout<<" read re why exit quietly? "<<endl;
		IntPairSet loops;
		bool status = gWrapper->read(filename);		// Read the g2o file
		if(!status) return status;
		gWrapper->getLoopClosures(loops);
		
		std::cerr<<"Number of loop closures : "<<loops.size()<<std::endl;

		clusterizer.clusterize_zihao(loops,filename);

				// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				// exit(0);
		return true;
	}

	bool setOptimizer(void* optimizerPtr, const char* filename)
	{
		cout<<"xx re why exit quietly? "<<endl;
		IntPairSet loops;
		gWrapper->setOptimizer(optimizerPtr);
		gWrapper->getLoopClosures(loops);
		
		std::cerr<<"Number of loop closures : "<<loops.size()<<std::endl;

		clusterizer.clusterize_zihao(loops,filename);

		return true;
	}


	
	bool intraClusterConsistent_goodSET_handleOneMenCluster(int clusterID, std::vector<double> & chiStatisVector, std::vector<std::pair<int, int> > & doubtLoopNumberVector, 
		int & originalLoopNum, std::vector<std::pair<double, std::pair<int, int> > > & chiStstis4eachLoop)
	{
		// cout<<" "<<endl;
		// cout<<"second is running"<<endl;
		// std::cerr<<" chi2 %95 2: "<<utils::chi2(2)<<std::endl;
		IntPairDoubleMap chi2LinkErrors, chi2LinkErrors2;
		IntPairSet& currentCluster = clusterizer.getClusterByID(clusterID);

		int originalTotalNum=0;
		originalTotalNum=currentCluster.size();
		cout<<"cluster "<<clusterID<<" has "<<originalTotalNum<<" loops"<<endl;
		int anotherNum = clusterizer._clustersFound[clusterID].positionserial.size();
		if(originalTotalNum != anotherNum)
		{
			cout<<"cluster "<<clusterID<<" has "<<anotherNum<<" loops"<<endl;
			cout<<"the member size in this cluster got form two ways are not equal "<<endl;
		}
		originalLoopNum = originalTotalNum;
		int realOrigianlNum = originalLoopNum;
		std::vector<double>  eachChi2Error, dis_clu;
		std::vector<std::pair<std::pair<int,int> , double> > chi2AndLoopPair;
		std::pair<std::pair<int,int> , double> eleffssxzh;
		std::vector<int>  mulCluster;
		//the third parameter should be set to negative things, it is used to skip some certain loop when debug
		int startRange, endRange;
		std::set<int> clusterNodeSet;

		IntPairSet::iterator findRange = currentCluster.begin(), findRangeEnd = currentCluster.end();
		for(; findRange != findRangeEnd; findRange++)
		{
			clusterNodeSet.insert((*findRange).first);
			clusterNodeSet.insert((*findRange).second);
		}

		startRange = *clusterNodeSet.begin();
		endRange   = *clusterNodeSet.rbegin();

        cout<<"node range from "<<startRange<<" to "<<endRange<<endl;
        cout<<"debug 3"<<endl;
		//pick out the subgraph corresponding to the node range
		string ba = "sdf", aa = "aa";
		// chi2LinkErrors =  gWrapper->optimize_active(currentCluster, small_numIterations, startRange, endRange, odoEdgeError, ba, aa);
		chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);

        cout<<"debug 4"<<endl;
		chiStstis4eachLoop.clear();
		doubtLoopNumberVector.clear();

		float activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
		int   activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];

		cout<<"active edge count: "<<activeEdgeCount<<endl;
		// if( currentCluster.find(std::pair<int,int>(1875, 0)) != currentCluster.end())

		double belieftoDelete = 0.95;
		double meanError = 0;
		int countLoop = 0;
		std::vector<std::pair<int,int> > badloop ;
		IntPairSet good_set, good_set2, good_set3;//bad_set
		std::vector<std::pair<int,int> > ID_IGNORE_vector;

		//check if there is someone loop whose chi2 error is bigger than the signgle edge limit error
		double biggestErr = 0, smallcHI = 200;
		std::pair<int, int> biggestErrorLoop(-1, -1);

		cout<<" "<<endl;
		IntPairDoubleMap::iterator eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		double goodLoop = 0;
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;
			eleffssxzh.first  = eIt->first;
			eleffssxzh.second = eIt->second;
			chi2AndLoopPair.push_back(eleffssxzh);
			dis_clu.push_back(eIt->second);
			if(eIt->second > biggestErr)
			{
				biggestErr = eIt->second;
				biggestErrorLoop = eIt->first;
			}
			if(eIt->second < smallcHI)
			{
				smallcHI = eIt->second;
			}	
			if(eIt->second < utils::chi2(edgeDimension))		
			{
				goodLoop=goodLoop+1;
			}

			cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
				" has error "<< eIt->second<<endl;
		}
		// if(clusterID == 85 and clusterizer.loopToClusterIDMap[std::pair<int,int> (9016, 3876)] != 85)
		// 	exit(0);

		// if(activeChi2Graph > utils::chi2_continuous(edgeDimension* floor(activeEdgeCount/10), 0.95))
		if(1 == originalTotalNum )
		{
			// if(activeChi2Graph > utils::chi2_continuous(edgeDimension* ceil(activeEdgeCount/10), 0.95) or 
			// 	activeChi2Graph > 50)
			// {
			// 	cout<<"this cluster has one loop, and the activeChi2Graph is bigger than the overall limit, so delete this cluster"<<endl;
			// 	cout<<"activeChi2Graph: "<<activeChi2Graph <<" activeEdgeCount: "<<activeEdgeCount<<" threshold: "<<
			// 		utils::chi2_continuous(edgeDimension* floor(activeEdgeCount/10), 0.1) <<endl;
			// 	cout<<"edgeCount: "<<activeEdgeCount<<" start range: "<<startRange <<" end range: "<<endRange<<" "<<endRange- startRange<<endl;
			// 	cout<<"gWrapper->odoSize: "<<gWrapper->odoSize<<endl;		
			// 	std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;

			// 	return false;
			// }

			if(activeChi2Graph > utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95) or 
				chi2LinkErrors[IntPair(-2,0)] > utils::chi2(edgeDimension))
			// if(activeChi2Graph > 500 or 
			// 	chi2LinkErrors[IntPair(-2,0)] > utils::chi2(edgeDimension))			
			{
				cout<<"this cluster has one loop, and the activeChi2Graph is bigger than the overall limit, so delete this cluster"<<endl;
				cout<<"activeChi2Graph: "<<activeChi2Graph <<" activeEdgeCount: "<<activeEdgeCount<<" threshold: "<<
					utils::chi2_continuous(edgeDimension* floor(activeEdgeCount/10), 0.1) <<endl;
				cout<<"edgeCount: "<<activeEdgeCount<<" start range: "<<startRange <<" end range: "<<endRange<<" "<<endRange- startRange<<endl;
				cout<<"gWrapper->odoSize: "<<gWrapper->odoSize<<endl;		
				std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;

				return false;
			}			

		} 

		if(activeChi2Graph > utils::chi2(edgeDimension* activeEdgeCount) and goodLoop == 0)
		{
			cout<<"activeChi2Graph: "<<activeChi2Graph<<" is bigger than threshold "<<endl;
			cout<<"empty cluster "<<clusterID<<" has "<<originalTotalNum<<" loop "<<endl;
			return 0;
		}

		//delete the loop withgiggest error that bigger than upperlimit then redo opeimize		
		if(biggestErr > upperLimitChi and smallcHI< upperLimitChi)
		{
			// clusterizer.setClusterID(biggestErrorLoop,ID_IGNORE);
				clusterizer.getClusterID_multiExist(biggestErrorLoop, mulCluster);
				if(mulCluster.size() != 0)
				{
					clusterizer.clusterIDtoLoopsMap[clusterID].erase(biggestErrorLoop);
					clusterizer.deleteClusterID_multiExist(biggestErrorLoop,  clusterID);
				}
				else
					clusterizer.setClusterID(biggestErrorLoop,ID_IGNORE); 

			cout<<"delete a large error edge that perhaps damage other loop"<<endl;
			// cin.get();
			// currentCluster.clear();
			currentCluster = clusterizer.getClusterByID(clusterID);
			originalTotalNum = currentCluster.size();
			//the third parameter should be set to negative things, it is used to skip some certain loop when debug
			chi2LinkErrors.clear();
			chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
			activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
		}



		eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		cout<<" "<<endl;
		cout<<"*********************** first time to find the bad loops ************************"<<endl;
		cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<
			utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;
			
			chiStstis4eachLoop.push_back(std::pair<double, std::pair<int, int> > (eIt->second, eIt->first));
			meanError = meanError + eIt->second;
			countLoop = countLoop+1;

			if(eIt->second > utils::chi2(edgeDimension, belieftoDelete))
			// if(eIt->second > 100)
			{
				// bad_set.insert(eIt->first);
				badloop.push_back(eIt->first);//save the bad loop nodes
				cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;

			}
			else
			{
				cout<<"good loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;
				good_set.insert(eIt->first);
			}
		}
		cout<<"**********************************************************************************"<<endl;
		cout<<" "<<endl;
		// if all loops are good loop, return true
		if(good_set.size() == originalTotalNum )
		{
			if( activeChi2Graph < utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95))
			{
				cout<<"the good set number equal to the total number, so suvived"<<endl;
				return true;
			}
			else
			{
				cout<<"the good set number equal to the total number, but whole graph error bigger than two loop edge "<<endl;
				return false;
			}	

		}

		cout<<" "<<endl;
		cout<<"*********************** optimise on good loop set 1 ************************"<<endl;
		if(good_set.size() > 0 and badloop.size() > 0 )
		{
			chi2LinkErrors.clear();
			chi2LinkErrors = gWrapper->optimize(good_set, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();

			activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
			activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];

			cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<
				utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;
			for(; eIt!=eEnd; eIt++)
			{
				if(eIt->first.first < 0)
					continue;
				
				chiStstis4eachLoop.push_back(std::pair<double, std::pair<int, int> > (eIt->second, eIt->first));
				meanError = meanError + eIt->second;
				countLoop = countLoop+1;

				if(eIt->second > utils::chi2(edgeDimension,belieftoDelete))
				// if(eIt->second > 100)
				{
					// bad_set.insert(eIt->first);
					badloop.push_back(eIt->first);//save the bad loop nodes
					cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;

				}
				else
				{
					cout<<"good loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;
					good_set2.insert(eIt->first);
				}
			}
			cout<<"**********************************************************************************"<<endl;
			cout<<" "<<endl;

			eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
			int goodLoopNum1 = good_set.size(),goodLoopNum2 = good_set2.size();

			while(goodLoopNum1 > goodLoopNum2)
			{
				goodLoopNum1 = goodLoopNum2;
				good_set3.clear();
				cout<<" "<<endl;
				cout<<"*********************** optimise on good loop set 2 ************************"<<endl;

					chi2LinkErrors.clear();
					chi2LinkErrors = gWrapper->optimize(good_set2, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
					activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
					activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];

				eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();

				cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<
					utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;
				for(; eIt!=eEnd; eIt++)
				{
					if(eIt->first.first < 0)
						continue;
					
					chiStstis4eachLoop.push_back(std::pair<double, std::pair<int, int> > (eIt->second, eIt->first));
					meanError = meanError + eIt->second;
					countLoop = countLoop+1;

					if(eIt->second > utils::chi2(edgeDimension,belieftoDelete))
					// if(eIt->second > 100)
					{
						// bad_set.insert(eIt->first);
						badloop.push_back(eIt->first);//save the bad loop nodes
						cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;

					}
					else
					{
						cout<<"good loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;
						good_set3.insert(eIt->first);
					}
				}
				cout<<"**********************************************************************************"<<endl;
				cout<<" "<<endl;
				goodLoopNum2 = good_set3.size();
				if(goodLoopNum2 == 0)
					break;
				if(goodLoopNum2 == goodLoopNum1)
					break;
				good_set2.clear();
				good_set2.insert(good_set3.begin(), good_set3.end());
			}
		}

		//although badloop equal to original elements number ,we shoudl not delete it without futher check
		if(badloop.size() == originalTotalNum)
		{
			cout<<"the bad loop number equals to the total loop number so delete the cluster"<<endl;
			std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
			return false;
		}


		if(badloop.size() != 0)
		{
			for(int i = 0; i < badloop.size(); i++)
			{
				
				clusterizer.getClusterID_multiExist(badloop[i], mulCluster);
				if(mulCluster.size() != 0)
				{
					clusterizer.clusterIDtoLoopsMap[clusterID].erase(badloop[i]);
					clusterizer.deleteClusterID_multiExist(badloop[i],  clusterID);
				}
				else
					clusterizer.setClusterID(badloop[i],ID_IGNORE); 
			}	
		}

		if(clusterizer.getClusterByID(clusterID).size() >= 1 )
		{
			std::cerr<<" Cluster "<<clusterID <<"survived with "<<clusterizer.getClusterByID(clusterID).size()<<" links "<<"from "<<realOrigianlNum<<std::endl;
		
			if(originalTotalNum == 1)
			{
				overallChi2SET4singleElementCluster[clusterID] = activeChi2Graph;	
			}
			if( activeChi2Graph < utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95))
			{
				return true;
			}
			else
			{
				cout<<"the good set number equal to the total number, but whole graph error bigger than two loop edge"<<endl;
				return false;
			}	
		}
		else
		{
			std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
			return false;
		}
		// // exit(0);
		return false;
	}
	void generateNewEdge_parallel(const std::pair<int, int> & loop_node_pair, 
						std::vector< g2o::HyperGraph::Edge* > & newEdgeVector,
						std::vector<std::pair<int,int > > & artificialLoops, int futileDis, int maxVertexID,
						bool adjoint)
	{
			int start = loop_node_pair.first, end = loop_node_pair.second;
			bool fultileBit1;
			g2o::SE2  TransOdoSeg, TransLoop, newTrans, new_trans_start, containerINverse, newTrans_final;
			std::pair<g2o::SE2, Matrix3d > result_synthetic_odo;
			Matrix3d  covOdoSeg, covOdoSeg_start, covLoop, newCov, j1, j2, newCov_final;
			std::array<std::array<double, 5>, 5>  length;
			int nodeNum = clusterizer.VertexInf.size();

			// construct the test loop SE2 MATRIX3D
			TransLoop = clusterizer.LP_Trans_Covar_Map[loop_node_pair].first;
			containerINverse = TransLoop.inverse();
			covLoop = clusterizer.LP_Trans_Covar_Map[loop_node_pair].second;
			cout<<"cov of "<<start <<" "<<end <<" is "<<endl<<covLoop<<" "<<endl;
			// cin.get();
         	cout<<"loop "<<start<<" "<<end<<endl;
         	if(abs(end - start) > 20)
         	{
	         	int nodeStep = (end - start)/19, nextNode, nextNode_1; 
         		// nodeStep = 1*nodeStep/abs(nodeStep);
	         	if(abs(nodeStep) > futileDis)
	         	{
	         		int mid = (nodeStep)/(abs(nodeStep))*futileDis;
	         		nodeStep = mid;
	         	}
	         	cout<<"node step is "<<nodeStep<<endl;
	         	for(int i = 0; i < 4; i++)
	         	{
	         		// create a new edge
	         		g2o::EdgeSE2* newEdge =  new g2o::EdgeSE2(); 
	         		// calculate the next node
	         		nextNode_1 = end - (i+1)*nodeStep;
	         		
	         		// cout<<"start "<<start<<" nextNode "<<nextNode<<endl;
	         		//initialize the vertex  g2o::VertexSE2, g2o::EdgeSE2

	         		// calculate the measurement and information matrix
	         			// cal the odo part
	         		if(end < nextNode_1)
					{
						cout<<"end odo segment "<<end<<" to "<<nextNode_1<<endl;

						if(adjoint)
							clusterizer.synthesize_odo_edges_adjoint(end, nextNode_1, clusterizer.OdoInf, result_synthetic_odo.first, 
								result_synthetic_odo.second, length[0], fultileBit1);
						else
							clusterizer.accelerated_synthetic_odo(end, nextNode_1, clusterizer.OdoInf, 
								result_synthetic_odo,  length[0], fultileBit1);

						TransOdoSeg = result_synthetic_odo.first;
						covOdoSeg   = result_synthetic_odo.second;				
					}
	         		else
					{
						cout<<"end odo segment re "<<nextNode_1<<" to "<<end<<endl;
						if(adjoint)
							clusterizer.synthesize_odo_edges_adjoint(nextNode_1, end, clusterizer.OdoInf, result_synthetic_odo.first, 
								result_synthetic_odo.second, length[0], fultileBit1);
						else
							clusterizer.accelerated_synthetic_odo(nextNode_1, end, clusterizer.OdoInf, 
								result_synthetic_odo,  length[0], fultileBit1);
						TransOdoSeg = result_synthetic_odo.first.inverse();
						covOdoSeg   = result_synthetic_odo.second;				
					}
					// cout<<"TransOdoSeg "<<TransOdoSeg[0]<<" "<<TransOdoSeg[1]<<" "<<TransOdoSeg[2]<<endl;
					// cout<<"TransLoop   "<<TransLoop[0]<<" "<<TransLoop[1]<<" "<<TransLoop[2]<<endl;
					newTrans = TransLoop * TransOdoSeg;
					// cout<<"synthesized trnas: "<<newTrans[0]<<" "<<newTrans[1]<<" "<<newTrans[2]<<endl;
					// cout<<

					clusterizer.Jacobian_4_edge_propagate(TransLoop, TransOdoSeg, j1, j2);
					clusterizer.covariance_propagate(covLoop, covOdoSeg, j1, j2, newCov);


					//calculate the start segment
					nextNode = start + (i+1)*nodeStep;
	         		if(start < nextNode)
					{
						cout<<"start odo segment re "<<start<<" to "<<nextNode<<endl;
						if(adjoint)
							clusterizer.synthesize_odo_edges_adjoint(start, nextNode, clusterizer.OdoInf, result_synthetic_odo.first, 
								result_synthetic_odo.second, length[0], fultileBit1);
						else
							clusterizer.accelerated_synthetic_odo(start, nextNode, clusterizer.OdoInf, 
								result_synthetic_odo,  length[0], fultileBit1);
						new_trans_start = result_synthetic_odo.first.inverse();
						covOdoSeg_start   = result_synthetic_odo.second;				
					}
	         		else
					{
						cout<<"start odo segment "<<nextNode<<" to "<<start<<endl;
						if(adjoint)
							clusterizer.synthesize_odo_edges_adjoint(nextNode, start, clusterizer.OdoInf, result_synthetic_odo.first, 
								result_synthetic_odo.second, length[0], fultileBit1);
						else						
							clusterizer.accelerated_synthetic_odo(nextNode, start, clusterizer.OdoInf, 
								result_synthetic_odo,  length[0], fultileBit1);
						new_trans_start = result_synthetic_odo.first;
						covOdoSeg_start   = result_synthetic_odo.second;				
					}
					// cout<<"TransOdoSeg "<<TransOdoSeg[0]<<" "<<TransOdoSeg[1]<<" "<<TransOdoSeg[2]<<endl;
					// cout<<"TransLoop   "<<TransLoop[0]<<" "<<TransLoop[1]<<" "<<TransLoop[2]<<endl;

					newTrans_final = new_trans_start * newTrans;
					cout<<"synthesized trnas: "<<newTrans_final[0]<<" "<<newTrans_final[1]<<" "<<newTrans_final[2]<<endl;
					// cout<<

					clusterizer.Jacobian_4_edge_propagate(new_trans_start, newTrans, j1, j2);
					clusterizer.covariance_propagate(covOdoSeg_start, newCov, j1, j2, newCov_final);
					cout<<"synthesized edge "<<nextNode <<" "<< nextNode_1<<" cov is "<<endl<<newCov_final<<" "<<endl;
					cout<<" "<<endl;

	         		// assign measurement and information to it
	         		newEdge->setVertex( 0, dynamic_cast<g2o::VertexSE2*> (gWrapper->optimizer->vertex(nextNode)) );
	         		newEdge->setVertex( 1, dynamic_cast<g2o::VertexSE2*> (gWrapper->optimizer->vertex(nextNode_1)) );	         		
	         		newEdge->setMeasurement( newTrans_final );
	         		newEdge->setInformation( newCov_final.inverse() );
	         		// store the new edge 
	         		newEdgeVector.push_back(newEdge);
	         		artificialLoops.push_back(std::pair<int, int>(nextNode, nextNode_1));
	         	}

         	}
         	else
         	{
         		cout<<"the node dis between start and end node <= 4, so do not make new edge for it"<<endl;
         		return;
         	}
	}
	void generateNewEdge(const std::pair<int, int> & loop_node_pair, 
						std::vector< g2o::HyperGraph::Edge* > & newEdgeVector,
						std::vector<std::pair<int,int > > & artificialLoops, int futileDis, int maxVertexID)
	{
			int start = loop_node_pair.first, end = loop_node_pair.second;
			bool fultileBit1;
			g2o::SE2  TransOdoSeg, TransLoop, newTrans, containerINverse;
			std::pair<g2o::SE2, Matrix3d > result_synthetic_odo;
			Matrix3d  covOdoSeg, covLoop, newCov, j1, j2;
			std::array<std::array<double, 5>, 5>  length;

			// construct the test loop SE2 MATRIX3D
			TransLoop = clusterizer.LP_Trans_Covar_Map[loop_node_pair].first;
			containerINverse = TransLoop.inverse();
			covLoop = clusterizer.LP_Trans_Covar_Map[loop_node_pair].second;
			// cout<<"covLoop is "<<covLoop<<endl;
			// cin.get();
         	cout<<"loop "<<start<<" "<<end<<endl;
         	if(abs(end - start) > 20)
         	{
	         	int nodeStep = (end - start)/20, nextNode;
         		nodeStep = nodeStep/abs(nodeStep);
	         	if(abs(nodeStep) > futileDis)
	         	{
	         		int mid = (nodeStep)/(abs(nodeStep))*futileDis;
	         		nodeStep = mid;
	         	}
	         	// cout<<"node step is "<<nodeStep<<endl;
	         	for(int i = 0; i < 5; i++)
	         	{
	         		// create a new edge
	         		g2o::EdgeSE2* newEdge =  new g2o::EdgeSE2(); 
	         		// calculate the next node
	         		nextNode = end - (i+1)*nodeStep;
	         		artificialLoops.push_back(std::pair<int, int>(start, nextNode));
	         		// cout<<"start "<<start<<" nextNode "<<nextNode<<endl;
	         		//initialize the vertex  g2o::VertexSE2, g2o::EdgeSE2
	         		newEdge->setVertex( 0, dynamic_cast<g2o::VertexSE2*> (gWrapper->optimizer->vertex(start)) );
	         		newEdge->setVertex( 1, dynamic_cast<g2o::VertexSE2*> (gWrapper->optimizer->vertex(nextNode)) );
	         		// calculate the measurement and information matrix
	         			// cal the odo part
	         		if(end < nextNode)
					{
						// cout<<"odo segment "<<end<<" to "<<nextNode<<endl;
						clusterizer.accelerated_synthetic_odo(end, nextNode, clusterizer.OdoInf, 
							result_synthetic_odo,  length[0], 
							fultileBit1);
						TransOdoSeg = result_synthetic_odo.first;
						covOdoSeg   = result_synthetic_odo.second;				
					}
	         		else
					{
						// cout<<"odo segment "<<nextNode<<" to "<<end<<endl;
						clusterizer.accelerated_synthetic_odo(nextNode, end, clusterizer.OdoInf, 
							result_synthetic_odo,  length[0], 
							fultileBit1);
						TransOdoSeg = result_synthetic_odo.first.inverse();
						covOdoSeg   = result_synthetic_odo.second;				
					}
					// cout<<"TransOdoSeg "<<TransOdoSeg[0]<<" "<<TransOdoSeg[1]<<" "<<TransOdoSeg[2]<<endl;
					// cout<<"TransLoop   "<<TransLoop[0]<<" "<<TransLoop[1]<<" "<<TransLoop[2]<<endl;
					newTrans = TransLoop * TransOdoSeg;
					// cout<<"synthesized trnas: "<<newTrans[0]<<" "<<newTrans[1]<<" "<<newTrans[2]<<endl;
					// cout<<

					clusterizer.Jacobian_4_edge_propagate(TransLoop, TransOdoSeg, j1, j2);
					clusterizer.covariance_propagate(covLoop, covOdoSeg, j1, j2, newCov);

	         		// assign measurement and information to it
	         		newEdge->setMeasurement( newTrans );
	         		newEdge->setInformation( newCov.inverse() );
	         		// store the new edge 
	         		newEdgeVector.push_back(newEdge);
	         	}

	  /*       	for(int i = 0; i < 4; i++)
	         	{

	         		// create a new edge
	         		g2o::EdgeSE2* newEdge =  new g2o::EdgeSE2(); 
	         		// calculate the next node
	         		nextNode = end + (i+1)*nodeStep;
	         		if(nextNode > maxVertexID)
	         			continue;
	         		artificialLoops.push_back(std::pair<int, int>(start, nextNode));
	         		cout<<"nextNode "<<start<<" end "<<nextNode<<endl;
	         		//initialize the vertex  g2o::VertexSE2, g2o::EdgeSE2
	         		newEdge->setVertex( 0, dynamic_cast<g2o::VertexSE2*> (gWrapper->optimizer->vertex(start )) );
	         		newEdge->setVertex( 1, dynamic_cast<g2o::VertexSE2*> (gWrapper->optimizer->vertex(nextNode)) );
	         		// calculate the measurement and information matrix
	         			// cal the odo part
	         		if(end < nextNode)
					{
						cout<<"odo segment "<<end<<" to "<<nextNode<<endl;
						clusterizer.accelerated_synthetic_odo(end, nextNode, clusterizer.OdoInf, 
							result_synthetic_odo,  length[0], 
							fultileBit1);
						TransOdoSeg = result_synthetic_odo.first;
						covOdoSeg   = result_synthetic_odo.second;				
					}
	         		else
					{
						cout<<"odo segment "<<nextNode<<" to "<<start<<endl;
						clusterizer.accelerated_synthetic_odo(nextNode, end, clusterizer.OdoInf, 
							result_synthetic_odo,  length[0], 
							fultileBit1);
						TransOdoSeg = result_synthetic_odo.first.inverse();
						covOdoSeg   = result_synthetic_odo.second;				
					}
					newTrans = TransLoop *TransOdoSeg;

					clusterizer.Jacobian_4_edge_propagate(TransLoop, TransOdoSeg, j1, j2);
					clusterizer.covariance_propagate(covLoop, covOdoSeg, j1, j2, newCov);

	         		// assign measurement and information to it
	         		newEdge->setMeasurement( newTrans );
	         		newEdge->setInformation( newCov.inverse() );
	         		// store the new edge 
	         		newEdgeVector.push_back(newEdge);
	         	}	         	*/
    //      		// create a new edge
    //      		g2o::EdgeSE2* newEdge =  new g2o::EdgeSE2(); 
    //      		// calculate the next node
    //      		nextNode = end + nodeStep;
    //      		if(nextNode > maxVertexID or nextNode < 0)
    //      			return;
    //      		artificialLoops.push_back(std::pair<int, int>(start, nextNode));
    //      		cout<<"nextNode  is "<<nextNode<<endl;
    //      		//initialize the vertex  g2o::VertexSE2, g2o::EdgeSE2
    //      		newEdge->setVertex( 0, dynamic_cast<g2o::VertexSE2*> (gWrapper->optimizer->vertex(start)) );
    //      		newEdge->setVertex( 1, dynamic_cast<g2o::VertexSE2*> (gWrapper->optimizer->vertex(nextNode)) );
    //      		// calculate the measurement and information matrix
    //      			// cal the odo part
    //      		if(end < nextNode)
				// {
				// 	clusterizer.accelerated_synthetic_odo(end, nextNode, clusterizer.OdoInf, 
				// 		result_synthetic_odo,  length[0], 
				// 		fultileBit1);
				// 	TransOdoSeg = result_synthetic_odo.first;
				// 	covOdoSeg   = result_synthetic_odo.second;				
				// }
    //      		else
				// {
				// 	clusterizer.accelerated_synthetic_odo(nextNode, end, clusterizer.OdoInf, 
				// 		result_synthetic_odo,  length[0], 
				// 		fultileBit1);
				// 	TransOdoSeg = result_synthetic_odo.first;
				// 	covOdoSeg   = result_synthetic_odo.second;				
				// }
				// newTrans = TransLoop * TransOdoSeg;

				// clusterizer.Jacobian_4_edge_propagate(TransLoop, TransOdoSeg, j1, j2);
				// clusterizer.covariance_propagate(covLoop, covOdoSeg, j1, j2, newCov);

    //      		// assign measurement and information to it
    //      		newEdge->setMeasurement( newTrans );
    //      		newEdge->setInformation( newCov.inverse() );
    //      		// store the new edge 
    //      		newEdgeVector.push_back(newEdge);

         	}
         	else
         	{
         		cout<<"the node dis between start and end node <= 4, so do not make new edge for it"<<endl;
         		return;
         	}
	}
	static bool cmp_reverse(  const std::pair<int, int> p,  const std::pair<int, int> q)
	{
		return p.second < q.second;
	}
	// optimize under robust kernel, and then select the error in non_robust_kernel sense larger than %95 as doubt one, i
	// if there are some loops can prove its genuity by transform distance and no one conflict with it , keep it ,else delete it.
	bool complementary_intra_third(int clusterID, std::vector<double> & chiStatisVector, std::vector<std::pair<int, int> > & doubtLoopNumberVector, 
		int & originalLoopNum, std::vector<std::pair<double, std::pair<int, int> > > & chiStstis4eachLoop,
		std::vector<std::pair<int, int> > & deleted_loops)
	{
		IntPairDoublePairMap chi2LinkErrors;
		IntPairSet& currentCluster = clusterizer.getClusterByID(clusterID);

		std::vector<std::pair<int,int> > goodLoops, badLoops,real_good_loops;
		int originalTotalNum = currentCluster.size();
		cout<<"cluster "<<clusterID<<" has "<<originalTotalNum<<" loops"<<endl;
		originalLoopNum = originalTotalNum;
		int realOrigianlNum = originalLoopNum;

		//the third parameter should be set to negative things, it is used to skip some certain loop when debug
		int startRange, endRange;
		std::set<int> clusterNodeSet;
		IntPairSet::iterator findRange = currentCluster.begin(), findRangeEnd = currentCluster.end();
		for(; findRange != findRangeEnd; findRange++)
		{
			clusterNodeSet.insert((*findRange).first);
			clusterNodeSet.insert((*findRange).second);
		}
		startRange = *clusterNodeSet.begin();
		endRange   = *clusterNodeSet.rbegin();
        cout<<"node range from "<<startRange<<" to "<<endRange<<endl;

        // if originalTotalNum = 1, 
        string ba = "sdf", aa = "aa";

        std::vector< g2o::HyperGraph::Edge* > newEdgeVector;
        if(originalTotalNum == 1)
        {
        	std::vector<std::pair<int, int>> artificialLoops;
        	int maxVertexID = gWrapper->optimizer->vertices().size();
        	int futileDis;
        	if(clusterizer.futilePairVector.size() != 0)
        		futileDis = abs(clusterizer.futilePairVector[0].first - clusterizer.futilePairVector[0].second);
        	else
        		futileDis = maxVertexID;
        	
        	// generateNewEdge(*(currentCluster.begin()), newEdgeVector, artificialLoops, futileDis, maxVertexID);
			chi2LinkErrors =  
				// gWrapper->optimize_oneMember(currentCluster, small_numIterations, newEdgeVector, artificialLoops);
				gWrapper->optimize_active_robustKernel_pair_return(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);
        }
        else
        {
			chi2LinkErrors =  		
				gWrapper->optimize_active_robustKernel_pair_return(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);
        }

		float activeChi2Graph = chi2LinkErrors[IntPair(-1,0)].first;
		int   activeEdgeCount = int(chi2LinkErrors[IntPair(-1,-1)].first);
		cout<<"active edge count: "<<activeEdgeCount<<endl;

					cout<<"activeChi2Graph: "<<chi2LinkErrors[IntPair(-1,0)].first <<" activeEdgeCount: "<<activeEdgeCount<<endl;
					cout<<"clsuter error: "<<chi2LinkErrors[IntPair(-2,0)].first <<" originalLoopNum: "<<originalLoopNum<<endl;


		//check if there is someone loop whose chi2 error is bigger than the signgle edge limit error
		double biggestErr = 0, smallcHI = 200;
		std::pair<int, int> biggestErrorLoop(-1, -1);

		cout<<" "<<endl;
		IntPairDoublePairMap::iterator eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		double goodLoop = 0;
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;

			// if(eIt->second < utils::chi2_continuous(3, 0.5))
			// {
			// 	real_good_loops.push_back(eIt->first);
			// }
			//optimieze with robust kernel, if the loop error under non-kernel smaller than %95, it is regard as a good one
			if(eIt->second.second < utils::chi2_continuous(3, 0.95))	
			{
				goodLoop=goodLoop+1;
				goodLoops.push_back(eIt->first);
			}
			else
			{
				badLoops.push_back(eIt->first);
				// cout<<"add bad loop "<<eIt->first.first<<" "<<eIt->first.second<<endl;
			}
			// cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
			// 	" has error "<< eIt->second<<endl;
		}

		// if(goodLoop == 0)
		// {
		// 	cout<<"good loop size = 0, abandon this cluster"<<endl;
		// 	return 0;
		// }

		//handle the bad loopps, through check their transform distance to good ones	
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		std::array<std::array<double, 5>, 5>  length;
		std::vector<double> disVector;
		std::pair<bool, double> reV;
		bool futileBit1, futileBit2, futileBit;
		double  transX_residual, transY_residual, transA_residual, covX, covY, beli = 0.95;
		int badNum, deletedNum = 0;
		bool conflict = 0, have_deleted ;
		std::vector<std::pair<int, int>> good_loop_node_distance_pair_vector;
		std::pair<int, int> ele_local;
		int serial_good;
		while(1)
		{
			badNum = badLoops.size() - 1;
			deletedNum = 0;
			int testTime;
			double average;
			for( ; badNum>=0; badNum--)
			{
				testTime = 9;
				disVector.clear();
				if(testTime > int(goodLoops.size()) -1){
					testTime = int(goodLoops.size()) -1;
				}
				// testTime = testTime-1;
				// cout<<"goodLoops.size() "<<goodLoops.size()<<endl;
				// cout<<"goodLoops.size() -1:  "<<goodLoops.size()-1<<endl;

				//sort the good loop based on node distacne to the test bad loop	
				good_loop_node_distance_pair_vector.clear();
				

				if(goodLoops.size() >= 1)
				{
					for(int local_i = 0; local_i < goodLoops.size(); local_i++)
					{
						ele_local.first = local_i;
						ele_local.second = abs(badLoops[badNum].first - goodLoops[local_i].first) + 
							abs(badLoops[badNum].second - goodLoops[local_i].second);
						good_loop_node_distance_pair_vector.push_back(ele_local);
					}	
					std::sort (good_loop_node_distance_pair_vector.begin(), good_loop_node_distance_pair_vector.end(), cmp_reverse);	
				}
				
				cout<<"testTime "<<testTime<<endl;

				cout<<" "<<endl;
				// cout<<"goodLoops size "<<goodLoops.size()<<" "<<" badLoops size "<<badLoops.size()<<endl;
				cout<<"bad loop "<<badNum<<"   "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<
					" has robust error "<<chi2LinkErrors[badLoops[badNum]].first<<"  and non-robust error "<<chi2LinkErrors[badLoops[badNum]].second <<endl;
				// beli = 1 - utils::chi2_p(3, chi2LinkErrors[badLoops[badNum]]);
				// cout<<"beli for this loop is "<<beli<<endl;
				if(testTime == -1)
				{
					// no good loop so do self check
					reV =  clusterizer.check_single_loop_odo(badLoops[badNum], clusterizer.LP_Trans_Covar_Map, futileBit1, beli);
					if((futileBit1) != 1)
					{
						if((reV.first == 1))
						{
							cout<<" * the dis is "<<reV.second<<endl;
							disVector.push_back(reV.second);
						}		
						else
						{
							cout<<"the dis is "<<reV.second<<endl;
							disVector.push_back(reV.second);
					
						}					
					}
					else{
						cout<<"futile dis is "<<reV.second<<endl;
					}					
				}
				else
				{
					bool done_single_check = 0;
					have_deleted = 0;
					int oritianl_test_times = testTime;
					for(; testTime >=0; testTime--){
						// cout<<"testTime "<<testTime<<endl;
						//if the error is too big, delete it directly
						if(chi2LinkErrors[badLoops[badNum]].first > utils::chi2_continuous(3, 0.95)){
							cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
							clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
							deleted_loops.push_back(badLoops[badNum]);
							deletedNum = deletedNum + 1;	
							have_deleted = 1;
							break;
						}
						//
						serial_good = good_loop_node_distance_pair_vector[oritianl_test_times - testTime].first;
						cout<<"god loop "<<goodLoops[serial_good].first<<" "<<goodLoops[serial_good].second<<" in cluster "<<
							clusterizer.getClusterID(goodLoops[serial_good])<<
							" with node dis "<<good_loop_node_distance_pair_vector[testTime].second<<endl;

						if((good_loop_node_distance_pair_vector[testTime].second > abs(badLoops[badNum].first - badLoops[badNum].second)) and
							done_single_check == 0)
						{
							cout<<"do single loop check"<<endl;

							reV =  clusterizer.check_single_loop_odo(badLoops[badNum], clusterizer.LP_Trans_Covar_Map, futileBit1, beli);//, double& statis
							done_single_check = 1;
							if((futileBit1) != 1)
							{
								if((reV.first == 1))
								{
									cout<<" * the dis is "<<reV.second<<endl;
									disVector.push_back(reV.second);
								}		
								else
								{
									cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
									clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
									deleted_loops.push_back(badLoops[badNum]);
									deletedNum = deletedNum + 1;
									cout<<"the dis is "<<reV.second<<endl;
									disVector.push_back(reV.second);
									have_deleted = 1;
									cout<<"self check result in delete"<<endl;
									cin.get();
									break;
							
								}					
							}
							else{
								cout<<"futile dis is "<<reV.second<<endl;
								cout<<"this segment is the smallest distance but it is futile, so we do not do further check"<<endl;
								break;
							}
						}
						else
						{
							clusterizer.prepare_to_signle_loop_pair_check(badLoops[badNum], goodLoops[serial_good], 
								FullInfo,  clusterizer.LP_Trans_Covar_Map, length, futileBit1, futileBit2);
							// cout<<"%95 : "<<utils::chi2_p(3, utils::chi2_continuous(3, 0.99) )<<endl;
							// cin.get();

							if((futileBit2 or futileBit1) != 1)
							{
								reV =  clusterizer.check_single_loop_inter_varying_belief(FullInfo, covX, covY, clusterizer.displayCov, 
									transX_residual, transY_residual, transA_residual, length[4][0], futileBit, beli);
								if((futileBit == 0))
								{
									if((reV.first == 1))
									{
										cout<<" * the dis is "<<reV.second<<endl;
										disVector.push_back(reV.second);
									}		
									else
									{
										cout<<"the dis is "<<reV.second<<endl;
										disVector.push_back(reV.second);
								
									}					
								}
								else{
									cout<<"futile dis is "<<reV.second<<endl;
								}
							}
							else
							{
								cout<<"futile dis is "<<reV.second<<endl;
								cout<<"current segment is futile, so we do not do further check"<<endl;
								break;
							}
						}
					}
				}
				if(have_deleted == 0)
				{
					// average = accumulate( disVector.begin(), disVector.end(), 0.0)/disVector.size();  				
				
					if(disVector.size() == 0)
					{
						// cout<<"no evidence can prove the genuity of this doublt loop, so delete it"<<endl;
						// cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
						// clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
						// deleted_loops.push_back(badLoops[badNum]);
						// deletedNum = deletedNum + 1;		
						if(chi2LinkErrors[badLoops[badNum]].first > utils::chi2_continuous(3, 0.5))
						{
							cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
							clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
							deleted_loops.push_back(badLoops[badNum]);
							deletedNum = deletedNum + 1;
						}
						cout<<"no valid nearby loop to check it, but error larger than %50"<<endl;		
					} 
					else
					{
						// average = *(std::max_element(disVector.begin(), disVector.end()));
						average = accumulate( disVector.begin(), disVector.end(), 0.0)/disVector.size();  
						cout<<"max dis is "<<average<<" threshold is "<<utils::chi2_continuous(3, beli)<<endl;
						if(average > utils::chi2_continuous(3, beli) )
						{
							cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
							clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
							deleted_loops.push_back(badLoops[badNum]);
							deletedNum = deletedNum + 1;				
						} 
					} 
				}	
			}
			  
			if(deletedNum == 0 )
			{
				currentCluster = clusterizer.getClusterByID(clusterID);
				double current_size = currentCluster.size();
				if(current_size == 0)
				{
					cout<<"all loops are deleted "<<endl;
					return 0;
				}

				//the third parameter should be set to negative things, it is used to skip some certain loop when debug
				clusterNodeSet.clear();
				findRange = currentCluster.begin();
				findRangeEnd = currentCluster.end();
				for(; findRange != findRangeEnd; findRange++)
				{
					clusterNodeSet.insert((*findRange).first);
					clusterNodeSet.insert((*findRange).second);
				}
				startRange = *clusterNodeSet.begin();
				endRange   = *clusterNodeSet.rbegin();

				chi2LinkErrors.clear();
				// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
				chi2LinkErrors = gWrapper->optimize_active_robustKernel_pair_return(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);

				activeEdgeCount = int(chi2LinkErrors[IntPair(-1,-1)].first);
				if(chi2LinkErrors[IntPair(-1,0)].first < utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95) and
					chi2LinkErrors[IntPair(-2,0)].first < utils::chi2_continuous(edgeDimension* int(chi2LinkErrors[IntPair(-2,-1)].first), 0.95))
				{
					std::cerr<<" Cluster "<<clusterID <<" survived with "<<
						clusterizer.getClusterByID(clusterID).size()<<" links "<<"from "<<realOrigianlNum<<std::endl;
					cout<<"activeChi2Graph: "<<chi2LinkErrors[IntPair(-1,0)].first <<" activeEdgeCount: "<<activeEdgeCount<<endl;
					cout<<"clsuter error: "<<chi2LinkErrors[IntPair(-2,0)].first <<" activeEdgeCount: "<<activeEdgeCount<<endl;
					// cout<<"cluster size > 1, and some menbers are deleted"<<endl;
					// std::cin.get();
					// if()
					return 1;
				}	
				else{
					std::cerr<<" Cluster "<<clusterID <<" being deleted with "<<realOrigianlNum<<" loops"<<std::endl;
					return 0;
				}		
			}
			else if(deletedNum == originalTotalNum or (clusterizer.getClusterByID(clusterID)).size() == 0)
			{
				cout<<"all loops are deleted "<<endl;
				return 0;
			}

			currentCluster = clusterizer.getClusterByID(clusterID);
			if(currentCluster.size() == 0)
			{
				cout<<"all loops are deleted"<<endl;
				return 0;
			}

			//the third parameter should be set to negative things, it is used to skip some certain loop when debug
			goodLoop = 0;
			// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			clusterNodeSet.clear();
			findRange = currentCluster.begin();
			findRangeEnd = currentCluster.end();
			for(; findRange != findRangeEnd; findRange++)
			{
				clusterNodeSet.insert((*findRange).first);
				clusterNodeSet.insert((*findRange).second);
			}
			startRange = *clusterNodeSet.begin();
			endRange   = *clusterNodeSet.rbegin();

			chi2LinkErrors.clear();
			// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			chi2LinkErrors = gWrapper->optimize_active_robustKernel_pair_return(currentCluster, small_numIterations, 
									  startRange, endRange, odoEdgeError);

			activeChi2Graph = chi2LinkErrors[IntPair(-1,0)].first;
			activeEdgeCount = int(chi2LinkErrors[IntPair(-1,-1)].second);
			eIt = chi2LinkErrors.begin();
			eEnd = chi2LinkErrors.end();
			goodLoops.clear();
			badLoops.clear();
			for(; eIt!=eEnd; eIt++)
			{
				if(eIt->first.first < 0)
					continue;

				// if(eIt->second < chi2_continuous(3, 0.5))
				// {
				// 	real_good_loops.push_back(eIt->first);
				// }
				if(eIt->second.second < utils::chi2_continuous(3, 0.95))	
				{
					goodLoop=goodLoop+1;
					goodLoops.push_back(eIt->first);
				}
				else
				{
					badLoops.push_back(eIt->first);
				}
				cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
					" has robust error "<< eIt->second.first<<" and non-robust err "<< eIt->second.second<<endl;
			}
			cout<<"bad loop size:"<<badLoops.size()<<" good size "<<goodLoop<<endl;
			cout<<"activeChi2Graph is "<<activeChi2Graph<<endl;
			if(badLoops.size() == 0)
			{
				if(activeChi2Graph > utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95))
				{
					std::cerr<<"deleted Cluster "<<" have "<<realOrigianlNum<<" links "<<std::endl;
					return 0;
				}
				std::cerr<<" Cluster "<<clusterID <<"survived with "<<
					clusterizer.getClusterByID(clusterID).size()<<" links "<<"from "<<realOrigianlNum<<std::endl;
				cout<<"cluster size > 1, and some menbers are deleted"<<endl;

					cout<<"activeChi2Graph: "<<chi2LinkErrors[IntPair(-1,0)].first <<" activeEdgeCount: "<<activeEdgeCount<<endl;
					cout<<"clsuter error: "<<chi2LinkErrors[IntPair(-2,0)].first <<" activeEdgeCount: "<<activeEdgeCount<<endl;	
								
				// std::cin.get();
				return 1;
			}
			if(goodLoops.size() == 0)
			{
				std::cerr<<"deleted Cluster "<<" have "<<realOrigianlNum<<" links "<<std::endl;
				return 0;
			}			
		}
	}	
	bool complementary_intra_second(int clusterID, std::vector<double> & chiStatisVector, std::vector<std::pair<int, int> > & doubtLoopNumberVector, 
		int & originalLoopNum, std::vector<std::pair<double, std::pair<int, int> > > & chiStstis4eachLoop,
		std::vector<std::pair<int, int> > & deleted_loops)
	{
		IntPairDoubleMap chi2LinkErrors;
		IntPairSet& currentCluster = clusterizer.getClusterByID(clusterID);

		std::vector<std::pair<int,int> > goodLoops, badLoops,real_good_loops;
		int originalTotalNum = currentCluster.size();
		cout<<"cluster "<<clusterID<<" has "<<originalTotalNum<<" loops"<<endl;
		originalLoopNum = originalTotalNum;
		int realOrigianlNum = originalLoopNum;

		//the third parameter should be set to negative things, it is used to skip some certain loop when debug
		int startRange, endRange;
		std::set<int> clusterNodeSet;
		IntPairSet::iterator findRange = currentCluster.begin(), findRangeEnd = currentCluster.end();
		for(; findRange != findRangeEnd; findRange++)
		{
			clusterNodeSet.insert((*findRange).first);
			clusterNodeSet.insert((*findRange).second);
		}
		startRange = *clusterNodeSet.begin();
		endRange   = *clusterNodeSet.rbegin();
        cout<<"node range from "<<startRange<<" to "<<endRange<<endl;

        // if originalTotalNum = 1, 
        string ba = "sdf", aa = "aa";

        std::vector< g2o::HyperGraph::Edge* > newEdgeVector;
        if(originalTotalNum == 1)
        {
        	std::vector<std::pair<int, int>> artificialLoops;
        	int maxVertexID = gWrapper->optimizer->vertices().size();
        	int futileDis;
        	if(clusterizer.futilePairVector.size() != 0)
        		futileDis = abs(clusterizer.futilePairVector[0].first - clusterizer.futilePairVector[0].second);
        	else
        		futileDis = maxVertexID;
        	
        	// generateNewEdge(*(currentCluster.begin()), newEdgeVector, artificialLoops, futileDis, maxVertexID);
			chi2LinkErrors =  
				// gWrapper->optimize_oneMember(currentCluster, small_numIterations, newEdgeVector, artificialLoops);
				gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);
        }
        else
        {
			chi2LinkErrors =  		
				gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);
        }

		float activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
		int   activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
		cout<<"active edge count: "<<activeEdgeCount<<endl;

					cout<<"activeChi2Graph: "<<chi2LinkErrors[IntPair(-1,0)] <<" activeEdgeCount: "<<activeEdgeCount<<endl;
					cout<<"clsuter error: "<<chi2LinkErrors[IntPair(-2,0)] <<" originalLoopNum: "<<originalLoopNum<<endl;


		//check if there is someone loop whose chi2 error is bigger than the signgle edge limit error
		double biggestErr = 0, smallcHI = 200;
		std::pair<int, int> biggestErrorLoop(-1, -1);

		cout<<" "<<endl;
		IntPairDoubleMap::iterator eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		double goodLoop = 0;
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;
			if(eIt->second > biggestErr)
			{
				biggestErr = eIt->second;
				biggestErrorLoop = eIt->first;
			}
			// if(eIt->second < utils::chi2_continuous(3, 0.5))
			// {
			// 	real_good_loops.push_back(eIt->first);
			// }
			if(eIt->second < utils::chi2_continuous(3, 0.5))	
			{
				goodLoop=goodLoop+1;
				goodLoops.push_back(eIt->first);
			}
			else
			{
				badLoops.push_back(eIt->first);
				// cout<<"add bad loop "<<eIt->first.first<<" "<<eIt->first.second<<endl;
			}
			// cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
			// 	" has error "<< eIt->second<<endl;
		}


		//handle the bad loopps, through check their transform distance to good ones	
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		std::array<std::array<double, 5>, 5>  length;
		std::vector<double> disVector;
		std::pair<bool, double> reV;
		bool futileBit1, futileBit2, futileBit;
		double  transX_residual, transY_residual, transA_residual, covX, covY, beli = 0.95;
		int badNum, deletedNum = 0;
		bool conflict = 0, have_deleted ;
		std::vector<std::pair<int, int>> good_loop_node_distance_pair_vector;
		std::pair<int, int> ele_local;
		int serial_good;
		while(1)
		{
			badNum = badLoops.size() - 1;
			deletedNum = 0;
			int testTime;
			double average;
			for( ; badNum>=0; badNum--)
			{
				testTime = 9;
				disVector.clear();
				if(testTime > int(goodLoops.size()) -1){
					testTime = int(goodLoops.size()) -1;
				}
				// testTime = testTime-1;
				// cout<<"goodLoops.size() "<<goodLoops.size()<<endl;
				// cout<<"goodLoops.size() -1:  "<<goodLoops.size()-1<<endl;

				//sort the good loop based on node distacne to the test bad loop	
				good_loop_node_distance_pair_vector.clear();
				

				if(goodLoops.size() >= 1)
				{
					for(int local_i = 0; local_i < goodLoops.size(); local_i++)
					{
						ele_local.first = local_i;
						ele_local.second = abs(badLoops[badNum].first - goodLoops[local_i].first) + 
							abs(badLoops[badNum].second - goodLoops[local_i].second);
						good_loop_node_distance_pair_vector.push_back(ele_local);
					}	
					std::sort (good_loop_node_distance_pair_vector.begin(), good_loop_node_distance_pair_vector.end(), cmp_reverse);	
				}
		
				cout<<"testTime "<<testTime<<endl;

				cout<<" "<<endl;
				// cout<<"goodLoops size "<<goodLoops.size()<<" "<<" badLoops size "<<badLoops.size()<<endl;
				cout<<"bad loop "<<badNum<<"   "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<
					" has error "<<chi2LinkErrors[badLoops[badNum]]<<endl;
				// beli = 1 - utils::chi2_p(3, chi2LinkErrors[badLoops[badNum]]);
				cout<<"beli for this loop is "<<beli<<endl;
				if(testTime == -1)
				{
					// no good loop so do self check
					reV =  clusterizer.check_single_loop_odo(badLoops[badNum], clusterizer.LP_Trans_Covar_Map, futileBit1, beli);
					if((futileBit1) != 1)
					{
						if((reV.first == 1))
						{
							cout<<" * the dis is "<<reV.second<<endl;
							disVector.push_back(reV.second);
						}		
						else
						{
							cout<<"the dis is "<<reV.second<<endl;
							disVector.push_back(reV.second);
					
						}					
					}
					else{
						cout<<"futile dis is "<<reV.second<<endl;
					}					
				}
				else
				{
					bool done_single_check = 0;
					have_deleted = 0;
					for(; testTime >=0; testTime--){
						// cout<<"testTime "<<testTime<<endl;
						//if the error is too big, delete it directly
						if(chi2LinkErrors[badLoops[badNum]] > utils::chi2_continuous(3, 0.95)){
							cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
							clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
							deleted_loops.push_back(badLoops[badNum]);
							deletedNum = deletedNum + 1;	
							have_deleted = 1;
							break;
						}
						//
						serial_good = good_loop_node_distance_pair_vector[testTime].first;
						cout<<"god loop "<<goodLoops[serial_good].first<<" "<<goodLoops[serial_good].second<<
							" with node dis "<<good_loop_node_distance_pair_vector[testTime].second<<endl;

						if((good_loop_node_distance_pair_vector[testTime].second > abs(badLoops[badNum].first - badLoops[badNum].second)) and
							done_single_check == 0)
						{
							cout<<"do single loop check"<<endl;

							reV =  clusterizer.check_single_loop_odo(badLoops[badNum], clusterizer.LP_Trans_Covar_Map, futileBit1, beli);//, double& statis
							done_single_check = 1;
							if((futileBit1) != 1)
							{
								if((reV.first == 1))
								{
									cout<<" * the dis is "<<reV.second<<endl;
									disVector.push_back(reV.second);
								}		
								else
								{
									cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
									clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
									deleted_loops.push_back(badLoops[badNum]);
									deletedNum = deletedNum + 1;
									cout<<"the dis is "<<reV.second<<endl;
									disVector.push_back(reV.second);
									have_deleted = 1;
									break;
							
								}					
							}
							else{
								cout<<"futile dis is "<<reV.second<<endl;
							}
						}
						else
						{
							clusterizer.prepare_to_signle_loop_pair_check(badLoops[badNum], goodLoops[serial_good], 
								FullInfo,  clusterizer.LP_Trans_Covar_Map, length, futileBit1, futileBit2);
							// cout<<"%95 : "<<utils::chi2_p(3, utils::chi2_continuous(3, 0.99) )<<endl;
							// cin.get();

							if((futileBit2 or futileBit1) != 1){
								reV =  clusterizer.check_single_loop_inter_varying_belief(FullInfo, covX, covY, clusterizer.displayCov, 
									transX_residual, transY_residual, transA_residual, length[4][0], futileBit, beli);
								if((futileBit == 0))
								{
									if((reV.first == 1))
									{
										cout<<" * the dis is "<<reV.second<<endl;
										disVector.push_back(reV.second);
									}		
									else
									{
										cout<<"the dis is "<<reV.second<<endl;
										disVector.push_back(reV.second);
								
									}					
								}
								else{
									cout<<"futile dis is "<<reV.second<<endl;
								}
							}
						}
					}
				}
				if(have_deleted == 0)
				{
					average = accumulate( disVector.begin(), disVector.end(), 0.0)/disVector.size();  
					cout<<"averange dis is "<<average<<" threshold is "<<utils::chi2_continuous(3, beli)<<endl;
					if(average > utils::chi2_continuous(3, beli))
					{
						cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
						clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
						deleted_loops.push_back(badLoops[badNum]);
						deletedNum = deletedNum + 1;				
					}  	
				}
				
			}
			  
			if(deletedNum == 0 )
			{
				currentCluster = clusterizer.getClusterByID(clusterID);
				double current_size = currentCluster.size();
				if(current_size == 0)
				{
					cout<<"all loops are deleted "<<endl;
					return 0;
				}

				//the third parameter should be set to negative things, it is used to skip some certain loop when debug
				clusterNodeSet.clear();
				findRange = currentCluster.begin();
				findRangeEnd = currentCluster.end();
				for(; findRange != findRangeEnd; findRange++)
				{
					clusterNodeSet.insert((*findRange).first);
					clusterNodeSet.insert((*findRange).second);
				}
				startRange = *clusterNodeSet.begin();
				endRange   = *clusterNodeSet.rbegin();

				chi2LinkErrors.clear();
				// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
				chi2LinkErrors = gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);

				activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
				if(chi2LinkErrors[IntPair(-1,0)] < utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95) and
					chi2LinkErrors[IntPair(-2,0)] < utils::chi2_continuous(edgeDimension* chi2LinkErrors[IntPair(-2,-1)], 0.95))
				{
					std::cerr<<" Cluster "<<clusterID <<" survived with "<<
						clusterizer.getClusterByID(clusterID).size()<<" links "<<"from "<<realOrigianlNum<<std::endl;
					cout<<"activeChi2Graph: "<<chi2LinkErrors[IntPair(-1,0)] <<" activeEdgeCount: "<<activeEdgeCount<<endl;
					cout<<"clsuter error: "<<chi2LinkErrors[IntPair(-2,0)] <<" activeEdgeCount: "<<activeEdgeCount<<endl;
					// cout<<"cluster size > 1, and some menbers are deleted"<<endl;
					// std::cin.get();
					// if()
					return 1;
				}	
				else{
					std::cerr<<" Cluster "<<clusterID <<" being deleted with "<<realOrigianlNum<<" loops"<<std::endl;
					return 0;
				}		
			}
			else if(deletedNum == originalTotalNum or (clusterizer.getClusterByID(clusterID)).size() == 0)
			{
				cout<<"all loops are deleted "<<endl;
				return 0;
			}

			currentCluster = clusterizer.getClusterByID(clusterID);
			if(currentCluster.size() == 0)
			{
				cout<<"all loops are deleted"<<endl;
				return 0;
			}

			//the third parameter should be set to negative things, it is used to skip some certain loop when debug
			goodLoop = 0;
			// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			clusterNodeSet.clear();
			findRange = currentCluster.begin();
			findRangeEnd = currentCluster.end();
			for(; findRange != findRangeEnd; findRange++)
			{
				clusterNodeSet.insert((*findRange).first);
				clusterNodeSet.insert((*findRange).second);
			}
			startRange = *clusterNodeSet.begin();
			endRange   = *clusterNodeSet.rbegin();

			chi2LinkErrors.clear();
			// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			chi2LinkErrors = gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
									  startRange, endRange, odoEdgeError);

			activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
			activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
			eIt = chi2LinkErrors.begin();
			eEnd = chi2LinkErrors.end();
			biggestErr = 0;
			goodLoops.clear();
			badLoops.clear();
			for(; eIt!=eEnd; eIt++)
			{
				if(eIt->first.first < 0)
					continue;
				if(eIt->second > biggestErr)
				{
					biggestErr = eIt->second;
					biggestErrorLoop = eIt->first;
				}

				// if(eIt->second < chi2_continuous(3, 0.5))
				// {
				// 	real_good_loops.push_back(eIt->first);
				// }
				if(eIt->second < utils::chi2_continuous(3, 0.5))	
				{
					goodLoop=goodLoop+1;
					goodLoops.push_back(eIt->first);
				}
				else
				{
					badLoops.push_back(eIt->first);
				}
				cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
					" has error "<< eIt->second<<endl;
			}
			cout<<"bad loop size:"<<badLoops.size()<<" good size "<<goodLoop<<endl;
			cout<<"activeChi2Graph is "<<activeChi2Graph<<endl;
			if(badLoops.size() == 0)
			{
				if(activeChi2Graph > utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95))
				{
					std::cerr<<"deleted Cluster "<<" have "<<realOrigianlNum<<" links "<<std::endl;
					return 0;
				}
				std::cerr<<" Cluster "<<clusterID <<"survived with "<<
					clusterizer.getClusterByID(clusterID).size()<<" links "<<"from "<<realOrigianlNum<<std::endl;
				cout<<"cluster size > 1, and some menbers are deleted"<<endl;

					cout<<"activeChi2Graph: "<<chi2LinkErrors[IntPair(-1,0)] <<" activeEdgeCount: "<<activeEdgeCount<<endl;
					cout<<"clsuter error: "<<chi2LinkErrors[IntPair(-2,0)] <<" activeEdgeCount: "<<activeEdgeCount<<endl;	
								
				// std::cin.get();
				return 1;
			}
			if(goodLoops.size() == 0)
			{
				std::cerr<<"deleted Cluster "<<" have "<<realOrigianlNum<<" links "<<std::endl;
				return 0;
			}			
		}
	}
	bool complementary_intra(int clusterID, std::vector<double> & chiStatisVector, std::vector<std::pair<int, int> > & doubtLoopNumberVector, 
		int & originalLoopNum, std::vector<std::pair<double, std::pair<int, int> > > & chiStstis4eachLoop,
		std::vector<std::pair<int, int> > & deleted_loops)
	{
		IntPairDoubleMap chi2LinkErrors;
		IntPairSet& currentCluster = clusterizer.getClusterByID(clusterID);

		std::vector<std::pair<int,int> > goodLoops, badLoops;
		int originalTotalNum = currentCluster.size();
		cout<<"cluster "<<clusterID<<" has "<<originalTotalNum<<" loops"<<endl;
		originalLoopNum = originalTotalNum;
		int realOrigianlNum = originalLoopNum;

		//the third parameter should be set to negative things, it is used to skip some certain loop when debug
		int startRange, endRange;
		std::set<int> clusterNodeSet;
		IntPairSet::iterator findRange = currentCluster.begin(), findRangeEnd = currentCluster.end();
		for(; findRange != findRangeEnd; findRange++)
		{
			clusterNodeSet.insert((*findRange).first);
			clusterNodeSet.insert((*findRange).second);
		}
		startRange = *clusterNodeSet.begin();
		endRange   = *clusterNodeSet.rbegin();
        cout<<"node range from "<<startRange<<" to "<<endRange<<endl;

        // if originalTotalNum = 1, 
        string ba = "sdf", aa = "aa";

        std::vector< g2o::HyperGraph::Edge* > newEdgeVector;
   //      if(originalTotalNum == 1)
   //      {
        	// std::vector<std::pair<int, int>> artificialLoops;
        	// int maxVertexID = gWrapper->optimizer->vertices().size();
        	// int futileDis;
        	// if(clusterizer.futilePairVector.size() != 0)
        	// 	futileDis = abs(clusterizer.futilePairVector[0].first - clusterizer.futilePairVector[0].second);
        	// else
        	// 	futileDis = maxVertexID;
        	
        	// generateNewEdge(*(currentCluster.begin()), newEdgeVector, artificialLoops, futileDis, maxVertexID);
			// chi2LinkErrors =  
			// 	gWrapper->optimize_oneMember(currentCluster, small_numIterations, newEdgeVector, artificialLoops);
			// 	gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
			// 							  startRange, endRange, odoEdgeError);
   //      }
   //      else
        cout<<"run complementary intra"<<endl;
        {
			chi2LinkErrors =  		
				// gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
				// 						  startRange, endRange, odoEdgeError, newEdgeVector, artificialLoops);
				gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);				
				// gWrapper->optimize_robustKernel(currentCluster, small_numIterations, odoEdgeError);
        }

		float activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
		int   activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)], loopCount;
		cout<<"active edge count: "<<activeEdgeCount<<endl;
		//check if there is someone loop whose chi2 error is bigger than the signgle edge limit error
		double biggestErr = 0, smallcHI = 200;
		std::pair<int, int> biggestErrorLoop(-1, -1);

		cout<<" "<<endl;
		IntPairDoubleMap::iterator eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		double goodLoop = 0;
		std::vector<std::pair<int, int>> good_loop_node_distance_pair_vector;
		std::pair<int, int> ele_local;
		int serial_good;

		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;
			if(eIt->second > biggestErr)
			{
				biggestErr = eIt->second;
				biggestErrorLoop = eIt->first;
			}
			if(eIt->second < utils::chi2(edgeDimension))		
			{
				goodLoop=goodLoop+1;
				goodLoops.push_back(eIt->first);
			}
			else
			{
				badLoops.push_back(eIt->first);
				// cout<<"add bad loop "<<eIt->first.first<<" "<<eIt->first.second<<endl;
			}
			// cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
			// 	" has error "<< eIt->second<<endl;
		}

		// if(clusterID == 85 and clusterizer.loopToClusterIDMap[std::pair<int,int> (9016, 3876)] != 85)
		// 	exit(0);


		// if(1 == originalTotalNum )
		// {
		// 	if(activeChi2Graph > utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95) or 
		// 		// badLoops.size() != 0)
		// 	// if(activeChi2Graph > 500 or 
		// 		chi2LinkErrors[IntPair(-2,0)] > utils::chi2(edgeDimension * (newEdgeVector.size()+1))
		// 		or badLoops.size() != 0)			
		// 	{
		// 		cout<<"this cluster has one loop, and the activeChi2Graph is bigger than the overall limit, so delete this cluster"<<endl;
		// 		cout<<"activeChi2Graph: "<<activeChi2Graph <<" activeEdgeCount: "<<activeEdgeCount<<" threshold: "<<
		// 			utils::chi2_continuous(edgeDimension* floor(activeEdgeCount), 0.1) <<endl;
		// 		cout<<"edgeCount: "<<activeEdgeCount<<" start range: "<<startRange <<" end range: "<<endRange<<" "<<endRange- startRange<<endl;
		// 		cout<<"gWrapper->odoSize: "<<gWrapper->odoSize<<endl;		
		// 		std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
		// 		return false;
		// 	}	
		// 	else{
		// 		cout<<"activeChi2Graph: "<<activeChi2Graph <<" activeEdgeCount: "<<activeEdgeCount<<endl;
		// 		cout<<"loop error "<<chi2LinkErrors[IntPair(-2,0)]<<endl;
		// 		cout<<"suvived , one menber clsuetr"<<endl;
		// 		// cin.get();
		// 		return true;
		// 	}		
		// } 

		// if(goodLoop == originalTotalNum)
		// {
		// 	if(activeChi2Graph > utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95))
		// 	{
		// 		cout<<"activeChi2Graph: "<<activeChi2Graph<<" is bigger than threshold "<<endl;
		// 		cout<<"empty cluster "<<clusterID<<" has "<<originalTotalNum<<" loop "<<endl;
		// 		return 0;
		// 	}
		// 	else
		// 	{
		// 		cout<<"activeChi2Graph: "<<activeChi2Graph<<", all loops are good"<<endl;
		// 		cout<<"suvived cluster "<<clusterID<<" has "<<originalTotalNum<<" loop "<<endl;
		// 		return 1;
		// 	}
		// }
		// else if(goodLoop == 0)
		// {
		// 	cout<<"activeChi2Graph: "<<activeChi2Graph<<" is bigger than threshold "<<endl;
		// 	cout<<"empty cluster "<<clusterID<<" has "<<originalTotalNum<<" loop "<<endl;
		// 	return 0;
		// }


		//handle the bad loopps, through check their transform distance to good ones	
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		std::array<std::array<double, 5>, 5>  length;
		std::vector<double> disVector;
		std::pair<bool, double> reV;
		bool futileBit1, futileBit2, futileBit;
		double  transX_residual, transY_residual, transA_residual, covX, covY, beli = 0.95, average;
		int badNum, deletedNum = 0;
		bool conflict = 0, have_deleted;
		while(1)
		{
			badNum = badLoops.size() - 1;
			deletedNum = 0;
			int testTime;
			for( ; badNum>=0; badNum--)
			{
				testTime = 9;
				disVector.clear();
				if(testTime > int(goodLoops.size()) -1){
					testTime = int(goodLoops.size()) -1;
				}
				// testTime = testTime-1;
				// cout<<"goodLoops.size() "<<goodLoops.size()<<endl;
				// cout<<"goodLoops.size() -1:  "<<goodLoops.size()-1<<endl;

				if(chi2LinkErrors[badLoops[badNum]] > utils::chi2_continuous(3, 0.95)){
					cout<<"error > % 0.95, direct delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
					clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
					deleted_loops.push_back(badLoops[badNum]);
					deletedNum = deletedNum + 1;	
					have_deleted = 1;
					break;
				}

				cout<<"chech if the program can run the next code"<<endl;
				cin.get();

				good_loop_node_distance_pair_vector.clear();
				
				if(goodLoops.size() >= 1)
				{
					for(int local_i = 0; local_i < goodLoops.size(); local_i++)
					{
						ele_local.first = local_i;
						ele_local.second = abs(badLoops[badNum].first - goodLoops[local_i].first) + 
							abs(badLoops[badNum].second - goodLoops[local_i].second);
						good_loop_node_distance_pair_vector.push_back(ele_local);
					}	
					std::sort (good_loop_node_distance_pair_vector.begin(), good_loop_node_distance_pair_vector.end(), cmp_reverse);	
				}


				cout<<" "<<endl;
				// cout<<"goodLoops size "<<goodLoops.size()<<" "<<" badLoops size "<<badLoops.size()<<endl;
				cout<<"bad loop "<<badNum<<"   "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<
					" has error "<<chi2LinkErrors[badLoops[badNum]]<<endl;

				have_deleted = 0;
				if(testTime == -1)
				{
					cout<<"no good loop to make a pair to test doubt loops"<<endl;
					// no good loop so do self check
					reV =  clusterizer.check_single_loop_odo(badLoops[badNum], clusterizer.LP_Trans_Covar_Map, futileBit1, beli);
					if((futileBit1) != 1)
					{
						if((reV.first == 1))
						{
							cout<<" * the dis is "<<reV.second<<endl;
							disVector.push_back(reV.second);
						}		
						else
						{
							cout<<"the dis is "<<reV.second<<endl;
							disVector.push_back(reV.second);
					
						}					
					}
					else{
						cout<<"futile dis is "<<reV.second<<endl;
					}					
				}
				else
				{
					bool done_single_check = 0;
					
					for(; testTime >=0; testTime--){
						// cout<<"testTime "<<testTime<<endl;
						//if the error is too big, delete it directly

						//
						serial_good = good_loop_node_distance_pair_vector[testTime].first;
						cout<<"god loop "<<goodLoops[serial_good].first<<" "<<goodLoops[serial_good].second<<
							" with node dis "<<good_loop_node_distance_pair_vector[testTime].second<<endl;

						if((good_loop_node_distance_pair_vector[testTime].second > abs(badLoops[badNum].first - badLoops[badNum].second)) and
							done_single_check == 0)
						{
							cout<<"do single loop check"<<endl;

							reV =  clusterizer.check_single_loop_odo(badLoops[badNum], clusterizer.LP_Trans_Covar_Map, futileBit1, beli);//, double& statis
							done_single_check = 1;
							if((futileBit1) != 1)
							{
								if((reV.first == 1))
								{
									cout<<" * the dis is "<<reV.second<<endl;
									disVector.push_back(reV.second);
								}		
								else
								{
									cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
									clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
									deleted_loops.push_back(badLoops[badNum]);
									deletedNum = deletedNum + 1;
									cout<<"the dis is "<<reV.second<<endl;
									disVector.push_back(reV.second);
									have_deleted = 1;
									break;
							
								}					
							}
							else{
								cout<<"futile dis is "<<reV.second<<endl;
							}
						}
						else
						{
							clusterizer.prepare_to_signle_loop_pair_check(badLoops[badNum], goodLoops[serial_good], 
								FullInfo,  clusterizer.LP_Trans_Covar_Map, length, futileBit1, futileBit2);
							// cout<<"%95 : "<<utils::chi2_p(3, utils::chi2_continuous(3, 0.99) )<<endl;
							// cin.get();

							if((futileBit2 or futileBit1) != 1){
								reV =  clusterizer.check_single_loop_inter_varying_belief(FullInfo, covX, covY, clusterizer.displayCov, 
									transX_residual, transY_residual, transA_residual, length[4][0], futileBit, beli);
								if((futileBit == 0))
								{
									if((reV.first == 1))
									{
										cout<<" * the dis is "<<reV.second<<endl;
										disVector.push_back(reV.second);
									}		
									else
									{
										cout<<"the dis is "<<reV.second<<endl;
										disVector.push_back(reV.second);
								
									}					
								}
								else{
									cout<<"futile dis is "<<reV.second<<endl;
								}
							}
						}
					}
				}
				if(have_deleted == 0)
				{
					if(disVector.size() != 0)
					{
						average = accumulate( disVector.begin(), disVector.end(), 0.0)/disVector.size();  
						cout<<"averange dis is "<<average<<" threshold is "<<utils::chi2_continuous(3, beli)<<endl;
						if(average > utils::chi2_continuous(3, beli))
						{
							cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
							clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
							deleted_loops.push_back(badLoops[badNum]);
							deletedNum = deletedNum + 1;				
						} 
					}
					else
					{
						cout<<"the dis vector size is 0"<<endl;
					}
 	
				}
			}

			if(deletedNum == 0 )
			{
				currentCluster = clusterizer.getClusterByID(clusterID);
				double current_size = currentCluster.size();
				if(current_size == 0)
				{
					cout<<"all loops are deleted "<<endl;
					return 0;
				}

				//the third parameter should be set to negative things, it is used to skip some certain loop when debug
				clusterNodeSet.clear();
				findRange = currentCluster.begin();
				findRangeEnd = currentCluster.end();
				for(; findRange != findRangeEnd; findRange++)
				{
					clusterNodeSet.insert((*findRange).first);
					clusterNodeSet.insert((*findRange).second);
				}
				startRange = *clusterNodeSet.begin();
				endRange   = *clusterNodeSet.rbegin();

				chi2LinkErrors.clear();

	        	std::vector<std::pair<int, int>> artificialLoops;
	        	int maxVertexID = gWrapper->optimizer->vertices().size();
	        	int futileDis;
	        	if(clusterizer.futilePairVector.size() != 0)
	        		futileDis = abs(clusterizer.futilePairVector[0].first - clusterizer.futilePairVector[0].second);
	        	else
	        		futileDis = maxVertexID;
	        	generateNewEdge(*(currentCluster.begin()), newEdgeVector, artificialLoops, futileDis, maxVertexID);
	        	// generateNewEdge_parallel(*(currentCluster.begin()), newEdgeVector, artificialLoops, futileDis, maxVertexID, 0);
	        	if(currentCluster.size() > 1)
	        	{
	        		generateNewEdge(*(currentCluster.rbegin()), newEdgeVector, artificialLoops, futileDis, maxVertexID);
	        		// generateNewEdge_parallel(*(currentCluster.rbegin()), newEdgeVector, artificialLoops, futileDis, maxVertexID, 1);
	        	}
	        	
				// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
				chi2LinkErrors = gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError, newEdgeVector, artificialLoops);	
								// gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
								// 		  startRange, endRange, odoEdgeError);

				activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
				if(chi2LinkErrors[IntPair(-1,0)] < utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95) and
					chi2LinkErrors[IntPair(-2,0)] < utils::chi2_continuous(edgeDimension* chi2LinkErrors[IntPair(-2,-1)], 0.95))
				{
					std::cerr<<" Cluster "<<clusterID <<" survived with "<<
						clusterizer.getClusterByID(clusterID).size()<<" links "<<"from "<<realOrigianlNum<<std::endl;
					// cout<<"cluster size > 1, and some menbers are deleted"<<endl;
					// std::cin.get();
					// if()
					return 1;
				}	
				else{
					std::cerr<<" Cluster "<<clusterID <<" being deleted with "<<realOrigianlNum<<" loops"<<std::endl;
					return 0;
				}		
			}
			else if(deletedNum == originalTotalNum or (clusterizer.getClusterByID(clusterID)).size() == 0)
			{
				cout<<"all loops are deleted "<<endl;
				return 0;
			}

			currentCluster = clusterizer.getClusterByID(clusterID);
			if(currentCluster.size() == 0)
			{
				cout<<"all loops are deleted"<<endl;
				return 0;
			}

			//the third parameter should be set to negative things, it is used to skip some certain loop when debug
			goodLoop = 0;
			// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			clusterNodeSet.clear();
			findRange = currentCluster.begin();
			findRangeEnd = currentCluster.end();
			for(; findRange != findRangeEnd; findRange++)
			{
				clusterNodeSet.insert((*findRange).first);
				clusterNodeSet.insert((*findRange).second);
			}
			startRange = *clusterNodeSet.begin();
			endRange   = *clusterNodeSet.rbegin();

			chi2LinkErrors.clear();

	        	std::vector<std::pair<int, int>> artificialLoops;
	        	int maxVertexID = gWrapper->optimizer->vertices().size();
	        	int futileDis;
	        	if(clusterizer.futilePairVector.size() != 0)
	        		futileDis = abs(clusterizer.futilePairVector[0].first - clusterizer.futilePairVector[0].second);
	        	else
	        		futileDis = maxVertexID;
	        	generateNewEdge(*(currentCluster.begin()), newEdgeVector, artificialLoops, futileDis, maxVertexID);
	        	// generateNewEdge_parallel(*(currentCluster.begin()), newEdgeVector, artificialLoops, futileDis, maxVertexID, 1);
	        	if(currentCluster.size() > 1)
	        	{
	        		generateNewEdge(*(currentCluster.rbegin()), newEdgeVector, artificialLoops, futileDis, maxVertexID);
	        		// generateNewEdge_parallel(*(currentCluster.rbegin()), newEdgeVector, artificialLoops, futileDis, maxVertexID, 1);
	        	}
				// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
				chi2LinkErrors = gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError, newEdgeVector, artificialLoops);	
								// gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
								// 		  startRange, endRange, odoEdgeError);			
			// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			// chi2LinkErrors = gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
			// 						  startRange, endRange, odoEdgeError);

			activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
			activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
			loopCount = chi2LinkErrors[IntPair(-2,-1)];
			eIt = chi2LinkErrors.begin();
			eEnd = chi2LinkErrors.end();
			biggestErr = 0;
			goodLoops.clear();
			badLoops.clear();
			for(; eIt!=eEnd; eIt++)
			{
				if(eIt->first.first < 0)
					continue;
				if(eIt->second > biggestErr)
				{
					biggestErr = eIt->second;
					biggestErrorLoop = eIt->first;
				}

				if(eIt->second < utils::chi2(edgeDimension))		
				{
					goodLoop=goodLoop+1;
					goodLoops.push_back(eIt->first);
				}
				else
				{
					badLoops.push_back(eIt->first);
				}
				cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
					" has error "<< eIt->second<<endl;
			}
			cout<<"bad loop size:"<<badLoops.size()<<" good size "<<goodLoop<<endl;
			cout<<"activeChi2Graph is "<<activeChi2Graph<<endl;
			if(badLoops.size() == 0)
			{
				if(activeChi2Graph > utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95) or 
					chi2LinkErrors[IntPair(-2,0)] > utils::chi2_continuous(edgeDimension* loopCount, 0.95) )
				{
					std::cerr<<"deleted Cluster "<<" have "<<realOrigianlNum<<" links "<<std::endl;
					return 0;
				}
				std::cerr<<" Cluster "<<clusterID <<"survived with "<<
					clusterizer.getClusterByID(clusterID).size()<<" links "<<"from "<<realOrigianlNum<<std::endl;
				cout<<"cluster size > 1, and some menbers are deleted"<<endl;
				// std::cin.get();
				return 1;
			}
			if(goodLoops.size() == 0)
			{
				std::cerr<<"deleted Cluster "<<" have "<<realOrigianlNum<<" links "<<std::endl;
				return 0;
			}			
		}
	}
	bool complementary_intra_multi(IntPairSet & currentCluster, std::vector<double> & chiStatisVector, 
		std::vector<std::pair<int, int> > & doubtLoopNumberVector, 
		int & originalLoopNum, std::vector<std::pair<double, std::pair<int, int> > > & chiStstis4eachLoop, 
		std::set<int> & clusters_set)
	{
		if(clusters_set.size()== 1)
			return 1;
		IntPairDoubleMap chi2LinkErrors;

		std::vector<std::pair<int,int> > goodLoops, badLoops;
		int originalTotalNum = currentCluster.size();
		originalLoopNum = originalTotalNum;
		int realOrigianlNum = originalLoopNum;
		std::vector<std::pair<int, int>> good_loop_node_distance_pair_vector;

		//the third parameter should be set to negative things, it is used to skip some certain loop when debug
		int startRange, endRange, serial_good;
		std::set<int> clusterNodeSet;
		IntPairSet::iterator findRange = currentCluster.begin(), findRangeEnd = currentCluster.end();
		for(; findRange != findRangeEnd; findRange++)
		{
			clusterNodeSet.insert((*findRange).first);
			clusterNodeSet.insert((*findRange).second);
		}
		startRange = *clusterNodeSet.begin();
		endRange   = *clusterNodeSet.rbegin();
        cout<<"node range from "<<startRange<<" to "<<endRange<<endl;

        // if originalTotalNum = 1, 
        string ba = "sdf", aa = "aa";

        std::vector< g2o::HyperGraph::Edge* > newEdgeVector;
        if(originalTotalNum == 1)
        {
        	std::vector<std::pair<int, int>> artificialLoops;
        	int maxVertexID = gWrapper->optimizer->vertices().size();
        	int futileDis;
        	if(clusterizer.futilePairVector.size() != 0)
        		futileDis = abs(clusterizer.futilePairVector[0].first - clusterizer.futilePairVector[0].second);
        	else
        		futileDis = maxVertexID;
        	
        	// generateNewEdge(*(currentCluster.begin()), newEdgeVector, artificialLoops, futileDis, maxVertexID);
			chi2LinkErrors =  
				// gWrapper->optimize_oneMember(currentCluster, small_numIterations, newEdgeVector, artificialLoops);
				// gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
				// 						  startRange, endRange, odoEdgeError);
				gWrapper->optimize_active_robustKernel_multi(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);			
        }
        else
        {
			chi2LinkErrors =  		
				// gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
				// 						  startRange, endRange, odoEdgeError);
				gWrapper->optimize_active_robustKernel_multi(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);				
        }

		float activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
		int   activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
		cout<<"active edge count: "<<activeEdgeCount<<endl;
		//check if there is someone loop whose chi2 error is bigger than the signgle edge limit error
		double biggestErr = 0, smallcHI = 200;
		std::pair<int, int> biggestErrorLoop(-1, -1);

		// find good clusters
		std::map<int, std::pair<int, int> > to_find_good_clusters;
		std::set<int> good_clusters;
		int cluster_id_;
		cout<<" "<<endl;
		IntPairDoubleMap::iterator eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		// collect each clsuter's bad loop
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;
			cluster_id_ = clusterizer.getClusterID(eIt->first);
			if(eIt->second > utils::chi2_continuous(edgeDimension, 0.95))		
			{
				if(to_find_good_clusters.find(cluster_id_) == to_find_good_clusters.end())
				{
					to_find_good_clusters[cluster_id_].second = 1;
					to_find_good_clusters[cluster_id_].first  = clusterizer.getClusterByID(cluster_id_).size();
				}
				else
					to_find_good_clusters[cluster_id_].second = to_find_good_clusters[cluster_id_].second + 1;
			}
			// cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
			// 	" has error "<< eIt->second<<endl;
		}
		// select those clusters whose bad loop propation is less than %50 as good ones
		std::set<int>::iterator eIt_find_good_cluster = clusters_set.begin(), 
			eEnd_find_good_cluster = clusters_set.end();
		// cout<<"clsuters num is "<<clusters_set.size()<<endl;
		// cin.get();
		
		for(; eIt_find_good_cluster!=eEnd_find_good_cluster; eIt_find_good_cluster++)
		{

			cluster_id_ = *eIt_find_good_cluster;
			// cout<<"cluster "<<cluster_id_<<endl;


			if (to_find_good_clusters.find(cluster_id_) == to_find_good_clusters.end())
			{
				// cout<<"add to good cluster set"<<endl;
				good_clusters.insert(cluster_id_);
			}
			else
			{
				// cout<<"all loop num = "<<to_find_good_clusters[cluster_id_].first<<"; bad num = "<<
				// 	to_find_good_clusters[cluster_id_].second<<endl;
				if(to_find_good_clusters[cluster_id_].second / to_find_good_clusters[cluster_id_].first <  0.5)
				{
					// cout<<"add to good cluster set"<<endl;
					good_clusters.insert(cluster_id_);
				}
				else
				{
					// cout<<"add to bad cluster set"<<endl;
				}
			}
		}

		eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		double goodLoop = 0;
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;
			if(eIt->second > biggestErr)
			{
				biggestErr = eIt->second;
				biggestErrorLoop = eIt->first;
			}
			if(eIt->second < utils::chi2_continuous(edgeDimension, 0.95))		
			{
				// good loops should come from good clusters
				if(good_clusters.find(clusterizer.getClusterID(eIt->first)) != good_clusters.end())
				{
					goodLoop=goodLoop+1;
					goodLoops.push_back(eIt->first);
				}
				else
				{
					// cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" from bad cluster "<<
					// 	clusterizer.getClusterID(eIt->first)<<endl;
				}
			}
			else
			{
				badLoops.push_back(eIt->first);
				// cout<<"add bad loop "<<eIt->first.first<<" "<<eIt->first.second<<endl;
			}
			// cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
			// 	" has error "<< eIt->second<<endl;
		}


		//handle the bad loopps, through check their transform distance to good ones	
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		std::array<std::array<double, 5>, 5>  length;
		std::vector<double> disVector;
		std::pair<bool, double> reV;
		bool futileBit1, futileBit2, futileBit;
		double  average, transX_residual, transY_residual, transA_residual, covX, covY, beli = 0.95;
		int badNum, deletedNum = 0;
		bool conflict = 0, have_deleted;
		std::pair<int, int> ele_local;
		while(1)
		{
			badNum = badLoops.size() - 1;
			deletedNum = 0;
			int testTime, test_time_const;

			test_time_const = int(3);
			if(test_time_const > int(goodLoops.size()) -1){
				test_time_const = int(goodLoops.size()) -1;
			}			
			for( ; badNum>=0; badNum--)
			{
				disVector.clear();
				testTime = test_time_const;
				// testTime = testTime-1;
				// cout<<"goodLoops.size() "<<goodLoops.size()<<endl;
				// cout<<"goodLoops.size() -1:  "<<goodLoops.size()-1<<endl;
				int bad_cluster_ID = clusterizer.getClusterID(badLoops[badNum]);

				cout<<" "<<endl;
				// cout<<"goodLoops size "<<goodLoops.size()<<" "<<" badLoops size "<<badLoops.size()<<endl;
				cout<<"cluster "<<bad_cluster_ID <<" bad loop: "<<badNum<<"   "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<
					"  error: "<<chi2LinkErrors[badLoops[badNum]]<<endl;

				good_loop_node_distance_pair_vector.clear();

				if(goodLoops.size() >= 1)
				{
					for(int local_i = 0; local_i < goodLoops.size(); local_i++)
					{
						ele_local.first = local_i;
						ele_local.second = abs(badLoops[badNum].first - goodLoops[local_i].first) + 
							abs(badLoops[badNum].second - goodLoops[local_i].second);
						good_loop_node_distance_pair_vector.push_back(ele_local);
					}	
					std::sort (good_loop_node_distance_pair_vector.begin(), good_loop_node_distance_pair_vector.end(), cmp_reverse);	
				}					

				// beli = 1 - utils::chi2_p(3, chi2LinkErrors[badLoops[badNum]]);
				cout<<"beli for this loop is "<<beli<<endl;

				if(testTime == -1)
				{
					cout<<"no good dis to complementary check as goodloop.size is "<<goodLoops.size()<<endl;
					// no good loop so do self check
					reV =  clusterizer.check_single_loop_odo(badLoops[badNum], clusterizer.LP_Trans_Covar_Map, futileBit1, beli);
					if((futileBit1) != 1)
					{
						if((reV.first == 1))
						{
							cout<<" * the dis is "<<reV.second<<endl;
							disVector.push_back(reV.second);
						}		
						else
						{
							cout<<"the dis is "<<reV.second<<endl;
							disVector.push_back(reV.second);
						}					
					}
					else{
						cout<<"futile dis is "<<reV.second<<endl;
					}					
				}
				else
				{
					have_deleted = 0;
					bool done_single_check = 0;
					int cluster_id_to_deleted;
					for(; testTime >=0; testTime--){
						cout<<"testTime "<<testTime<<endl;
						//if the error is too big, delete it directly
						if(chi2LinkErrors[badLoops[badNum]] > upperLimitChi){
							cout<<"delete big error loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;

							if(clusterizer.getClusterByID(clusterizer.getClusterID(badLoops[badNum])).size() == 1)
							{
								clusters_set.erase(clusterizer.getClusterID(badLoops[badNum]));
							}
				
							clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
							currentCluster.erase(badLoops[badNum]);
							deletedNum = deletedNum + 1;	
							have_deleted = 1;								
								
							break;						

						}
						// if(bad_cluster_ID == clusterizer.getClusterID(goodLoops[testTime])
						// {

						// }

						// cout<<"testTime "<<testTime<<endl;
						serial_good = good_loop_node_distance_pair_vector[testTime].first;
						cout<<"god loop "<<goodLoops[serial_good].first<<" "<<goodLoops[serial_good].second<<
							" with node dis "<<good_loop_node_distance_pair_vector[testTime].second<<endl;

						if((good_loop_node_distance_pair_vector[testTime].second > abs(badLoops[badNum].first - badLoops[badNum].second)) and
							done_single_check == 0)
						{
							cout<<"do single loop check"<<endl;

							reV =  clusterizer.check_single_loop_odo(badLoops[badNum], clusterizer.LP_Trans_Covar_Map, futileBit1, beli);//, double& statis
							done_single_check = 1;
							if((futileBit1) != 1)
							{
								if((reV.first == 1))
								{
									cout<<" * the dis is "<<reV.second<<endl;
									disVector.push_back(reV.second);
								}		
								else
								{
									cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
									if(clusterizer.getClusterByID(clusterizer.getClusterID(badLoops[badNum])).size() == 1)
									{
										clusters_set.erase(clusterizer.getClusterID(badLoops[badNum]));
										currentCluster.erase(badLoops[badNum]);	
										deletedNum = deletedNum + 1;	
										have_deleted =1;
										break;										
									}
									else
									{

										clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
										currentCluster.erase(badLoops[badNum]);									
										// deleted_loops.push_back(badLoops[badNum]);
										deletedNum = deletedNum + 1;
										cout<<"the dis is "<<reV.second<<endl;
										disVector.push_back(reV.second);
										have_deleted =1;
										break;
									}
							
								}					
							}
							else{
								cout<<"futile dis is "<<reV.second<<endl;
								disVector.push_back(reV.second);
							}
						}
						else
						{
							clusterizer.prepare_to_signle_loop_pair_check(badLoops[badNum], goodLoops[serial_good], 
								FullInfo,  clusterizer.LP_Trans_Covar_Map, length, futileBit1, futileBit2);

							if((futileBit2 or futileBit1) != 1){
								reV =  clusterizer.check_single_loop_inter_varying_belief(FullInfo, covX, covY, clusterizer.displayCov, 
									transX_residual, transY_residual, transA_residual, length[4][0], futileBit, beli);
								if((futileBit == 0))
								{
									if((reV.first == 1))
									{
										cout<<" * the dis is "<<reV.second<<endl;
										disVector.push_back(reV.second);
									}		
									else{
										disVector.push_back(reV.second);
										cout<<"the dis is "<<reV.second<<endl;
										// cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
										// clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
										// currentCluster.erase(badLoops[badNum]);
										// deletedNum = deletedNum + 1;
										// break;								
									}					
								}
								else{
									cout<<"futile dis is "<<reV.second<<endl;
									disVector.push_back(reV.second);
								}
							}
						}

					}
				}					

				if(have_deleted == 0)
				{
					average = accumulate( disVector.begin(), disVector.end(), 0.0)/disVector.size();  
					cout<<"averange dis is "<<average<<" threshold is "<<utils::chi2_continuous(3, beli)<<endl;
					cout<<"nearest loop dis is "<<disVector.back()<<endl;
					if(average > utils::chi2_continuous(3, beli) or disVector.back() > utils::chi2_continuous(3, 0.99))
					{
						cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<" from cluster "<<
							clusterizer.getClusterID(badLoops[badNum])<<endl;
						if(clusterizer.getClusterByID(clusterizer.getClusterID(badLoops[badNum])).size() == 1)
						{
							clusters_set.erase(clusterizer.getClusterID(badLoops[badNum]));
							currentCluster.erase(badLoops[badNum]);
							deletedNum = deletedNum + 1;	
						}
						else
						{
							clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
							int be = currentCluster.size();
							currentCluster.erase(badLoops[badNum]);	
							int af = currentCluster.size();
							if(be == af)
							{
								cout<<"failed to remove bad loop"<<endl;
								cin.get();
							}
							// deleted_loops.push_back(badLoops[badNum]);
							deletedNum = deletedNum + 1;	
						}
			
					}  	
				}	
				cout<<"deletedNum is "<<deletedNum<<endl;	
			}

			if(deletedNum == 0 )
			{
				// currentCluster = clusterizer.getClusterByID(clusterID);
				double current_size = currentCluster.size();
				if(current_size == 0)
				{
					cout<<"all loops are deleted "<<endl;
					return 0;
				}

				//the third parameter should be set to negative things, it is used to skip some certain loop when debug
				clusterNodeSet.clear();
				findRange = currentCluster.begin();
				findRangeEnd = currentCluster.end();
				for(; findRange != findRangeEnd; findRange++)
				{
					clusterNodeSet.insert((*findRange).first);
					clusterNodeSet.insert((*findRange).second);
				}
				startRange = *clusterNodeSet.begin();
				endRange   = *clusterNodeSet.rbegin();

				chi2LinkErrors.clear();
				// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
				chi2LinkErrors = 
								// gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
								// 		  startRange, endRange, odoEdgeError);
								gWrapper->optimize_active_robustKernel_multi(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);	

				activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
				if(chi2LinkErrors[IntPair(-1,0)] < utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95) and
					chi2LinkErrors[IntPair(-2,0)] < utils::chi2_continuous(edgeDimension* chi2LinkErrors[IntPair(-2,-1)], 0.95))
				{
					// cout<<"cluster size > 1, and some menbers are deleted"<<endl;
					// std::cin.get();
					// if()
					cout<<"return __ 1 "<<endl;
					return 1;
				}	
				else{
					cout<<"all graph error: "<<chi2LinkErrors[IntPair(-1,0)]<<endl<<"cluster error: "<<chi2LinkErrors[IntPair(-2,0)]<<endl; 
					cout<<"return __ 0 "<<endl;
					return 0;
				}		
			}
			else if(deletedNum == originalTotalNum or (currentCluster.size() == 0))
			{
				cout<<"all loops are deleted "<<endl;
				return 0;
			}

			// currentCluster = clusterizer.getClusterByID(clusterID);
			if(currentCluster.size() == 0)
			{
				cout<<"all loops are deleted"<<endl;
				return 0;
			}

			//the third parameter should be set to negative things, it is used to skip some certain loop when debug
			goodLoop = 0;
			// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			clusterNodeSet.clear();
			findRange = currentCluster.begin();
			findRangeEnd = currentCluster.end();
			for(; findRange != findRangeEnd; findRange++)
			{
				clusterNodeSet.insert((*findRange).first);
				clusterNodeSet.insert((*findRange).second);
			}
			startRange = *clusterNodeSet.begin();
			endRange   = *clusterNodeSet.rbegin();

			chi2LinkErrors.clear();
			// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			chi2LinkErrors = 
							// gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
							// 		  startRange, endRange, odoEdgeError);
							gWrapper->optimize_active_robustKernel_multi(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);	

			activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
			activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
			eIt = chi2LinkErrors.begin();
			eEnd = chi2LinkErrors.end();
			biggestErr = 0;
			goodLoops.clear();
			badLoops.clear();
			for(; eIt!=eEnd; eIt++)
			{
				if(eIt->first.first < 0)
					continue;
				if(eIt->second > biggestErr)
				{
					biggestErr = eIt->second;
					biggestErrorLoop = eIt->first;
				}

				if(eIt->second < utils::chi2(edgeDimension))		
				{
					goodLoop=goodLoop+1;
					goodLoops.push_back(eIt->first);
				}
				else
				{
					badLoops.push_back(eIt->first);
				}
				cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
					" has error "<< eIt->second<<endl;
			}
			cout<<"bad loop size:"<<badLoops.size()<<" good size "<<goodLoop<<endl;
			cout<<"activeChi2Graph is "<<activeChi2Graph<<endl;
			if(badLoops.size() == 0)
			{
				if(activeChi2Graph > utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95))
				{
					std::cerr<<"deleted Cluster "<<" have "<<realOrigianlNum<<" links "<<std::endl;
					return 0;
				}
				// std::cin.get();
				cout<<"bad loops num is 0, return multi 1"<<endl;
				return 1;
			}
			if(goodLoops.size() == 0)
			{
				std::cerr<<"deleted Cluster "<<" have "<<realOrigianlNum<<" links "<<std::endl;
				return 0;
			}			
		}
	}

	bool complementary_intra_multi_non_robust(IntPairSet & currentCluster, std::vector<double> & chiStatisVector, std::vector<std::pair<int, int> > & doubtLoopNumberVector, 
		int & originalLoopNum, std::vector<std::pair<double, std::pair<int, int> > > & chiStstis4eachLoop, std::set<int> & clusters_set)
	{
		IntPairDoubleMap chi2LinkErrors;
		std::vector<std::pair<std::pair<int, int>, double> >  odoEdgeRelateLC_Error;

		std::vector<std::pair<int,int> > goodLoops, badLoops;
		int originalTotalNum = currentCluster.size();
		originalLoopNum = originalTotalNum;
		int realOrigianlNum = originalLoopNum;
		std::vector<std::pair<int, int>> good_loop_node_distance_pair_vector;

		//the third parameter should be set to negative things, it is used to skip some certain loop when debug
		int startRange, endRange, serial_good;
		std::set<int> clusterNodeSet;
		IntPairSet::iterator findRange = currentCluster.begin(), findRangeEnd = currentCluster.end();
		for(; findRange != findRangeEnd; findRange++)
		{
			clusterNodeSet.insert((*findRange).first);
			clusterNodeSet.insert((*findRange).second);
		}
		startRange = *clusterNodeSet.begin();
		endRange   = *clusterNodeSet.rbegin();
        cout<<"node range from "<<startRange<<" to "<<endRange<<endl;

        // if originalTotalNum = 1, 
        string ba = "sdf", aa = "aa";

        std::vector< g2o::HyperGraph::Edge* > newEdgeVector;
        if(originalTotalNum == 1)
        {
        	std::vector<std::pair<int, int>> artificialLoops;
        	int maxVertexID = gWrapper->optimizer->vertices().size();
        	int futileDis;
        	if(clusterizer.futilePairVector.size() != 0)
        		futileDis = abs(clusterizer.futilePairVector[0].first - clusterizer.futilePairVector[0].second);
        	else
        		futileDis = maxVertexID;
        	
        	// generateNewEdge(*(currentCluster.begin()), newEdgeVector, artificialLoops, futileDis, maxVertexID);
			chi2LinkErrors =  
				// gWrapper->optimize_oneMember(currentCluster, small_numIterations, newEdgeVector, artificialLoops);
				gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
        }
        else
        {
			chi2LinkErrors =  		
				gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
        }

		float activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
		int   activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
		cout<<"active edge count: "<<activeEdgeCount<<endl;
		//check if there is someone loop whose chi2 error is bigger than the signgle edge limit error
		double biggestErr = 0, smallcHI = 200;
		std::pair<int, int> biggestErrorLoop(-1, -1);

		cout<<" "<<endl;
		IntPairDoubleMap::iterator eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		double goodLoop = 0;
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;
			if(eIt->second > biggestErr)
			{
				biggestErr = eIt->second;
				biggestErrorLoop = eIt->first;
			}
			if(eIt->second < utils::chi2_continuous(edgeDimension, 0.5))		
			{
				goodLoop=goodLoop+1;
				goodLoops.push_back(eIt->first);
			}
			else
			{
				badLoops.push_back(eIt->first);
				// cout<<"add bad loop "<<eIt->first.first<<" "<<eIt->first.second<<endl;
			}
			// cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
			// 	" has error "<< eIt->second<<endl;
		}


		//handle the bad loopps, through check their transform distance to good ones	
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		std::array<std::array<double, 5>, 5>  length;
		std::vector<double> disVector;
		std::pair<bool, double> reV;
		bool futileBit1, futileBit2, futileBit;
		double  average, transX_residual, transY_residual, transA_residual, covX, covY, beli = 0.95;
		int badNum, deletedNum = 0;
		bool conflict = 0, have_deleted;
		std::pair<int, int> ele_local;
		while(1)
		{
			badNum = badLoops.size() - 1;
			deletedNum = 0;
			int testTime;
			for( ; badNum>=0; badNum--)
			{
				testTime = 9;
				disVector.clear();
				if(testTime > int(goodLoops.size()) -1){
					testTime = int(goodLoops.size()) -1;
				}
				// testTime = testTime-1;
				// cout<<"goodLoops.size() "<<goodLoops.size()<<endl;
				// cout<<"goodLoops.size() -1:  "<<goodLoops.size()-1<<endl;
				int bad_cluster_ID = clusterizer.getClusterID(badLoops[badNum]);

				cout<<" "<<endl;
				// cout<<"goodLoops size "<<goodLoops.size()<<" "<<" badLoops size "<<badLoops.size()<<endl;
				cout<<"bad loop "<<badNum<<"   "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<
					" has error "<<chi2LinkErrors[badLoops[badNum]]<<endl;

				good_loop_node_distance_pair_vector.clear();
				

				if(goodLoops.size() >= 1)
				{
					for(int local_i = 0; local_i < goodLoops.size(); local_i++)
					{
						ele_local.first = local_i;
						ele_local.second = abs(badLoops[badNum].first - goodLoops[local_i].first) + 
							abs(badLoops[badNum].second - goodLoops[local_i].second);
						good_loop_node_distance_pair_vector.push_back(ele_local);
					}	
					std::sort (good_loop_node_distance_pair_vector.begin(), good_loop_node_distance_pair_vector.end(), cmp_reverse);	
				}					

				// beli = 1 - utils::chi2_p(3, chi2LinkErrors[badLoops[badNum]]);
				cout<<"beli for this loop is "<<beli<<endl;

				if(testTime == -1)
				{
					// no good loop so do self check
					reV =  clusterizer.check_single_loop_odo(badLoops[badNum], clusterizer.LP_Trans_Covar_Map, futileBit1, beli);
					if((futileBit1) != 1)
					{
						if((reV.first == 1))
						{
							cout<<" * the dis is "<<reV.second<<endl;
							disVector.push_back(reV.second);
						}		
						else
						{
							cout<<"the dis is "<<reV.second<<endl;
							disVector.push_back(reV.second);
					
						}					
					}
					else{
						cout<<"futile dis is "<<reV.second<<endl;
					}					
				}
				else
				{
					have_deleted = 0;
					bool done_single_check = 0;
					int cluster_id_to_deleted;
					for(; testTime >=0; testTime--){
						// cout<<"testTime "<<testTime<<endl;
						//if the error is too big, delete it directly
						if(chi2LinkErrors[badLoops[badNum]] > upperLimitChi){
							cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;

							if(clusterizer.getClusterByID(clusterizer.getClusterID(badLoops[badNum])).size() == 1)
							{
								clusters_set.erase(clusterizer.getClusterID(badLoops[badNum]));
							}
							else
							{
								clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
								currentCluster.erase(badLoops[badNum]);
								deletedNum = deletedNum + 1;	
								have_deleted = 1;
								break;
							}							

						}
						// if(bad_cluster_ID == clusterizer.getClusterID(goodLoops[testTime])
						// {

						// }

						// cout<<"testTime "<<testTime<<endl;
						serial_good = good_loop_node_distance_pair_vector[testTime].first;
						cout<<"god loop "<<goodLoops[serial_good].first<<" "<<goodLoops[serial_good].second<<
							" with node dis "<<good_loop_node_distance_pair_vector[testTime].second<<endl;

						if((good_loop_node_distance_pair_vector[testTime].second > abs(badLoops[badNum].first - badLoops[badNum].second)) and
							done_single_check == 0)
						{
							cout<<"do single loop check"<<endl;

							reV =  clusterizer.check_single_loop_odo(badLoops[badNum], clusterizer.LP_Trans_Covar_Map, futileBit1, beli);//, double& statis
							done_single_check = 1;
							if((futileBit1) != 1)
							{
								if((reV.first == 1))
								{
									cout<<" * the dis is "<<reV.second<<endl;
									disVector.push_back(reV.second);
								}		
								else
								{
									cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
									if(clusterizer.getClusterByID(clusterizer.getClusterID(badLoops[badNum])).size() == 1)
									{
										clusters_set.erase(clusterizer.getClusterID(badLoops[badNum]));
									}
									else
									{

										clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
										currentCluster.erase(badLoops[badNum]);									
										// deleted_loops.push_back(badLoops[badNum]);
										deletedNum = deletedNum + 1;
										cout<<"the dis is "<<reV.second<<endl;
										disVector.push_back(reV.second);
										have_deleted =1;
										break;
									}
							
								}					
							}
							else{
								cout<<"futile dis is "<<reV.second<<endl;
							}
						}
						else
						{
							clusterizer.prepare_to_signle_loop_pair_check(badLoops[badNum], goodLoops[serial_good], 
								FullInfo,  clusterizer.LP_Trans_Covar_Map, length, futileBit1, futileBit2);

							if((futileBit2 or futileBit1) != 1){
								reV =  clusterizer.check_single_loop_inter_varying_belief(FullInfo, covX, covY, clusterizer.displayCov, 
									transX_residual, transY_residual, transA_residual, length[4][0], futileBit, beli);
								if((futileBit == 0))
								{
									if((reV.first == 1))
									{
										cout<<" * the dis is "<<reV.second<<endl;
										disVector.push_back(reV.second);
									}		
									else{
										disVector.push_back(reV.second);
										cout<<"the dis is "<<reV.second<<endl;
										// cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
										// clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
										// currentCluster.erase(badLoops[badNum]);
										// deletedNum = deletedNum + 1;
										// break;								
									}					
								}
								else{
									cout<<"futile dis is "<<reV.second<<endl;
								}
							}
						}

					}
				}					

				if(have_deleted == 0)
				{
					average = accumulate( disVector.begin(), disVector.end(), 0.0)/disVector.size();  
					cout<<"averange dis is "<<average<<" threshold is "<<utils::chi2_continuous(3, beli)<<endl;
					if(average > utils::chi2_continuous(3, beli))
					{
						cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
						if(clusterizer.getClusterByID(clusterizer.getClusterID(badLoops[badNum])).size() == 1)
						{
							clusters_set.erase(clusterizer.getClusterID(badLoops[badNum]));
						}
						else
						{
							clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
							int be = currentCluster.size();
							currentCluster.erase(badLoops[badNum]);	
							int af = currentCluster.size();
							if(be == af)
							{
								cout<<"failed to remove bad loop"<<endl;
								cin.get();
							}
							// deleted_loops.push_back(badLoops[badNum]);
							deletedNum = deletedNum + 1;	
						}
			
					}  	
				}		
			}

			if(deletedNum == 0 )
			{
				// currentCluster = clusterizer.getClusterByID(clusterID);
				double current_size = currentCluster.size();
				if(current_size == 0)
				{
					cout<<"all loops are deleted "<<endl;
					return 0;
				}

				//the third parameter should be set to negative things, it is used to skip some certain loop when debug
				clusterNodeSet.clear();
				findRange = currentCluster.begin();
				findRangeEnd = currentCluster.end();
				for(; findRange != findRangeEnd; findRange++)
				{
					clusterNodeSet.insert((*findRange).first);
					clusterNodeSet.insert((*findRange).second);
				}
				startRange = *clusterNodeSet.begin();
				endRange   = *clusterNodeSet.rbegin();

				chi2LinkErrors.clear();
				// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
				chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);

				activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
				if(chi2LinkErrors[IntPair(-1,0)] < utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95) and
					chi2LinkErrors[IntPair(-2,0)] < utils::chi2_continuous(edgeDimension* chi2LinkErrors[IntPair(-2,-1)], 0.95))
				{
					// cout<<"cluster size > 1, and some menbers are deleted"<<endl;
					// std::cin.get();
					// if()
					return 1;
				}	
				else{
					return 0;
				}		
			}
			else if(deletedNum == originalTotalNum or (currentCluster.size() == 0))
			{
				cout<<"all loops are deleted "<<endl;
				return 0;
			}

			// currentCluster = clusterizer.getClusterByID(clusterID);
			if(currentCluster.size() == 0)
			{
				cout<<"all loops are deleted"<<endl;
				return 0;
			}

			//the third parameter should be set to negative things, it is used to skip some certain loop when debug
			goodLoop = 0;
			// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			clusterNodeSet.clear();
			findRange = currentCluster.begin();
			findRangeEnd = currentCluster.end();
			for(; findRange != findRangeEnd; findRange++)
			{
				clusterNodeSet.insert((*findRange).first);
				clusterNodeSet.insert((*findRange).second);
			}
			startRange = *clusterNodeSet.begin();
			endRange   = *clusterNodeSet.rbegin();

			chi2LinkErrors.clear();
			// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
			activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
			eIt = chi2LinkErrors.begin();
			eEnd = chi2LinkErrors.end();
			biggestErr = 0;
			goodLoops.clear();
			badLoops.clear();
			for(; eIt!=eEnd; eIt++)
			{
				if(eIt->first.first < 0)
					continue;
				if(eIt->second > biggestErr)
				{
					biggestErr = eIt->second;
					biggestErrorLoop = eIt->first;
				}

				if(eIt->second < utils::chi2(edgeDimension))		
				{
					goodLoop=goodLoop+1;
					goodLoops.push_back(eIt->first);
				}
				else
				{
					badLoops.push_back(eIt->first);
				}
				cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
					" has error "<< eIt->second<<endl;
			}
			cout<<"bad loop size:"<<badLoops.size()<<" good size "<<goodLoop<<endl;
			cout<<"activeChi2Graph is "<<activeChi2Graph<<endl;
			if(badLoops.size() == 0)
			{
				if(activeChi2Graph > utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95))
				{
					std::cerr<<"deleted Cluster "<<" have "<<realOrigianlNum<<" links "<<std::endl;
					return 0;
				}
				// std::cin.get();
				return 1;
			}
			if(goodLoops.size() == 0)
			{
				std::cerr<<"deleted Cluster "<<" have "<<realOrigianlNum<<" links "<<std::endl;
				return 0;
			}			
		}
	}
	// the same idea with complementary intra third,
	bool complementary_intra_multi_third(IntPairSet & currentCluster, std::vector<double> & chiStatisVector, 
		std::vector<std::pair<int, int> > & doubtLoopNumberVector, 
		int & originalLoopNum, std::vector<std::pair<double, std::pair<int, int> > > & chiStstis4eachLoop, std::set<int> & clusters_set)
	{
		IntPairDoublePairMap chi2LinkErrors;

		std::vector<std::pair<int,int> > goodLoops, badLoops;
		int originalTotalNum = currentCluster.size();
		originalLoopNum = originalTotalNum;
		int realOrigianlNum = originalLoopNum;
		std::vector<std::pair<int, int>> good_loop_node_distance_pair_vector;

		//the third parameter should be set to negative things, it is used to skip some certain loop when debug
		int startRange, endRange, serial_good;
		std::set<int> clusterNodeSet;
		IntPairSet::iterator findRange = currentCluster.begin(), findRangeEnd = currentCluster.end();
		for(; findRange != findRangeEnd; findRange++)
		{
			clusterNodeSet.insert((*findRange).first);
			clusterNodeSet.insert((*findRange).second);
		}
		startRange = *clusterNodeSet.begin();
		endRange   = *clusterNodeSet.rbegin();
        cout<<"node range from "<<startRange<<" to "<<endRange<<endl;

        // if originalTotalNum = 1, 
        string ba = "sdf", aa = "aa";

        std::vector< g2o::HyperGraph::Edge* > newEdgeVector;
        if(originalTotalNum == 1)
        {
        	std::vector<std::pair<int, int>> artificialLoops;
        	int maxVertexID = gWrapper->optimizer->vertices().size();
        	int futileDis;
        	if(clusterizer.futilePairVector.size() != 0)
        		futileDis = abs(clusterizer.futilePairVector[0].first - clusterizer.futilePairVector[0].second);
        	else
        		futileDis = maxVertexID;
        	
        	// generateNewEdge(*(currentCluster.begin()), newEdgeVector, artificialLoops, futileDis, maxVertexID);
			chi2LinkErrors =  
				// gWrapper->optimize_oneMember(currentCluster, small_numIterations, newEdgeVector, artificialLoops);
				gWrapper->optimize_active_robustKernel_pair_return(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);
        }
        else
        {
			chi2LinkErrors =  		
				gWrapper->optimize_active_robustKernel_pair_return(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);
        }

		float activeChi2Graph = chi2LinkErrors[IntPair(-1,0)].first;
		int   activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)].first;
		cout<<"active edge count: "<<activeEdgeCount<<endl;
		//check if there is someone loop whose chi2 error is bigger than the signgle edge limit error
		double smallcHI = 200;
		std::pair<int, int> biggestErrorLoop(-1, -1);

		cout<<" "<<endl;
		IntPairDoublePairMap::iterator eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		double goodLoop = 0;
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;

			if(eIt->second.second < utils::chi2_continuous(edgeDimension, 0.95))		
			{
				goodLoop=goodLoop+1;
				goodLoops.push_back(eIt->first);
			}
			else
			{
				badLoops.push_back(eIt->first);
				// cout<<"add bad loop "<<eIt->first.first<<" "<<eIt->first.second<<endl;
			}
			// cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
			// 	" has error "<< eIt->second<<endl;
		}

		// if(goodLoops.size() == 0)
		// {
		// 	return 0;
		// }
		//handle the bad loopps, through check their transform distance to good ones	
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		std::array<std::array<double, 5>, 5>  length;
		std::vector<double> disVector;
		std::pair<bool, double> reV;
		bool futileBit1, futileBit2, futileBit;
		double  average, transX_residual, transY_residual, transA_residual, covX, covY, beli = 0.95;
		int badNum, deletedNum = 0;
		bool conflict = 0, have_deleted;
		std::pair<int, int> ele_local;
		while(1)
		{
			badNum = badLoops.size() - 1;
			deletedNum = 0;
			int testTime;
			for( ; badNum>=0; badNum--)
			{
				testTime = 9;
				disVector.clear();
				if(testTime > int(goodLoops.size()) -1){
					testTime = int(goodLoops.size()) -1;
				}
				// testTime = testTime-1;
				// cout<<"goodLoops.size() "<<goodLoops.size()<<endl;
				// cout<<"goodLoops.size() -1:  "<<goodLoops.size()-1<<endl;
				int bad_cluster_ID = clusterizer.getClusterID(badLoops[badNum]);

				cout<<" "<<endl;
				// cout<<"goodLoops size "<<goodLoops.size()<<" "<<" badLoops size "<<badLoops.size()<<endl;
				cout<<"bad loop "<<badNum<<"   "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<
					" has robust error "<<chi2LinkErrors[badLoops[badNum]].first
					<<" and non-robust err "<<chi2LinkErrors[badLoops[badNum]].second <<endl;

				good_loop_node_distance_pair_vector.clear();
				

				if(goodLoops.size() >= 1)
				{
					for(int local_i = 0; local_i < goodLoops.size(); local_i++)
					{
						ele_local.first = local_i;
						ele_local.second = abs(badLoops[badNum].first - goodLoops[local_i].first) + 
							abs(badLoops[badNum].second - goodLoops[local_i].second);
						good_loop_node_distance_pair_vector.push_back(ele_local);
					}	
					std::sort (good_loop_node_distance_pair_vector.begin(), good_loop_node_distance_pair_vector.end(), cmp_reverse);	
				}					

				// beli = 1 - utils::chi2_p(3, chi2LinkErrors[badLoops[badNum]]);
				// cout<<"beli for this loop is "<<beli<<endl;

				if(testTime == -1)
				{
					// no good loop so do self check
					reV =  clusterizer.check_single_loop_odo(badLoops[badNum], clusterizer.LP_Trans_Covar_Map, futileBit1, beli);
					if((futileBit1) != 1)
					{
						if((reV.first == 1))
						{
							cout<<" * the dis is "<<reV.second<<endl;
							disVector.push_back(reV.second);
						}		
						else
						{
							cout<<"the dis is "<<reV.second<<endl;
							disVector.push_back(reV.second);
					
						}					
					}
					else{
						cout<<"futile dis is "<<reV.second<<endl;
					}					
				}
				else
				{
					have_deleted = 0;
					bool done_single_check = 0;
					int cluster_id_to_deleted;
					int reverse_test_time = testTime;
					for(; testTime >=0; testTime--){
						// cout<<"testTime "<<testTime<<endl;
						//if the error is too big, delete it directly
						if(chi2LinkErrors[badLoops[badNum]].first > utils::chi2_continuous(3, 0.95)){
							cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;

							if(clusterizer.getClusterByID(clusterizer.getClusterID(badLoops[badNum])).size() == 1)
							{
								clusters_set.erase(clusterizer.getClusterID(badLoops[badNum]));
							}
							else
							{
								clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
								currentCluster.erase(badLoops[badNum]);
								deletedNum = deletedNum + 1;	
								have_deleted = 1;
								break;
							}							

						}
						// if(bad_cluster_ID == clusterizer.getClusterID(goodLoops[testTime])
						// {

						// }

						// cout<<"testTime "<<testTime<<endl;
						serial_good = good_loop_node_distance_pair_vector[reverse_test_time - testTime].first;
						cout<<"god loop "<<goodLoops[serial_good].first<<" "<<goodLoops[serial_good].second<<
							" with node dis "<<good_loop_node_distance_pair_vector[testTime].second<<endl;

						if((good_loop_node_distance_pair_vector[testTime].second > abs(badLoops[badNum].first - badLoops[badNum].second)) and
							done_single_check == 0)
						{
							cout<<"do single loop check"<<endl;

							reV =  clusterizer.check_single_loop_odo(badLoops[badNum], clusterizer.LP_Trans_Covar_Map, futileBit1, beli);//, double& statis
							done_single_check = 1;
							if((futileBit1) != 1)
							{
								if((reV.first == 1))
								{
									cout<<" * the dis is "<<reV.second<<endl;
									disVector.push_back(reV.second);
								}		
								else
								{
									cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
									if(clusterizer.getClusterByID(clusterizer.getClusterID(badLoops[badNum])).size() == 1)
									{
										clusters_set.erase(clusterizer.getClusterID(badLoops[badNum]));
									}
									else
									{

										clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
										currentCluster.erase(badLoops[badNum]);									
										// deleted_loops.push_back(badLoops[badNum]);
										deletedNum = deletedNum + 1;
										cout<<"the dis is "<<reV.second<<endl;
										disVector.push_back(reV.second);
										have_deleted =1;
										break;
									}
							
								}					
							}
							else{
								cout<<"futile dis is "<<reV.second<<endl;
							}
						}
						else
						{
							clusterizer.prepare_to_signle_loop_pair_check(badLoops[badNum], goodLoops[serial_good], 
								FullInfo,  clusterizer.LP_Trans_Covar_Map, length, futileBit1, futileBit2);

							if((futileBit2 or futileBit1) != 1){
								reV =  clusterizer.check_single_loop_inter_varying_belief(FullInfo, covX, covY, clusterizer.displayCov, 
									transX_residual, transY_residual, transA_residual, length[4][0], futileBit, beli);
								if((futileBit == 0))
								{
									if((reV.first == 1))
									{
										cout<<" * the dis is "<<reV.second<<endl;
										disVector.push_back(reV.second);
									}		
									else{
										disVector.push_back(reV.second);
										cout<<"the dis is "<<reV.second<<endl;
										// cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
										// clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
										// currentCluster.erase(badLoops[badNum]);
										// deletedNum = deletedNum + 1;
										// break;								
									}					
								}
								else{
									cout<<"futile dis is "<<reV.second<<endl;
								}
							}
						}

					}
				}					

				if(have_deleted == 0)
				{
					if(disVector.size() == 0)
					{
						cout<<"this is intra multi function, no valid nearby loop to check this double loop"<<endl;
						if(chi2LinkErrors[badLoops[badNum]].first > utils::chi2_continuous(3, 0.5))
						{
							cout<<"do you want to delete this loop?"<<endl<<"it has no valid neighbour but its error is larger than %50"<<endl;
						}
					}
					else
					{
						average = accumulate( disVector.begin(), disVector.end(), 0.0)/disVector.size();  
						cout<<"averange dis is "<<average<<" threshold is "<<utils::chi2_continuous(3, beli)<<endl;
						if(average > utils::chi2_continuous(3, beli))
						{
							cout<<"delete loop "<<badLoops[badNum].first<<" "<<badLoops[badNum].second<<endl;
							if(clusterizer.getClusterByID(clusterizer.getClusterID(badLoops[badNum])).size() == 1)
							{
								clusters_set.erase(clusterizer.getClusterID(badLoops[badNum]));
							}
							else
							{
								clusterizer.setClusterID(badLoops[badNum],ID_IGNORE); 
								int be = currentCluster.size();
								currentCluster.erase(badLoops[badNum]);	
								int af = currentCluster.size();
								if(be == af)
								{
									cout<<"failed to remove bad loop"<<endl;
									cin.get();
								}
								// deleted_loops.push_back(badLoops[badNum]);
								deletedNum = deletedNum + 1;	
							}
						} 
					}
				}		
			}

			if(deletedNum == 0 )
			{
				// currentCluster = clusterizer.getClusterByID(clusterID);
				double current_size = currentCluster.size();
				if(current_size == 0)
				{
					cout<<"all loops are deleted "<<endl;
					return 0;
				}

				//the third parameter should be set to negative things, it is used to skip some certain loop when debug
				clusterNodeSet.clear();
				findRange = currentCluster.begin();
				findRangeEnd = currentCluster.end();
				for(; findRange != findRangeEnd; findRange++)
				{
					clusterNodeSet.insert((*findRange).first);
					clusterNodeSet.insert((*findRange).second);
				}
				startRange = *clusterNodeSet.begin();
				endRange   = *clusterNodeSet.rbegin();

				chi2LinkErrors.clear();
				// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
				chi2LinkErrors = gWrapper->optimize_active_robustKernel_pair_return(currentCluster, small_numIterations, 
										  startRange, endRange, odoEdgeError);

				activeEdgeCount = int(chi2LinkErrors[IntPair(-1,-1)].first);
				if(chi2LinkErrors[IntPair(-1,0)].first < utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95) and
					chi2LinkErrors[IntPair(-2,0)].first < utils::chi2_continuous(edgeDimension* (int(chi2LinkErrors[IntPair(-2,-1)].first)), 0.95))
				{
					// cout<<"cluster size > 1, and some menbers are deleted"<<endl;
					// std::cin.get();
					// if()
					return 1;
				}	
				else{
					return 0;
				}		
			}
			else if(deletedNum == originalTotalNum or (currentCluster.size() == 0))
			{
				cout<<"all loops are deleted "<<endl;
				return 0;
			}

			// currentCluster = clusterizer.getClusterByID(clusterID);
			if(currentCluster.size() == 0)
			{
				cout<<"all loops are deleted"<<endl;
				return 0;
			}

			//the third parameter should be set to negative things, it is used to skip some certain loop when debug
			goodLoop = 0;
			// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			clusterNodeSet.clear();
			findRange = currentCluster.begin();
			findRangeEnd = currentCluster.end();
			for(; findRange != findRangeEnd; findRange++)
			{
				clusterNodeSet.insert((*findRange).first);
				clusterNodeSet.insert((*findRange).second);
			}
			startRange = *clusterNodeSet.begin();
			endRange   = *clusterNodeSet.rbegin();

			chi2LinkErrors.clear();
			// chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			chi2LinkErrors = gWrapper->optimize_active_robustKernel_pair_return(currentCluster, small_numIterations, 
									  startRange, endRange, odoEdgeError);

			activeChi2Graph = chi2LinkErrors[IntPair(-1,0)].first;
			activeEdgeCount = int(chi2LinkErrors[IntPair(-1,-1)].first);
			eIt = chi2LinkErrors.begin();
			eEnd = chi2LinkErrors.end();
			goodLoops.clear();
			badLoops.clear();
			for(; eIt!=eEnd; eIt++)
			{
				if(eIt->first.first < 0)
					continue;

				if(eIt->second.second < utils::chi2(edgeDimension))		
				{
					goodLoop=goodLoop+1;
					goodLoops.push_back(eIt->first);
				}
				else
				{
					badLoops.push_back(eIt->first);
				}
				cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" in cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<
					" robust_ERR: "<< eIt->second.first<<" non-robust: "<< eIt->second.second<<endl;
			}
			cout<<"bad loop size:"<<badLoops.size()<<" good size "<<goodLoop<<endl;
			cout<<"activeChi2Graph is "<<activeChi2Graph<<endl;
			if(badLoops.size() == 0)
			{
				if(activeChi2Graph > utils::chi2_continuous(edgeDimension* activeEdgeCount, 0.95))
				{
					std::cerr<<"deleted Cluster "<<" have "<<realOrigianlNum<<" links "<<std::endl;
					return 0;
				}
				// std::cin.get();
				return 1;
			}
			if(goodLoops.size() == 0)
			{
				std::cerr<<"deleted Cluster "<<" have "<<realOrigianlNum<<" links "<<std::endl;
				return 0;
			}			
		}
	}	
	bool gatherLinks(const IntSet& clusterList, IntPairSet& links)
	{
		IntSet::const_iterator it = clusterList.begin(), end = clusterList.end();
		// cout<<"gather links form cluster: ";
		for( ; it!=end; it++)
		{
			// cout<<(*it)<<" ";
			IntPairSet& currentCluster = clusterizer.getClusterByID(*it);
			links.insert(currentCluster.begin(),currentCluster.end());
		}
		// cout<<endl;
		return true;
	}
	bool gatherLinks(const std::vector<int> & clusterList, IntPairSet& links)
	{
		// cout<<"gather links form cluster: ";
		for(int gather = 0 ; gather < clusterList.size(); gather++)
		{
			// cout<<(clusterList[gather])<<" ";
			IntPairSet& currentCluster = clusterizer.getClusterByID(clusterList[gather]);
			links.insert(currentCluster.begin(),currentCluster.end());
		}
		// cout<<endl;
		return true;
	}
	bool gatherLinks_equalElement(const IntSet& clusterList, IntPairSet& links)
	{
		if(clusterList.size() == 0)
		{
			cout<<"the cluster set is empty"<<endl;
			return 0;
		}

		IntSet::const_iterator it = clusterList.begin(), end = clusterList.end();
		int clusterID = *it;
		// int minElementNum = (clusterizer.getClusterByID(clusterID)).size();
		int minElementNum = clusterizer.getClusterSize(clusterID);

		it++;

		for( ; it!=end; it++)
		{
			if(clusterizer.getClusterByID(*it).size() < minElementNum)
				minElementNum= clusterizer.getClusterByID(*it).size();
		}

		cout<<"the biggest share ele number is "<<minElementNum<<endl;
		it = clusterList.begin(), end = clusterList.end();
		for( ; it!=end; it++)
		{
			IntPairSet& currentCluster = clusterizer.getClusterByID(*it);
			IntPairSet::const_iterator  pointer = currentCluster.begin();
			for(int toINsert = 0; toINsert < minElementNum; toINsert++)
			{
				if(pointer == currentCluster.end())
				{
					printf("This error about start node and end node is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
				links.insert(*(pointer));
				pointer++;
			}
		}
		return true;
	}



	int getNumberIn(const IntSet& clusterList)
	{

		IntPairSet links;
		gatherLinks(clusterList,  links);
		{
			IntSet::const_iterator it = clusterList.begin(), end = clusterList.end();
			// cout<<"gather links form cluster: ";
			for( ; it!=end; it++)
			{
				// cout<<(*it)<<" ";
				IntPairSet& currentCluster = clusterizer.getClusterByID(*it);
				links.insert(currentCluster.begin(),currentCluster.end());
			}
			// cout<<endl;
			return links.size();
		}
	}	

	bool interClusterConsistent(IntSet& H, IntSet& goodSet, IntSet& RejectSet)
	{

		if(H.empty()) return true;

		IntPairSet activeLinks;

		gatherLinks(goodSet,activeLinks);
		gatherLinks(H,activeLinks);

		IntPairDoubleMap
			// linkErrors = gWrapper->optimize(activeLinks,small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			linkErrors = gWrapper->optimize_robustKernel(activeLinks,small_numIterations, odoEdgeError);
   //      	int maxVertexID = gWrapper->optimizer->vertices().size();
   //      	int futileDis;
   //      	if(clusterizer.futilePairVector.size() != 0)
   //      		futileDis = abs(clusterizer.futilePairVector[0].first - clusterizer.futilePairVector[0].second);
   //      	else
   //      		futileDis = maxVertexID;
        	
   //      	// generateNewEdge(*(currentCluster.begin()), newEdgeVector, artificialLoops, futileDis, maxVertexID);
			// chi2LinkErrors =  
			// 	// gWrapper->optimize_oneMember(currentCluster, small_numIterations, newEdgeVector, artificialLoops);
			// 	gWrapper->optimize_active_robustKernel(currentCluster, small_numIterations, 
			// 							  startRange, endRange, odoEdgeError);


		float activeChi2Graph = linkErrors[IntPair(-1,0)];
		int   activeEdgeCount = linkErrors[IntPair(-1,-1)];

		double allLinksError = 0;
		for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >= 0) // Special case for the two values we returned with a -1 as the first
			{
				cout<<"loop "<<it->first.first<<" "<<it->first.second<<
					" in cluster "<<clusterizer.getClusterID(it->first)<<" has error "<<it->second<<endl;
				allLinksError += it->second;
			}

		}


		IntPairSet::iterator it = activeLinks.begin(), end = activeLinks.end();

		std::map< int, double > errorMap;
		bool all_pass = 1;

		for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >=0)
			{
				int thisLinkClusterID = clusterizer.getClusterID(it->first);
				errorMap[thisLinkClusterID]	+= it->second;
			}
		}

		IntSet::iterator sIt = H.begin(), sEnd = H.end();

		double min_CI = 0;
		int rejectID = -1;

		for ( ; sIt!=sEnd; sIt++)
		{
			cout<<"all pass = "<<all_pass<<endl;
			// double CI = errorMap[*sIt]/utils::chi2(edgeDimension*clusterizer.getClusterByID(*sIt).size());
			if(errorMap[*sIt] > utils::chi2(edgeDimension*clusterizer.getClusterByID(*sIt).size()))
			{
				all_pass = 0;
				cout<<"errorMap[*sIt] = "<<errorMap[*sIt]<<"   "<<"cluster.size = "<<clusterizer.getClusterByID(*sIt).size()<<endl;
				// break;
			}
			double CI = errorMap[*sIt]/utils::chi2(edgeDimension*clusterizer.getClusterByID(*sIt).size());

			if( CI >= min_CI ) // Just looking for the ones in H
			{
				min_CI = CI;
				rejectID = *sIt;
			}
			cout<<"cluster "<<*sIt<<" has relative error "<<CI<<endl;
		}

		if(all_pass == 1)
			cout<<"exit for, all pass = "<<all_pass<<endl;
		else
			cout<<"***lower exit for, all pass = "<<all_pass<<endl;
	
		cout<<activeChi2Graph<<"  "<<utils::chi2(edgeDimension*activeEdgeCount)<<endl;
		cout<<allLinksError<<"  "<<utils::chi2(edgeDimension*(activeLinks.size()))<<endl;
		// if(
		// 		activeChi2Graph 	< 	utils::chi2(edgeDimension*activeEdgeCount)
		// 		and
		// 		allLinksError 		< 	utils::chi2(edgeDimension*(activeLinks.size()))
		// 	)
		if(
				activeChi2Graph 	< 	utils::chi2(edgeDimension*activeEdgeCount)
				and
				allLinksError 		< 	utils::chi2(edgeDimension*(activeLinks.size()))
				// and
				// all_pass == 1
			)			
		{
			cout<<"whole graoh chi2 error: "<<activeChi2Graph<<" all links error: "<<allLinksError<<endl;
			goodSet.insert(H.begin(),H.end());
			return true; // all done .. we found a consistent solution
		}
		else
		{
			// Find which cluster is causing the problems
			// Iterate over everything
			// Find the error for each cluster
			// Sum cluster-wise

			// IntPairSet::iterator it = activeLinks.begin(), end = activeLinks.end();

			// std::map< int, double > errorMap;

			// for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
			// 		it!=end;
			// 		it++)
			// {
			// 	if(it->first.first >=0)
			// 	{
			// 		int thisLinkClusterID = clusterizer.getClusterID(it->first);
			// 		errorMap[thisLinkClusterID]	+= it->second;
			// 	}
			// }

			// IntSet::iterator sIt = H.begin(), sEnd = H.end();

			// double min_CI = 0;
			// int rejectID = -1;

			// for ( ; sIt!=sEnd; sIt++)
			// {
			// 	double CI = errorMap[*sIt]/utils::chi2(edgeDimension*clusterizer.getClusterByID(*sIt).size());

			// 	if( CI >= min_CI ) // Just looking for the ones in H
			// 	{
			// 		min_CI = CI;
			// 		rejectID = *sIt;
			// 	}
			// 	cout<<"cluster "<<*sIt<<" has relative error "<<CI<<endl;
			// }
			H.erase(rejectID);
			cout<<"delete cluster "<<rejectID<<endl;
			RejectSet.insert(rejectID);

			return interClusterConsistent(H,goodSet,RejectSet);
		}
		return true;
	}


	void findGoodLoopRegin(std::vector<bool> &BitVector, std::vector< std::pair<int,int> > &goodLoopRegin)
	{
		int  start = -1, middle, times=0;
		std::pair<int, int>  midd;
		for(int iterator = 0; iterator < BitVector.size(); iterator++)
		{
			// if(VectorLoopError[iterator].second > utils::chi2(3, 0.95))
			// cout<<iterator<<" "<<BitVector[iterator] <<endl;
			if(BitVector[iterator] == 1 and iterator!= BitVector.size()-1)
			{
				if(times == 0)
				{
					start = iterator;
					times = 1;
				}
				else
				{
					times =times +1;
				}
			}
			else
			{
				if(times > 0)
				{
					midd.first = start;
					midd.second = times;
//					if(times > 0)
                    goodLoopRegin.push_back(midd);
					
//					 cout<<"start = "<<goodLoopRegin.back().first<<" times = "<<goodLoopRegin.back().second<<endl;
					times = 0;
				}
				else if(times == 0 and BitVector[iterator] == 1)
                {
                    midd.first = iterator;
                    midd.second = 1;
                    goodLoopRegin.push_back(midd);
                }
			}
		}
	}

	void handleBigErrorInInter(IntSet& H, IntSet& goodSet, IntSet& RejectSet, double beliefResidualError, double beliefTransError,
                               double goodLoopBelief, IntPairDoubleMap & linkErrors,  bool  & deleteSome, IntPairSet & activeLinks)
    {
        if(H.empty()) return ;
        cout<<" "<<endl;
        cout<<"three belief:"<<beliefResidualError<<" "<<beliefTransError<<" "<<goodLoopBelief<<endl;
		activeLinks.clear();
        std::vector<int> multiClusterExist;
        std::map<int, std::vector<std::pair<int, int> > >  mapCluserIDLoops;
        std::map<int, std::vector< double> >  mapCluserIDLoopErrors;
        std::vector<std::pair<std::pair<int,int> , double> > VectorLoopError;
        std::pair<std::pair<int,int> , double> pairLoopErr;
        std::set<int> good_clusters;

        good_clusters.insert(H.begin(),H.end());
        if(good_clusters.size() != H.size())
        {
        	cout<<"two set elements not equal"<<endl;
        	printf("This error about start node and end node is in %s on line %d\n",  __FILE__, __LINE__);
        	exit(0);
        }
        gatherLinks(goodSet,activeLinks);
        gatherLinks(H,activeLinks);


        // linkErrors = gWrapper->optimize(activeLinks, overOptimizeIteration, odoEdgeRelateLC_Error, odoEdgeError);
        linkErrors = gWrapper-> optimize_robustKernel(activeLinks, overOptimizeIteration,odoEdgeError)	;

        float activeChi2Graph = linkErrors[IntPair(-1,0)];
        int   activeEdgeCount = linkErrors[IntPair(-1,-1)];

        double allLinksError = linkErrors[IntPair(-2,0)];
        if(activeLinks.size() != linkErrors[IntPair(-2,-1)])
        {
            cout<<"activeLinks.size(): "<<activeLinks.size()<<" linkErrors[IntPair(-2,0)]: "<<linkErrors[IntPair(-2,0)]<<endl;
            printf("This error about start node and end node is in %s on line %d\n",  __FILE__, __LINE__);
            exit(0);
        }

        IntPairSet::iterator it = activeLinks.begin(), end = activeLinks.end();
        std::vector<bool> GoodBitForAllLoop;

        int thisLinkClusterID;
        int allPass = 1;
        //put the error to map
        for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
             it!=end;
             it++)
        {
            if(it->first.first >=0)
            {
                //build the good loop vertor , to find the good region
                if(it->second < utils::chi2(3, goodLoopBelief))
				{
					// GoodBitForAllLoop.push_back(1);
//					cout<<"good loop"<<endl;
				}
                else
				{
					// delete the cluster id from good clusters set
					if(clusterizer.getClusterID(it->first) != clusterizer.clusterSizeVector[0].first)
					{
						good_clusters.erase(clusterizer.getClusterID(it->first));
					}
				}

            }
        }

        double ID_of_good_loop;
        for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
             it!=end;
             it++)
        {
            if(it->first.first >=0)
            {
                //put all loop and its error into a vector, to check the loop with big error
                pairLoopErr.first = it->first;
                pairLoopErr.second = it->second;
                VectorLoopError.push_back(pairLoopErr);
                //build the good loop vertor , to find the good region
                if(it->second < utils::chi2(3, goodLoopBelief))
				{
					ID_of_good_loop = clusterizer.getClusterID(it->first);
					//a good cluster should from 
					if(good_clusters.find(ID_of_good_loop) != good_clusters.end())
						GoodBitForAllLoop.push_back(1);
//					cout<<"good loop"<<endl;
				}

                else
				{
//                	cout<<"bad looop"<<endl;
					GoodBitForAllLoop.push_back(0);
				}

                clusterizer.getClusterID_multiExist(it->first, multiClusterExist);
                if(multiClusterExist.size() == 0)
                {
                    thisLinkClusterID = clusterizer.getClusterID(it->first);

                    cout<<"loop "<<(it->first).first<<" "<<(it->first).second<<" in cluster "<<thisLinkClusterID<<" has error "<<it->second<<endl;

                    mapCluserIDLoops[thisLinkClusterID].push_back(it->first);
                    mapCluserIDLoopErrors[thisLinkClusterID].push_back(it->second);
                }
                else
                {
                    cout<<"loop "<<(it->first).first<<" "<<(it->first).second<<" in cluster ** has error "<<it->second<<" "<<" multiExist:";
                    for(int multiExist = 0; multiExist < multiClusterExist.size(); multiExist++)
                    {
                        thisLinkClusterID = multiClusterExist[multiExist];
                        cout<<thisLinkClusterID<<" ";

                        mapCluserIDLoops[thisLinkClusterID].push_back(it->first);
                        mapCluserIDLoopErrors[thisLinkClusterID].push_back(it->second);
                    }
                    cout<<endl;
                }
            }
        }

        int nextGood, previousGood;
        double dis, dis1, dis2;
        std::vector< std::pair<int,int> > goodLoopRegin;
        findGoodLoopRegin(GoodBitForAllLoop, goodLoopRegin);
        if(goodLoopBelief <= 0.95)
        {
        	if(goodLoopRegin.size() == 0)
        	{
        		cout<<"goodLoopBelief == "<<goodLoopBelief<<" but the good region size is 0,so delete all"<<endl;
        		H.clear();
        	}
        }
        //handle each big error based on its cluster value
        deleteSome = 0;
        
        if(goodLoopRegin.size() != 0)
        {
            for(int iterator = 0; iterator < VectorLoopError.size(); iterator++)
            {
                if(VectorLoopError[iterator].second > utils::chi2(3, beliefResidualError))
                {
                	double smallestDis = 200;
                    int NumIteratorForDis = 0;
                    //find the nearest good looop
                    int middleIndex, intindexIndex;
                    if(iterator < goodLoopRegin[0].first)
                        intindexIndex = 0;
                    else if(iterator > goodLoopRegin.back().first)
                        intindexIndex = goodLoopRegin.size()-1;
                    else
                    {
                        for(int loopRegion = 0; loopRegion < goodLoopRegin.size()-1; loopRegion++)
                        {
                            if(iterator > goodLoopRegin[loopRegion].first and
                               iterator < goodLoopRegin[loopRegion+1].first)
                            {
                                intindexIndex = loopRegion;
                                break;
                            }
                        }
                    }


                    cout<<" "<<endl;
                    cout<<"loop "<<(VectorLoopError[iterator].first).first<<" "
                        <<(VectorLoopError[iterator].first).second
                        <<" has error "<<VectorLoopError[iterator].second<<endl;
                    dis = 0.0;

                    if(intindexIndex == goodLoopRegin.size()-1)
                    {
                        int GodRepeatTimes = goodLoopRegin[intindexIndex].second;
                        int loopTimes =  min(5, GodRepeatTimes);
                        for(int iter = 0; iter < loopTimes; iter++)
                        {
                            middleIndex = goodLoopRegin[intindexIndex].first+GodRepeatTimes-1-iter;
                            dis1 = clusterizer.calTransformDis_loop_loop(
                                    VectorLoopError[iterator].first, VectorLoopError[middleIndex].first);
                            if(dis1 > dis)
                                dis = dis1;
                            if(dis1 > 0 and dis1 < smallestDis)
                            	smallestDis = dis1;
                            cout<<"the last good region start from "<<goodLoopRegin[intindexIndex].first<<", repeats "<<GodRepeatTimes<<" times"<<endl;
                            int sdf = goodLoopRegin[intindexIndex].first;
                            cout<<"start loop is "<<VectorLoopError[sdf].first.first<<" "<<VectorLoopError[sdf].first.second<<endl;
                            cout<<"next good loop is "<<VectorLoopError[middleIndex].first.first<<" "<<
                                VectorLoopError[middleIndex].first.second<<endl;
                            cout<<"transform distacne is "<<dis1<<endl;
                            cout<<"smallest  distacne is "<<dis1<<endl;
                        }
                        if(GodRepeatTimes > 5)
                        {
	                        for(int iter = 0; iter < loopTimes; iter++)
	                        {
	                            middleIndex = goodLoopRegin[intindexIndex].first+iter;
	                            dis1 = clusterizer.calTransformDis_loop_loop(
	                                    VectorLoopError[iterator].first, VectorLoopError[middleIndex].first);
	                            if(dis1 > dis)
	                                dis = dis1;
	                            if(dis1 > 0 and dis1 < smallestDis)
	                            	smallestDis = dis1;
	                            cout<<"the last good region start from "<<goodLoopRegin[intindexIndex].first<<", repeats "<<GodRepeatTimes<<" times"<<endl;
	                            int sdf = goodLoopRegin[intindexIndex].first;
	                            cout<<"start loop is "<<VectorLoopError[sdf].first.first<<" "<<VectorLoopError[sdf].first.second<<endl;
	                            cout<<"next good loop is "<<VectorLoopError[middleIndex].first.first<<" "<<
	                                VectorLoopError[middleIndex].first.second<<endl;
	                            cout<<"transform distacne is "<<dis1<<endl;
	                            cout<<"smallest  distacne is "<<dis1<<endl;
	                        } 
                        }
                       
                    }
                    else
                    {
                        int GodRepeatTimes = goodLoopRegin[intindexIndex].second;
                        int loopTimes =  min(5, GodRepeatTimes);
                        for(int iter = 0; iter < loopTimes; iter++)
                        {
                            middleIndex = goodLoopRegin[intindexIndex].first+GodRepeatTimes-1-iter;
                            dis1 = clusterizer.calTransformDis_loop_loop(
                                    VectorLoopError[iterator].first, VectorLoopError[middleIndex].first);
                            if(dis1 > dis)
                                dis = dis1;
                            cout<<"next good loop is "<<VectorLoopError[middleIndex].first.first<<" "<<
                                VectorLoopError[middleIndex].first.second<<endl;
                            cout<<"transform distacne is "<<dis1<<endl;
                            if(dis1 > 0 and dis1 < smallestDis)
                            	smallestDis = dis1;    
                            cout<<"smallest  distacne is "<<dis1<<endl;                        
                        }
                        if(GodRepeatTimes > 5)
                        {
	                        for(int iter = 0; iter < loopTimes; iter++)
	                        {
	                            middleIndex = goodLoopRegin[intindexIndex].first+iter;
	                            dis1 = clusterizer.calTransformDis_loop_loop(
	                                    VectorLoopError[iterator].first, VectorLoopError[middleIndex].first);
	                            if(dis1 > dis)
	                                dis = dis1;
	                            if(dis1 > 0 and dis1 < smallestDis)
	                            	smallestDis = dis1;
	                            cout<<"the last good region start from "<<goodLoopRegin[intindexIndex].first<<", repeats "<<GodRepeatTimes<<" times"<<endl;
	                            int sdf = goodLoopRegin[intindexIndex].first;
	                            cout<<"start loop is "<<VectorLoopError[sdf].first.first<<" "<<VectorLoopError[sdf].first.second<<endl;
	                            cout<<"next good loop is "<<VectorLoopError[middleIndex].first.first<<" "<<
	                                VectorLoopError[middleIndex].first.second<<endl;
	                            cout<<"transform distacne is "<<dis1<<endl;
	                            cout<<"smallest  distacne is "<<dis1<<endl;
	                        } 
                        }


                        GodRepeatTimes = goodLoopRegin[intindexIndex+1].second;
                        loopTimes =  min(5, GodRepeatTimes);
                        for(int iter = 0; iter < loopTimes; iter++)
                        {
                            middleIndex = goodLoopRegin[intindexIndex+1].first+iter;
                            dis1 = clusterizer.calTransformDis_loop_loop(
                                    VectorLoopError[iterator].first, VectorLoopError[middleIndex].first);
                            if(dis1 > dis)
                                dis = dis1;
                            cout<<"next good loop is "<<VectorLoopError[middleIndex].first.first<<" "<<
                                VectorLoopError[middleIndex].first.second<<endl;
                            cout<<"transform distacne is "<<dis1<<endl;
                            if(dis1 > 0 and dis1 < smallestDis)
                            	smallestDis = dis1;       
                            cout<<"smallest  distacne is "<<dis1<<endl;                     
                        }
                        if(GodRepeatTimes > 5)
						{
	                        for(int iter = 0; iter < loopTimes; iter++)
	                        {
	                            middleIndex = goodLoopRegin[intindexIndex+1].first+iter;
	                            dis1 = clusterizer.calTransformDis_loop_loop(
	                                    VectorLoopError[iterator].first, VectorLoopError[middleIndex].first);
	                            if(dis1 > dis)
	                                dis = dis1;
	                            if(dis1 > 0 and dis1 < smallestDis)
	                            	smallestDis = dis1;
	                            cout<<"the last good region start from "<<goodLoopRegin[intindexIndex].first<<", repeats "<<GodRepeatTimes<<" times"<<endl;
	                            int sdf = goodLoopRegin[intindexIndex].first;
	                            cout<<"start loop is "<<VectorLoopError[sdf].first.first<<" "<<VectorLoopError[sdf].first.second<<endl;
	                            cout<<"next good loop is "<<VectorLoopError[middleIndex].first.first<<" "<<
	                                VectorLoopError[middleIndex].first.second<<endl;
	                            cout<<"transform distacne is "<<dis1<<endl;
	                            cout<<"smallest  distacne is "<<dis1<<endl;
	                        }  
						}                       
                    }


//                    //calculate the transform distance
//                    cout<<"next good loop is "<<VectorLoopError[middleIndex].first.first<<" "<<
//                        VectorLoopError[middleIndex].first.second<<endl;
//
//                    dis1 = clusterizer.calTransformDis_loop_loop(
//                            VectorLoopError[iterator].first, VectorLoopError[middleIndex].first);
//                    dis2 = clusterizer.calTransformDis_loop_loop(
//                            VectorLoopError[iterator].first, VectorLoopError[middleIndex+1].first);
//                    dis = std::max(dis1, dis2);

                    cout<<"transform distacne is "<<dis<<endl;

                    if(dis < utils::chi2(3, beliefTransError))
                    {
                    	cout<<"the biggest distance is "<<dis<<" so we continue"<<endl;
                    	continue;
                    }
                    //get the cluster ID of this loop
                    clusterizer.getClusterID_multiExist(VectorLoopError[iterator].first, multiClusterExist);
                    if(multiClusterExist.size() == 0 or multiClusterExist.size() == 1)
                    {
                    	if(multiClusterExist.size() == 0)
                        	thisLinkClusterID = clusterizer.getClusterID(VectorLoopError[iterator].first);
						else
							thisLinkClusterID = multiClusterExist[0];
                        cout<<"clusterID: "<<thisLinkClusterID<<endl;
                        //get the minimun error in the cluster
                        double minError = *std::min_element(mapCluserIDLoopErrors[thisLinkClusterID].begin(),
                                                            mapCluserIDLoopErrors[thisLinkClusterID].end());

                        //if min_err > chi(0.95), continue
                        if (minError > utils::chi2(3, 0.999) and 
                        	clusterizer.getClusterByID(thisLinkClusterID).size()==1 and
                        	dis > utils::chi2(3, beliefTransError))
                        {
                        	cout<<"the transform dis to the nearest good loop is "<<dis<<endl;
                            cout<<"the min error in this cluster is "<<minError<<" and its a one loop cluster, so delete cluster "<<thisLinkClusterID<<endl;
                            H.erase(thisLinkClusterID);
                            deleteSome = 1;
                            RejectSet.insert(thisLinkClusterID);
                            continue;
                        }
                            //if min_err < chi(0.95)
                        else
                        {
                            if(dis > utils::chi2(3, beliefTransError))
                            {
                                cout<<"dis "<<dis<<" bigger than belief "<<beliefTransError<<", so delete this loop; sigle exist"<<endl;
                                if(clusterizer.getClusterByID(thisLinkClusterID).size() == 1)
                                {
                                    deleteSome = 1;
                                    H.erase(thisLinkClusterID);
                                    RejectSet.insert(thisLinkClusterID);
                                }
                                else
                                {
                                    deleteSome = 1;
                                    clusterizer.setClusterID(VectorLoopError[iterator].first,ID_IGNORE);
                                    //erase form cluster id to loop map
                                    clusterizer.clusterIDtoLoopsMap[thisLinkClusterID].erase(VectorLoopError[iterator].first);
                                }
                            }
                        }
                    }
                    else
                    {
                        //print the clusters contain this loop
                        // cout<<"loop "<<(it->first).first<<" "<<(it->first).second<<" in cluster ** has error "<<it->second<<" "<<" multiExist:";
                        cout<<"clusterID:";
                        for(int multiExist = 0; multiExist < multiClusterExist.size(); multiExist++)
                        {
                            thisLinkClusterID = multiClusterExist[multiExist];
                            cout<<thisLinkClusterID<<" ";
                        }
                        cout<<endl;
                        //iterator the clusters constain this loop to determin if this loop should be deleted away
                        for(int multiExist = 0; multiExist < multiClusterExist.size(); multiExist++)
                        {
                            thisLinkClusterID = multiClusterExist[multiExist];
                            //get the minimun error in the cluster
                            double minError = *std::min_element(mapCluserIDLoopErrors[thisLinkClusterID].begin(),
                                                                mapCluserIDLoopErrors[thisLinkClusterID].end());

//                            //if min_err > chi(0.95), continue
//                            if (minError > utils::chi2(3, 0.95))
//                            {
//
//                                cout<<"the mininum error "<<minError<<" in this cluster is bigger than %95, so continue"<<endl;
//                                continue;
//                            }
//
//                                //if min_err < chi(0.95)
//                            else
                            {
                                if(dis > utils::chi2(3, beliefTransError))
                                {
                                    cout<<"dis bigger than "<<beliefTransError<<", so delete this loop; multi exist"<<endl;
                                    if(clusterizer.getClusterByID(thisLinkClusterID).size() == 1)
                                    {
                                        deleteSome = 1;
                                        H.erase(thisLinkClusterID);
                                        RejectSet.insert(thisLinkClusterID);
                                    }
                                    else
                                    {
                                        deleteSome = 1;
                                        clusterizer.clusterIDtoLoopsMap[thisLinkClusterID].erase(VectorLoopError[iterator].first);
                                        clusterizer.deleteClusterID_multiExist(VectorLoopError[iterator].first, thisLinkClusterID);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    void showClusterIDandNum(std::set<int> & consistentClusters)
    {
        //print the clusters that suvived
        cout<<"suvived clusters:"<<endl;
        if(consistentClusters.size() == 0)
        {
        	cout<<"no cluster suvived"<<endl;
        }
        else
        {
	        // for(std::set<int>::iterator  conflict_check_iter = consistentClusters.begin(); conflict_check_iter != consistentClusters.end(); conflict_check_iter++)
	        for(std::set<int>::iterator  conflict_check_iter = consistentClusters.begin(); conflict_check_iter != consistentClusters.end(); conflict_check_iter++)
	        {
	            //get the cluster serial survived
	            //use  all_conflict_cluster to check
	            cout<<(*conflict_check_iter)<<"("<<clusterizer.getClusterByID(*conflict_check_iter).size()<<") "<<endl;
	        }
	        cout<<endl;
        }

    }

	bool robustify_simplify(bool eraseIncorrectLinks=false)
	{
		std::set<int> potentialGoodCluster;
		double clusterSize = clusterizer.clusterCount();

		for(int i = 0; i < clusterSize; i++)
		{
			potentialGoodCluster.insert(i);
		}

		if(potentialGoodCluster.size() != clusterizer.clusterSizeVector.size())
		{
			cout<<"two sizes not equal, exit"<<endl;
			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		if(!gWrapper->isInitialized())
		{
			std::cerr<<" Please read in a graph file with read() or set the optimizer with setOptimizer()"<<std::endl;
			return false;
		}
		// First look for intra-cluster consistency

		IntSet
			// consistentClusters,
			hypotheses,
			hypotheses_largeCluster,
			hypotheses_oneEleCluster,
			goodSet,goodSetForLarge,
			rejectSet,rejectSetForLarge,
			tempRejectSet;
		std::vector<int> consistentClusters;
			
		std::vector<std::pair<int, int> > loopNUmbersInConsistentClusters;
		std::vector<std::pair<int,int> > reAccept;
		
		std::cout<<"Number of Clusters found : "<<clusterSize<<std::endl;
		std::cout<<"Checking Intra cluster consistency : "<<std::endl;
		//save all loops information once the cluster survive
		ofstream fileStreamr; 
		std::pair<g2o::SE2, Matrix3d> tSave;
		// std::pair<int,int> ty;
		int  originalLoopNum;
		bool passIntrachenck_second, passIntrachenck_active = 1;
		std::vector<std::pair<double, std::pair<int, int> > >  chiStstis4eachLoop;
		std::vector<std::pair<int, int> >  doubtLoopVector;
		std::vector<int> clusterIDofDoubtLoops, ToDeleteFromConflict;

		std::vector<int> acitveEffectCount, secondEffectCount;


		struct timeval t1, t2;
		gettimeofday(&t1, NULL);

// 		for(size_t i=0; i< clusterizer.clusterCount(); i++)
// 		{
// 			cout<<" "<<endl;
// 			std::cout<<i<<endl;
// 			std::vector<double>  chiStatisVector;
// 			std::vector<std::pair<int, int> >  doubtLoopNumberVector;
		
// 			//this one find good set, and have stritc for one menber cluster and strict belief for single loop, previous &95, now is %70;	
			// passIntrachenck_second = intraClusterConsistent_goodSET_handleOneMenCluster(i, chiStatisVector, doubtLoopNumberVector, originalLoopNum, chiStstis4eachLoop);
			
// // clusterizer.clusterSizeVector
// 		// std::set<std::pair<int,int > >::iterator sta = conflictPairSet.begin(), enD = conflictPairSet.end();
// 		// for(int num = 0; num < clusterSizeVector.size();num++)
// 		// {
// 		// 	sta = conflictPairSet.begin(), enD = conflictPairSet.end();
// 		// 	for(; sta != enD; sta++)
// 		// 	{
// 		// 		if((*sta).first == clusterSizeVector[num].first)
// 		// 			cout<<"cluster "<<(*sta).first<<" has conflict "<<(*sta).second<<endl;
// 		// 	}
			
// 		// }

// 			passIntrachenck_second = complementary_intra(i, chiStatisVector, doubtLoopNumberVector, originalLoopNum, chiStstis4eachLoop);
			
// 			// std::cin.get();

// 			//save the suvived clusters that pass intra check
// 			if(passIntrachenck_second)
// 			{

// 				consistentClusters.push_back(i);
// 				loopNUmbersInConsistentClusters.push_back(std::pair<int, int> (i, originalLoopNum));
// 			}
// 		}


		double clusterToTest;
		std::vector<double>  chiStatisVector;
		std::vector<std::pair<int, int> > deleted_loops;
		std::vector<std::pair<int, int> >  doubtLoopNumberVector;
		// if (clusterizer.clusterSizeVector.size() > 1)
		// 	clusterToTest = clusterizer.clusterSizeVector[1].first;
		// else
		// 	clusterToTest = clusterizer.clusterSizeVector[0].first;

		// 	passIntrachenck_second = complementary_intra_second(clusterToTest, chiStatisVector, 
		// 		doubtLoopNumberVector, originalLoopNum, chiStstis4eachLoop, deleted_loops);			
			
		// 	cin.get();

		
		for(size_t i=0; i< clusterSize; i++)
		{
			cout<<" "<<endl;
			std::cout<<i<<endl;
			chiStatisVector.clear();
			deleted_loops.clear();
			doubtLoopNumberVector.clear();
			clusterToTest = clusterizer.clusterSizeVector[i].first;
			cout<<"cluster "<<clusterToTest<<endl;
		
			//this one find good set, and have stritc for one menber cluster and strict belief for single loop, previous &95, now is %70;	
			// passIntrachenck_second = intraClusterConsistent_goodSET_handleOneMenCluster(i, chiStatisVector, doubtLoopNumberVector, originalLoopNum, chiStstis4eachLoop);
			
// clusterizer.clusterSizeVector
		// std::set<std::pair<int,int > >::iterator sta = conflictPairSet.begin(), enD = conflictPairSet.end();
		// for(int num = 0; num < clusterSizeVector.size();num++)
		// {
		// 	sta = conflictPairSet.begin(), enD = conflictPairSet.end();
		// 	for(; sta != enD; sta++)
		// 	{
		// 		if((*sta).first == clusterSizeVector[num].first)
		// 			cout<<"cluster "<<(*sta).first<<" has conflict "<<(*sta).second<<endl;
		// 	}
			
		// // }


			if(potentialGoodCluster.find(clusterToTest) == potentialGoodCluster.end())
			{
				// cout<<"can not find cluster "<< clusterToTest <<" in potentialGoodCluster, so continue"<<endl;
				continue;
			}

			// the complementary_intra() use robust kernel to optimize, and delete the loops whose non-robust error is larger than %95
			// cout<<"begin to run complementary_intra"<<endl;
			passIntrachenck_second = complementary_intra(clusterToTest, chiStatisVector, 
				doubtLoopNumberVector, originalLoopNum, chiStstis4eachLoop, deleted_loops);
			// cout<<"if you want to check the next cluster, press any key to continue"<<endl;
			// cin.get();

			// the third version use robust kernel to opeimize, then select the loops whose non-robust error is alrge than %95
			// as double one, then check their tranform distance.
			// passIntrachenck_second = complementary_intra_third(clusterToTest, chiStatisVector, 
			// 	doubtLoopNumberVector, originalLoopNum, chiStstis4eachLoop, deleted_loops);


			// passIntrachenck_second = complementary_intra_second(clusterToTest, chiStatisVector, 
			// 	doubtLoopNumberVector, originalLoopNum, chiStstis4eachLoop, deleted_loops);			
			
			// std::cin.get();

			//save the suvived clusters that pass intra check
			if(passIntrachenck_second)
			{
				// if this cluster size is bigger than 1
				if(clusterizer.clusterSizeVector[i].second > 1)
				{
					std::pair<int, int> Pair, Pair_reverse;
					Pair.first = clusterToTest;
					Pair_reverse.second = clusterToTest;
					// test each one member clsuter
					for(int findOne = 0; findOne < clusterSize; findOne++)
					{
						// if(clusterizer.clusterSizeVector[findOne].second >= 2)
						// 	continue;
						Pair.second = clusterizer.clusterSizeVector[findOne].first;
						Pair_reverse.first = clusterizer.clusterSizeVector[findOne].first;
						// delete those one member clsuter that is conflict with this good clsuter
						if(clusterizer.conflictPairSet.find(Pair) != clusterizer.conflictPairSet.end())
						{

							// std::set<std::pair<int, int> >::iterator sta_l = clusterizer.conflict_cause[Pair].begin(),
							// 	end_l = clusterizer.conflict_cause[Pair].end();
							// for(; 
							// 	sta_l != end_l;
							// 	std_l++)

							// bool wrong_conflict = 0;

							// for(int check_wrong_conflict = 0; 
							// 		check_wrong_conflict < int(deleted_loops.size()); 
							// 		check_wrong_conflict++)
							// {
							// 	// if(clusterizer.conflit_clsuter_pair_map_cause_loops.find(Pair) == 
							// 	// 	clusterizer.conflit_clsuter_pair_map_cause_loops.end())
							// 	// {
							// 	// 	cout<<"shold find this key but not"<<endl;
							// 	// 	exit(0);
							// 	// }
							// 	// cout<<"erase wrong cause loop "<<deleted_loops[check_wrong_conflict].first<<" "<<
							// 	// 	deleted_loops[check_wrong_conflict].second<<endl;
							// 	clusterizer.conflit_clsuter_pair_map_cause_loops[Pair].erase(deleted_loops[check_wrong_conflict]);
							// 	clusterizer.conflit_clsuter_pair_map_cause_loops[Pair_reverse].erase(deleted_loops[check_wrong_conflict]);
							// }
							// if(clusterizer.conflit_clsuter_pair_map_cause_loops[Pair].size() != 0)
							{
								cout<<"cluster "<<Pair.second<<
									" is deleted as it is conflict with cluster "<<Pair.first<<endl;
								potentialGoodCluster.erase(Pair.second);
							}

						}
					}
				}
				cout<<"clsuter "<<clusterToTest<<" suvived "<<endl;
				// cin.get();
				consistentClusters.push_back(clusterToTest);
				loopNUmbersInConsistentClusters.push_back(std::pair<int, int> (clusterToTest, originalLoopNum));
			}
		}

		std::cout<<"consistentClusters size: "<<consistentClusters.size()<<std::endl;
		// sleep(2);
		// save the clusters after intra
		fileStreamr.open("clusters_after_intra.txt",ios::trunc);
		std::array<double,6> ty={1,1,1,1,1,1};
		std::pair<int, int> ty_pair;
		for(int save_cluster = 0; save_cluster < clusterizer.clusterSizeVector.size(); save_cluster++)
		{
			int cluster_id = clusterizer.clusterSizeVector[save_cluster].first;
			if(find_ele.find(consistentClusters, cluster_id) < 0)
				continue;
			if(clusterizer.mixed_clusters.find(cluster_id) != clusterizer.mixed_clusters.end())
				fileStreamr<<"clusterID "<<cluster_id<<" "<<cluster_id+1<<" mixed"<<"\n";
			else
				fileStreamr<<"clusterID "<<cluster_id<<" "<<cluster_id+1<<"\n";

			IntPairSet& currentCluster = clusterizer.getClusterByID(cluster_id);
			IntPairSet::const_iterator itVertex = currentCluster.begin(), 
				lendVertex = currentCluster.end();
			for(;itVertex!=lendVertex;itVertex++)
			{
				ty_pair = *itVertex;

				tSave = clusterizer.LP_Trans_Covar_Map[ty_pair];

				bool bad  = 0;

				for(int i = 0; i <clusterizer.bad_loops_vector.size(); i++)
				{
					if(clusterizer.bad_loops_vector[i].first == ty[0] and clusterizer.bad_loops_vector[i].second == ty[3])
					{
						bad = 1;
						break;
					}
				}

				fileStreamr<<"EDGE_SE2 "<<ty_pair.first<<" "<<ty_pair.second<<" "<<tSave.first[0]<<" "<<tSave.first[1]<<" "<<tSave.first[2]
					<<" "<<1.0/tSave.second(0,0)<<" "<<0<<" "<<0<<" "<<1.0/tSave.second(1,1)<<" "<<0<<" "<<1.0/tSave.second(2,2)<<" "<<bad<<"\n";
				// fileStream<<trystdarray[0]<<"\n";
			}	
		}
		fileStreamr.close();		

		gettimeofday(&t2, NULL);
		//那么函数f运行所花的时间为
		deltaT = (t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec) / 1000000;// 微秒
		cout<<"*********************"<<endl;
		cout<<"*********************"<<endl;
		cout<<"*********************"<<endl;
		cout<<"intra time consumed: "<<deltaT<<endl;
		cout<<"*********************"<<endl;
		cout<<"*********************"<<endl;
		cout<<"*********************"<<endl;

		fileStreamr.close();


		std::set<int> suvivedTointerSet(consistentClusters.begin(), consistentClusters.end());
		std::set<int> suvivedTointerSet_backup(consistentClusters.begin(), consistentClusters.end());
		//pick out the cluster whose has more then one members
		// std::set<std::set<int> > mergeSet;
		IntPairSet activeLoops;
		// IntSet::iterator
		std::vector<int>::iterator
			cIt 	= consistentClusters.begin(),
			cEnd 	= consistentClusters.end();

		// cout<<"big consistentClusters (whos loop number is bigger than 1) : "<<endl;

		for( ; cIt!=cEnd; cIt++)
		{
			if(clusterizer.getClusterByID(*cIt).size() > 1)
				hypotheses_largeCluster.insert(*cIt);
			else
			{
				hypotheses_oneEleCluster.insert(*cIt);
			}
		}

		//chenck conflict info before inter
		std::vector<std::pair<int, int> > confPairVector;
		std::set<int > confClusterSet;

		// cout<<" "<<endl;
		// cout<<"start to check if cluseter pair in good set has confict"<<endl;

		cout<<"after intra"<<endl;
        showClusterIDandNum(suvivedTointerSet);
		if(suvivedTointerSet.size() > 1)
		{
			bool conf =0;
			std::vector<int > temporary_vector_int(suvivedTointerSet.begin(), suvivedTointerSet.end());
			std::pair<int, int> pair;
			for(int i = 0; i < temporary_vector_int.size()-1; i++)
			{
				for(int j = i+1; j < temporary_vector_int.size(); j++)
				{
					pair.first = temporary_vector_int[i];
					pair.second = temporary_vector_int[j];
					if(clusterizer.conflictPairSet.find(pair) != clusterizer.conflictPairSet.end())
					{
						cout<<"conflict clusters pair: "<<pair.first<<" "<<pair.second<<endl;
						conf = 1;
					}
				}
				
			}
			if (conf == 0)
			{
				cout<<"no conflict pair "<<endl;
			}
		}
		// cout<<"intra is just finished"<<endl;
		// cin.get();
		// // std::set<int> suvivedTointerSet_backup_sec(suvivedToInter.begin(), suvivedToInter.end());

		gettimeofday(&t1, NULL);
		std::map<int,double>  errorEachCluster;
		IntPairSet  activeLinks;
	

		if(suvivedTointerSet.size() == 0)
		{
			cout<<"no loop suvived"<<endl;
			fileStreamr.open("clusters_after_intra_multi.txt",ios::trunc);
			fileStreamr.close();
			if(clusterizer.NumOfRealLoops != -1)
			{		
				IntPairSet::iterator iit = clusterizer.set4GTL.begin(), iiend = clusterizer.set4GTL.end();
			
				for(; iit != iiend; iit++)
				{

					lostGood = lostGood+1;
					cout<<"true loop "<<(*iit).first <<" "<<(*iit).second<< " in cluster "<< clusterizer.getClusterID(*iit)<< " is abandened"<<endl;

				}
			}
			return 1;
		}
		activeLinks.clear();
		gatherLinks(suvivedTointerSet,activeLinks);
		bool final_suvived = complementary_intra_multi(activeLinks, chiStatisVector, 
						doubtLoopNumberVector, originalLoopNum, chiStstis4eachLoop, suvivedTointerSet);	

		// save clsuters after complementary multi
		std::set<int>::iterator
			cItSet	= suvivedTointerSet.begin(),
			cEndSet = suvivedTointerSet.end();
	
		fileStreamr.open("clusters_after_intra_multi.txt",ios::trunc);
		for(int save_cluster = 0; save_cluster < clusterizer.clusterSizeVector.size(); save_cluster++)
		{

			int cluster_id = clusterizer.clusterSizeVector[save_cluster].first;
			if(suvivedTointerSet.find(cluster_id) == suvivedTointerSet.end())
				continue;
			if(clusterizer.mixed_clusters.find(cluster_id) != clusterizer.mixed_clusters.end())
				fileStreamr<<"clusterID "<<cluster_id<<" "<<cluster_id+1<<" mixed"<<"\n";
			else
				fileStreamr<<"clusterID "<<cluster_id<<" "<<cluster_id+1<<"\n";
			IntPairSet& currentCluster = clusterizer.getClusterByID(cluster_id);
			IntPairSet::const_iterator itVertex = currentCluster.begin(), 
				lendVertex = currentCluster.end();
			for(;itVertex!=lendVertex;itVertex++)
			{
				tSave = clusterizer.LP_Trans_Covar_Map[*itVertex];

				bool bad  = 0;

				for(int i = 0; i <clusterizer.bad_loops_vector.size(); i++)
				{
					if(clusterizer.bad_loops_vector[i].first == (*itVertex).first and 
						clusterizer.bad_loops_vector[i].second == (*itVertex).second)
					{
						bad = 1;
						break;
					}
				}

				fileStreamr<<"EDGE_SE2 "<<(*itVertex).first<<" "<<(*itVertex).second<<" "<<tSave.first[0]<<" "<<tSave.first[1]<<" "<<tSave.first[2]
					<<" "<<1.0/tSave.second(0,0)<<" "<<0<<" "<<0<<" "<<1.0/tSave.second(1,1)<<" "<<0<<" "<<1.0/tSave.second(2,2)<<" "<<bad<<"\n";
				// fileStream<<trystdarray[0]<<"\n";
			}	
		}
		fileStreamr.close();	

		cout<<"after multi "<<endl;
		showClusterIDandNum(suvivedTointerSet);				
		if(final_suvived = 0)
		{
			suvivedTointerSet.clear();
		}

	
		std::cout<<" GoodSet before inter:";
		showClusterIDandNum(suvivedTointerSet);		

		interClusterConsistent( suvivedTointerSet, goodSet, rejectSet);
		// goodSet.insert(suvivedTointerSet.begin(),suvivedTointerSet.end());
		
		if(goodSet.size() > 1)
		{
			bool conf =0;
			std::vector<int > temporary_vector_int(goodSet.begin(), goodSet.end());
			std::pair<int, int> pair;
			for(int i = 0; i < temporary_vector_int.size()-1; i++)
			{
				for(int j = i+1; j < temporary_vector_int.size(); j++)
				{
					pair.first = temporary_vector_int[i];
					pair.second = temporary_vector_int[j];
					if(clusterizer.conflictPairSet.find(pair) != clusterizer.conflictPairSet.end())
					{
						cout<<"conflict clusters pair: "<<pair.first<<" "<<pair.second<<endl;
						conf = 1;
					}
				}
				
			}
			if (conf == 0)
			{
				cout<<"no conflict pair "<<endl;
			}
		}


		
		// int lastDcLUSTER = -1;
		// double meanErrorChi = 0;
		// interClusterConsistent_retrace( suvivedTointerSet, goodSet, rejectSet, lastDcLUSTER, meanErrorChi);

		gettimeofday(&t2, NULL);

		std::cout<<" GoodSet final:";
		for( IntSet::iterator it= goodSet.begin(), end = goodSet.end(); it!=end; it++)
		{
			std::cout<<(*it)<<" ";
		}
		std::cout<<std::endl;
		std::cerr<<" total loop closure in goodSet is: "<<getNumberIn(goodSet)<<std::endl;

		//那么函数f运行所花的时间为
		deltaT = (t2.tv_sec-t1.tv_sec) * 1000000 + t2.tv_usec-t1.tv_usec;// 微秒
		cout<<"*********************"<<endl;
		cout<<"*********************"<<endl;
		cout<<"*********************"<<endl;
		cout<<"inter time consumed: "<<deltaT<<endl;
		cout<<"*********************"<<endl;
		cout<<"*********************"<<endl;
		cout<<"*********************"<<endl;

		// interClusterConsistent_second( consistentClusters, goodSet, rejectSet);
		cout<<"number of clusters put into backend: "<<suvivedTointerSet_backup.size()<<endl;
		cout<<"number of clusters in good set: "<<goodSet.size()<<endl;
		cout<<"cluster in goodSET: ";

		cItSet	= goodSet.begin(),
		cEndSet = goodSet.end();
		for( ; cItSet!=cEndSet; cItSet++)//traverse all survived cluster
		{
			//*cIt is the serial number of the test cluster
			cout<<(*cItSet)<<" ";
		}
		cout<<endl;

		cout<<"cluster in suvivedTointerSet: ";

		cItSet	= suvivedTointerSet.begin(),
		cEndSet = suvivedTointerSet.end();
		for( ; cItSet!=cEndSet; cItSet++)//traverse all survived cluster
		{
			//*cIt is the serial number of the test cluster
			cout<<(*cItSet)<<" ";
		}
		cout<<endl;


		cout<<" "<<endl;
//		printf("This is in %s on line %d\n",  __FILE__, __LINE__);
		cout<<"final elemtns in good set:";		
		for(std::set<int>::iterator start = goodSet.begin(); start != goodSet.end(); start++)
		{
			cout<<*start<<" ";
		}
		cout<<endl;
		//intra and inter
		activeLoops.clear();


		gatherLinks(goodSet,activeLoops);
		//save the result
		const char *g2f_intra_inter=OptimizedResultIntraInter.data();
		// gWrapper->saveOptimizedResult(g2f_intra_inter, activeLoops);
		int numofGoodS = activeLoops.size();

        _goodSet.clear();
        _goodSet.insert(goodSet.begin(),goodSet.end());
        cout<<"_goood set size:"<<_goodSet.size()<<endl;
        cout<<"good set size:"<<goodSet.size()<<endl;

		//if there is no ground truth file, exit
		// if(clusterizer.set4GTL.size()==0)
		// 	return 1;

		int  forRecall=0;// number of the true loops that has been successfully picked out
		int goodLN = 0;  // number of good loops in all the picked out loops
		IntPairSet::iterator iit = activeLoops.begin(), iiend = activeLoops.end();
		

		//accuracy
		if(clusterizer.NumOfRealLoops != -1)
		{
			for(; iit != iiend; iit++)
			{
				if(clusterizer.set4GTL.find(*iit) == clusterizer.set4GTL.end())
				{
					allfalseAccepted = allfalseAccepted+1;
					cout<<"suvived loop "<<(*iit).first <<" "<<(*iit).second<< " in cluster "<< clusterizer.getClusterID(*iit)<<"is a flase one"<<endl;
				}
				else
					goodLN = goodLN+1;
			}
			if(allfalseAccepted != 0)
				cout<<" there are "<< allfalseAccepted<<" bad loops"<<endl;
			else
				cout<<"all loops suvived are true "<<endl;
		}

		//recall: find out the true loop that was abandened
		if(clusterizer.NumOfRealLoops != -1)
		{		
			// if(clusterizer.set4GTL.size() > activeLoops.size())
			// {
				iit = clusterizer.set4GTL.begin(), iiend = clusterizer.set4GTL.end();
			
				for(; iit != iiend; iit++)
				{
					if(activeLoops.find(*iit) == activeLoops.end())
					{
						lostGood = lostGood+1;
						cout<<"true loop "<<(*iit).first <<" "<<(*iit).second<< " in cluster "<< clusterizer.getClusterID(*iit)<< " is abandened"<<endl;
					}
					else
						forRecall = forRecall+1;
				}
			// }
		}

		//if there exists an ground truth date, then print out the relevent information
		cout<<" "<<endl;
		// if(clusterizer.NumOfRealLoops != -1)
		{
			cout<<" ******** information without recheck doubt loops ******** "<<endl;
			cout<<"the number of     suvived good loops    is: "<<activeLoops.size()<<std::endl;
			cout<<" "<<endl;
			cout<<"the number of         true  loops       is: "<< clusterizer.NumOfRealLoops<<" "<<endl;
			allTrue = clusterizer.NumOfRealLoops;
			cout<<" "<<endl;
			cout<<"the number of      loops in testfile    is: "<< clusterizer.LP_Trans_Covar_Map.size()<<" "<<endl;
			allTest = clusterizer.LP_Trans_Covar_Map.size();
			cout<<" "<<endl;
			if(activeLoops.size() > 0)
			{
				cout<<"the                  accuracy rate      is: "<<
					(goodLN)*1.0/(activeLoops.size())<<" "<<endl;
			}
			else
			{
				if(goodLN == 0)
				{				
					cout<<"the                  accuracy rate      is: "<< 1<<" "<<endl;
				}
				else
					cout<<"????????????"<<"goodLN is "<<goodLN<<" active loops nunber is "<<activeLoops.size()<<endl;
			} 
			// cout<<"goodLN "<< goodLN <<endl;
			cout<<" "<<endl;
			if(clusterizer.set4GTL.size() > 0)
			{
				cout<<"the                  recall  rate       is: "<< 
					double(forRecall*1.0)/(clusterizer.set4GTL.size())<<" "<<endl;
			}
			else
			{
				if(forRecall == 0)
				{
					cout<<"the                  recall  rate       is: "<< 1<<" "<<endl;					
				}
				else
					cout<<"????????????"<<"forRecall is "<<forRecall<<" active loops nunber is "<<clusterizer.set4GTL.size()<<endl;				
			}

				cout<<"forRecall: "<<forRecall<<endl;
			cout<<" ********************************************************* "<<endl;
			cout<<" "<<endl;
		}


		// if(allfalseAccepted != 0 or lostGood != 0)
		// {
		// 	cout<<"allfalseAccepted = "<<allfalseAccepted<<"   "<<"lostGood = "<<lostGood<<endl;
		// 	cin.get();
		// }

		return 1;


		// //recheck the doubt loops
		// std::vector<int> trueBitOfReacceptDoubtLoop;
		// if (doubtLoopVector.size() > 0)
		// {
		// 	if(clusterIDofDoubtLoops.size() != doubtLoopVector.size())
		// 	{
		// 			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
		// 			exit(0);
		// 	}
			
		// 	for(int doubt = 0; doubt < doubtLoopVector.size(); doubt++)
		// 	{
		// 		cout<<"doubt chenck "<<doubtLoopVector[doubt].first<<" "<<doubtLoopVector[doubt].second;
		// 		std::vector< std::pair<double, std::pair<std::pair<int, int>, int> > > & middCheck = 
		// 			clusterizer._clustersFound[clusterIDofDoubtLoops[doubt]].chiStatisWhenAdd2Cluster;;
		// 		// middCheck = clusterizer._clustersFound[clusterIDofDoubtLoops[doubt]].chiStatisWhenAdd2Cluster;
		// 		bool find = 0;
		// 		int j=0;
		// 		for(j =0; j< middCheck.size(); j++)
		// 		{
		// 			if(middCheck[j].second.first == doubtLoopVector[doubt])
		// 			{
		// 				find = 1;
		// 				cout<<" has transdis "<<middCheck[j].first<<endl;
		// 				if(middCheck[j].first < chi2_continuous(edgeDimension, 0.6))
		// 				{
		// 					reAccept.push_back(doubtLoopVector[doubt]);
		// 					//reput the loop into the cluster
		// 					clusterizer.setClusterID(doubtLoopVector[doubt],clusterIDofDoubtLoops[doubt]);
		// 					//check the reaccept are true or not

		// 					if(clusterizer.set4GTL.find(doubtLoopVector[doubt]) == clusterizer.set4GTL.end())
		// 						trueBitOfReacceptDoubtLoop.push_back(0);
		// 					else
		// 						trueBitOfReacceptDoubtLoop.push_back(1);

		// 				}
		// 				else
		// 				{
		// 					cout<<" this loop "<<doubtLoopVector[doubt].first<<" "<<doubtLoopVector[doubt].second<<" has been deleted"<<endl;
		// 					clusterizer.setClusterID(doubtLoopVector[doubt],ID_IGNORE);
		// 				}
		// 				break;
		// 			}
		// 		}
		// 		if(find == 0)
		// 		{
		// 			cout<<"the loop in doubt should be find out in the previous cluster, however, falied!"<<endl;
		// 			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
		// 			exit(0);
		// 		}

		// 		cout<<"doubt loop "<<doubt<<" is "<<doubtLoopVector[doubt].first<<" "<<doubtLoopVector[doubt].second<<
		// 			" ,  the chi statis is "<<middCheck[j].first<< " when it is added to the cluster "<<
		// 			clusterIDofDoubtLoops[doubt]<<endl;
		// 	}	
		// }
		// else
		// 	cout<<"no doube loop need to be checked"<<endl;

		// //print out the information of reaccept loops,
		// int addTrue = 0, addFalse = 0;
		// for(int reacce = 0; reacce < reAccept.size(); reacce++)
		// {
		// 	cout<<"reaccept doubt loop "<<reAccept[reacce].first<<" "<<reAccept[reacce].second;

		// 	if(trueBitOfReacceptDoubtLoop[reacce])
		// 	{
		// 		addTrue =addTrue +1;
		// 		cout<<" is true"<<endl;
		// 	}
		// 	else
		// 	{
		// 		addFalse =addFalse +1;
		// 		cout<<" is false"<<endl;
		// 	}
		// }


		// //if there exists an ground truth date, then print out the relevent information
		// cout<<" "<<endl;
		// if((clusterizer.NumOfRealLoops != -1) and (reAccept.size()>0))
		// {
		// 	cout<<" ***************  after recheck doubt loops *************** "<<endl;
		// 	cout<<"the number of     suvived good loops    is: "<<activeLoops.size()+addTrue+addFalse<<std::endl;
		// 	cout<<" "<<endl;
		// 	cout<<"the number of         true  loops       is: "<< clusterizer.NumOfRealLoops<<" "<<endl;
		// 	cout<<" "<<endl;
		// 	cout<<"the number of      loops in testfile    is: "<< clusterizer.LP_Trans_Covar_Map.size()<<" "<<endl;
		// 	cout<<" "<<endl;
		// 	cout<<"the                  accuracy rate      is: "<<
		// 		(goodLN+addTrue)*1.0/(activeLoops.size()+addTrue+addFalse)<<" "<<endl;
		// 	cout<<" "<<endl;
		// 	cout<<"the                  recall  rate       is: "<< 
		// 		(forRecall+addTrue)*1.0/(clusterizer.set4GTL.size())<<" "<<endl;
		// 	cout<<" ********************************************************* "<<endl;
		// 	cout<<" "<<endl;
		// }
		// else
		// 	cout<<"no change happen after recheck the doubt loops "<<endl;

		// //re get the loops in activeLoops after boubt check
		// activeLoops.clear();
		// gatherLinks(hypotheses,activeLoops);
		// 		//save the result file before doube check 
		// // gWrapper->optimize(activeLoops,nIterations, odoEdgeRelateLC_Error, odoEdgeError);

		// if((reAccept.size()>0))
		// {
		// 	const char *g2fre=OptimizedResultWithRecheck.data();
		// 	gWrapper->saveOptimizedResult(g2fre, activeLoops);
		// }
		// //save the cluster info after doubt check
		// // ofstream fileStreamr; 
		// fileStreamr.open("loops_info_of_all_suvived_clusters_after_intracheck_after_doubtcheck.txt",ios::trunc);

		// cItSet 	= hypotheses.begin(),
		// cEndSet 	= hypotheses.end();
		// cout<<"remained clusters in hypotheses: "<<endl;
		// for(; cItSet != cEndSet; cItSet++)
		// 	cout<<*cItSet<<" ";
		// cout<<endl;
		// cout<<" "<<endl;
		// cout<<"clusters in consistentClusters"<<endl;
		// cIt 	= consistentClusters.begin(),
		// cEnd 	= consistentClusters.end();
		// for(; cIt != cEnd; cIt++)
		// {
		// 	std::cout<<*cIt<<" "; std::cout.flush();
		// 	std::vector<double>  chiStatisVector;
		// 	std::vector<std::pair<int, int> >  doubtLoopNumberVector;

		// 		fileStreamr<<*cIt<<"\n";
		// 		//save loops information in i to a txt file
		// 		IntPairSet::iterator it = clusterizer.getClusterByID(*cIt).begin(), itend = clusterizer.getClusterByID(*cIt).end();
		// 		for(; it != itend; it++)
		// 		{
		// 			int first = (*it).first, toNode = (*it).second;
		// 			int j;
		// 			for(j = 0; j<clusterizer._clustersFound[*cIt].positionserial.size();j++)
		// 			{
		// 				if( first == clusterizer._clustersFound[*cIt].positionserial[j][0] and
		// 					toNode == clusterizer._clustersFound[*cIt].positionserial[j][3])
		// 				{
		// 					ty.first = clusterizer._clustersFound[*cIt].positionserial[j][0];
		// 					ty.second = clusterizer._clustersFound[*cIt].positionserial[j][3];

		// 					tSave = clusterizer.LP_Trans_Covar_Map[ty];
		// 					Matrix3d fsave;
		// 					fsave = tSave.second;

		// 					fileStreamr<<"EDGE_SE2 "<<ty.first<<" "<<ty.second<<" "<<tSave.first[0]<<" "<<tSave.first[1]<<" "<<tSave.first[2]
		// 						<<" "<<fsave.inverse()(0,0)<<" "<<fsave.inverse()(0,1)<<" "<<fsave.inverse()(0,2)<<" "<<
		// 						fsave.inverse()(1,1)<<" "<<fsave.inverse()(1,2)<<" "<<fsave.inverse()(2,2)<<" "<<"\n";	
		// 					break;
		// 				}
		// 				if (j == clusterizer._clustersFound[*cIt].positionserial.size()-1)
		// 				{
		// 					cout<<"can not find the postion info from the cluster for the suvived loop "<<first<<" "<<toNode <<endl;
		// 					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
		// 					exit(0);
		// 				}
	
		// 			}
		// 		}				
			
		// }

		// fileStreamr.close();
		// cout<<" "<<endl;
		// if((clusterizer.NumOfRealLoops != -1) and (reAccept.size()>0))
		// {

		// 	cout<<" ***************  after recheck doubt loops, to see if the reaccept loops have been reput into cluster *************** "<<endl;
		// 	cout<<"the number of     suvived good loops    is: "<<activeLoops.size()<<std::endl;
		// 	cout<<" "<<endl;
		// 	cout<<"the number of         true  loops       is: "<< clusterizer.NumOfRealLoops<<" "<<endl;
		// 	cout<<" "<<endl;
		// 	cout<<"the number of      loops in testfile    is: "<< clusterizer.LP_Trans_Covar_Map.size()<<" "<<endl;
		// 	cout<<" "<<endl;
		// 	cout<<"the                  accuracy rate      is: "<<
		// 		(goodLN+addTrue)*1.0/(activeLoops.size())<<" "<<endl;
		// 	cout<<" "<<endl;
		// 	cout<<"the                  recall  rate       is: "<< 
		// 		(forRecall+addTrue)*1.0/(clusterizer.set4GTL.size())<<" "<<endl;
		// 	cout<<" ********************************************************* "<<endl;
		// 	cout<<" "<<endl;
		// }
		// else
		// 	cout<<"no change happen after recheck the doubt loops "<<endl;
		// // sleep(3);
		// return true;
	}

    bool checkConflict(std::set<int> & goodSet, std::set<int > & confClusterSet, std::vector<int> &clustersDeletedInConflictCheck)
    {
		std::vector<std::pair<int, int> > confPairVector;
		

		cout<<"********************************************* "<<endl;
		cout<<"run checkConflict program"<<endl;
		cout<<"start to check if cluseter pair in good set has confict"<<endl;
		std::set<int>::iterator
			cItSet	= goodSet.begin(),
			cEndSet = goodSet.end();
		for( ; cItSet!=cEndSet; cItSet++)//traverse all survived cluster
		{
			//*cIt is the serial number of the test cluster
			cout<<"base cluster: "<< *cItSet<<endl;

			int serialLocal = find_ele.find(clusterizer.all_conflict_cluster, *cItSet);
			if(serialLocal != -1)
			{
				
				cout<<"has conflict clusters: ";
				for(int ITuncSET = 0; ITuncSET < clusterizer.all_conflict_cluster[serialLocal].second.size(); ITuncSET++)//tranverse each element in the uncertain set of the current cluster
				{	
					cout<<clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]<<" ";

				}
				cout<<endl;				
			}
			else
			{
				cout<<"has no conflict cluster"<<endl;
				cout<<" "<<endl;
			}

			if(serialLocal != -1)
			{			
				
				int haveSuvivedConflict, haveBit = 0;
				std::pair<int,int> toFindClusterPair, pairToCheck;
				for(int ITuncSET = 0; ITuncSET < clusterizer.all_conflict_cluster[serialLocal].second.size(); ITuncSET++)//tranverse each element in the uncertain set of the current cluster
				{	
					haveSuvivedConflict = find_ele.find(goodSet, clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]);
					if(haveSuvivedConflict != -1)
					{
						haveBit = 1;
						cout<<"suvived conflict clusters:"<<clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]<<" "<<endl;;

						if(*cItSet > clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET])
						{
							pairToCheck.first  = *cItSet;
							pairToCheck.second = clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET];
							if(find_ele.find(confPairVector, pairToCheck) == -1)
							{
								confPairVector.push_back(pairToCheck);
								confClusterSet.insert(pairToCheck.first);
								confClusterSet.insert(pairToCheck.second);
							}

						}

						//it is a map that contains conflict clsuter pair and  its distance when add to cluter
						toFindClusterPair = std::pair<int , int > (*cItSet, clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]);
						if(clusterizer.conflictClusterPairDistanceMap.find(toFindClusterPair) != 
							clusterizer.conflictClusterPairDistanceMap.end())
						{
							cout<<"the distance when form this conflict is "<<clusterizer.conflictClusterPairDistanceMap[toFindClusterPair]<<endl;
						}
						else
						{
							cout<<"can not find conflict info in previous pair map"<<endl;
							cout<<"the total number in the map is "<<clusterizer.conflictClusterPairDistanceMap.size()<<endl;
							printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);

							for(int ITuncSET_second = ITuncSET+1; ITuncSET_second < clusterizer.all_conflict_cluster[serialLocal].second.size(); ITuncSET_second++)
							{
								toFindClusterPair = std::pair<int , int > (*cItSet, clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET_second]);
								cout<<clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET_second]<<" conflict cluster"<<endl;;

								if(clusterizer.conflictClusterPairDistanceMap.find(toFindClusterPair) != 
									clusterizer.conflictClusterPairDistanceMap.end())
								{
									cout<<"the distance when form this conflict is "<<clusterizer.conflictClusterPairDistanceMap[toFindClusterPair]<<endl;
								}
								else
								{
									cout<<"can not find conflict info in previous pair map"<<endl;
									cout<<"the total number in the map is "<<clusterizer.conflictClusterPairDistanceMap.size()<<endl;
									printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
								}	
							}
							// exit(0);
						}
						cout<<"conflict exist after inter check"<<endl;
						// exit(0);
					}
				}
				if(haveBit == 0)
					cout<<"has no suvived conflict cluster"<<endl;
				cout<<endl;				
			}
				//check if the uncertian cluster also exit in the survived cluster
		}

		std::vector<std::pair<int, double> > conflictTimes;		
		std::pair<int, double> confTimeEle;		
		//calculate the conflict times for each conflict cluster
		for(std::set<int>::iterator its = confClusterSet.begin(); its != confClusterSet.end(); its++)
		{
			cout<<"conflict clsuter "<<*its;
			confTimeEle.first = *its;
			int times = 0;
			for(int i = confPairVector.size()-1; i >= 0; i--)
			{
				if(confPairVector[i].first == *its or confPairVector[i].second == *its)
					times =times+1;
			}
			confTimeEle.second = times;
			cout<<" has "<<times<<" times conflict"<<endl;;
			if(times == 0)
			{
				cout<<"times == 0"<<endl;
				exit(0);
			}
			conflictTimes.push_back(confTimeEle);
		}	
		if(conflictTimes.size() == 0)
		{
			cout<<"no conflict in good set***"<<endl;
			return 1;
		}
		if(conflictTimes.size() < 2)
		{
			cout<<"conflictTimes.size < 0*****"<<endl;
			exit(0);
		}

		int conflictTime;
		while(confPairVector.size() != 0)
		{
			//find the cluster has biggest conflict times
			conflictTime = 0;
			for(int i= 1; i < conflictTimes.size(); i++)
			{
				if(conflictTimes[i].second > conflictTimes[conflictTime].second)
					conflictTime = i;
			}

			if(conflictTimes[conflictTime].second == 1)
				break;
			//delete the cluster pair that has the biggest times clsuterID
			if(conflictTimes[conflictTime].second > 1)
			{
				//deletet the cluster that has biggest conflict times from good set
				clustersDeletedInConflictCheck.push_back(conflictTimes[conflictTime].first);
				goodSet.erase(conflictTimes[conflictTime].first);
				cout<<"1 delete clsuter "<<conflictTimes[conflictTime].first<<endl;
				for(int i = confPairVector.size()-1; i >= 0; i--)
				{
					if(confPairVector[i].first == conflictTimes[conflictTime].first or confPairVector[i].second == conflictTimes[conflictTime].first)
						confPairVector.erase(confPairVector.begin()+i);
					cout<<"confPairVector.size: "<<confPairVector.size()<<endl;
				}
				// cout<<"confPairVector.size: "<<confPairVector.size()<<endl;
			}
			if(confPairVector.size() == 0)
				return 1;
			//update the conflict clsuter set
			confClusterSet.clear();
			for(int i = confPairVector.size()-1; i >= 0; i--)
			{
				confClusterSet.insert(confPairVector[i].first);
				confClusterSet.insert(confPairVector[i].second);
			}
			//re calcualtel the conflict times for each conflict cluster
			conflictTimes.clear();
			for(std::set<int>::iterator its = confClusterSet.begin(); its != confClusterSet.end(); its++)
			{
				confTimeEle.first = *its;
				int times = 0;
				for(int i = confPairVector.size()-1; i >= 0; i--)
				{
					if(confPairVector[i].first == *its or confPairVector[i].second == *its)
						times =times+1;
				}
				confTimeEle.second = times;
				if(times == 0)
				{
					cout<<"times == 0"<<endl;
					exit(0);
				}
				conflictTimes.push_back(confTimeEle);
			}

		}

		if(confPairVector.size()!= 0 )
		{
			for(int i= 1; i < conflictTimes.size(); i++)
			{
				if(conflictTimes[i].second != 1)
				{
					cout<<"the conflict time should equal to 1,but not"<<endl;
					exit(0);
				}
			}
			//deletet the cluster if it has only one loop member
			int clusterNum = 1;
			for(std::set<int>::iterator its = confClusterSet.begin(); its != confClusterSet.end(); its++)
			{
				if(clusterizer.getClusterByID(*its).size() > clusterNum)
					// clusterizer.getClusterByID(*it).size();
				{
					cout<<"cluster "<<*its<<" has "<<clusterizer.getClusterByID(*its).size() <<" loops"<<endl;
					clusterNum = clusterizer.getClusterByID(*its).size();
				}

			}
			if (clusterNum == 1)
			{
				cout<<"**************only one loop cluster conflict exist"<<endl;
				return 0;
			}

			std::vector<int> deleted;
			for(std::set<int>::iterator its = confClusterSet.begin(); its != confClusterSet.end(); its++)
			{
				cout<<"see if run this"<<endl;
				cout<<"cluster "<<*its<<" has "<<clusterizer.getClusterByID(*its).size() <<" loops"<<endl;
		
				if(clusterizer.getClusterByID(*its).size() == 1)
				{
					clustersDeletedInConflictCheck.push_back(*its);
					goodSet.erase(*its);
					deleted.push_back(*its);
					cout<<"2 delete one ele clsuter "<<*its<<endl;
				}
				for(int i = confPairVector.size()-1; i >= 0; i--)
				{
					if(confPairVector[i].first == *its or confPairVector[i].second == *its)
						confPairVector.erase(confPairVector.begin()+ i);
				}
			}
			for(int i = 0; i < deleted.size(); i++)
			{
				confClusterSet.erase(deleted[i]);
				for(int i = confPairVector.size()-1; i >= 0; i--)
				{
					if(confPairVector[i].first == deleted[i] or confPairVector[i].second == deleted[i])
						confPairVector.erase(confPairVector.begin()+ i);
				}
				confClusterSet.clear();
				for(int i = confPairVector.size()-1; i >= 0; i--)
				{
					confClusterSet.insert(confPairVector[i].first);
					confClusterSet.insert(confPairVector[i].second);
				}
			}
			cout<<"clusterNum is "<<clusterNum<<endl;
			cout<<"end of conflict check program*****"<<endl;
			if(confClusterSet.size() == 0)
				return 1;
			else
				return 0;
		}
	}

	bool write(const char* filename)
	{
		IntPairSet correctLinks ;
		gatherLinks(_goodSet, correctLinks);
		cout<<"size of remained links: "<<correctLinks.size()<<" size of clusters in goodSET:"<<_goodSet.size()<<endl;
		gWrapper->write(filename,correctLinks);

		return true;
	}

	bool write_resolved_result(const char* filename, const char* resultFileName)
	{
		IntPairSet correctLinks ;
		gWrapper->saveOptimizedResult(resultFileName);
		gatherLinks(_goodSet, correctLinks);
		gWrapper->write(filename,correctLinks);
		return true;
	}

	bool removeIncorrectLoops()
	{
		size_t count = clusterizer.clusterCount();

		IntSet rejected;

		for(size_t i = 0 ; i<count ; i++)
		{
			if(_goodSet.find(i)==_goodSet.end())
			{
				rejected.insert(i);
			}
		}
		rejected.insert(ID_IGNORE);

		IntPairSet rejectedLinks;

		gatherLinks(rejected,rejectedLinks);
		gWrapper->removeEdges(rejectedLinks);

		// Clear up the clusters from Clusterizer

		for( IntSet::const_iterator it = rejected.begin(), end = rejected.end(); it!=end; it++)
		{
			clusterizer.deleteCluster(*it);
		}

		return true;
	}

};




#endif /* RRR_HPP_ */
