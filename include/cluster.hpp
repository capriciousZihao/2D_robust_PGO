// RRR - Robust Loop Closing over Time
// Copyright (C) 2014 Y.Latif, C.Cadena, J.Neira


#ifndef CLUSTER_HPP_
#define CLUSTER_HPP_

#include <iostream>
#include <sstream>
#include <fstream> 
#include <string.h>
#include <vector>
#include <list>
#include <cstdlib>

#include <math.h>
#include <unistd.h>
#include <algorithm>  
#include <assert.h>

#include "g2o/core/eigen_types.h"
#include "g2o/types/slam2d/se2.h"
#include "g2o/types/slam2d/edge_se2.h"
#include "g2o_interface_se2_zihao.hpp"


#include <Eigen/Dense>
#include "utils.hpp"

using namespace g2o;
using namespace Eigen; 
using namespace std;
using namespace utils;
// using namespace g2o::tutorial;

G2O_USE_TYPE_GROUP(slam2d); 

struct cluster
{
	int startLow, startHigh;
	int endLow, endHigh, nearestLCID=-1;
	int size, id_cluster;
	std::vector<std::array<double,6>> positionserial;
	std::vector<std::vector<std::pair<int, double> > > consistentGroup;
	std::vector<std::pair<int,int > > uncert_pair_group;
	std::vector< std::pair<double, std::pair<std::pair<int, int>, int> > > chiStatisWhenAdd2Cluster;
	// std::vector< std::pair<int,int>> nodesInCluster;

	cluster(): startLow(-1), startHigh(-1), endLow(-1), endHigh(-1), size(0) {}
	cluster(int start, int end, int idofcluster) : startLow(start), startHigh(start), endLow(end), endHigh(end), size(1), id_cluster(idofcluster){}

};

class Clusterizer
{
	typedef IntPairIDMap 			LoopToClusterIDMap;//loop is a pair of Vertex ID, this map is a int pair,LC, to a int, cluster ID 
	typedef IDintPairSetMap 		ClusterIDtoLoopsMap;
	double fiveNineBelief = 28, nineNineBelief = 44.8413, twelveBelief = 58.9413, upperLimitChi = 77.4, fiveupperLimitChi = 387, tenUpperLimit = 774,
		twoNineBelief = 11.345, threeNineBelief = 16.27, fourNineBelief = 21.11, Belief95 = 7.8147, twoupperLimitChi = 154.8;
	
		// double snrThres = 10;

	std::vector<std::vector< std::vector<int> > >  node_sequence;
	std::vector<std::vector<std::vector<std::pair<g2o::SE2, Matrix3d> > > > transSequence_whole;
	

	std::map<std::pair<int, int>, std::array<double,9>> LC_Inf;	
	std::vector<double> distance;

	std::vector<std::vector<int> > merge_vector;
	std::vector<int> killed_cluster;
	
	int xx = 0;
	int wholeLoopN = 0;

public:
	std::vector<std::pair<int ,int>> good_loops_vector, bad_loops_vector;
	std::set<int> mixed_clusters;
	std::vector<std::array<double,4>> VertexInf;
	std::set<std::pair<int, int > > bad_loops_set_self_check;
	double snrThres = 1, index_dispersion; //3.16; 10
	bool futile_two_segments;
	std::vector<std::pair<int, int> > node_distance_cons;
	std::vector<std::pair<int,int> > nodeDisVector;//for each test loop , find the nearest node dis
	std::vector<std::pair<int, int> > futilePairVector;
	Matrix3d displayCov, Cov_two_odo;
	ClusterIDtoLoopsMap	clusterIDtoLoopsMap;
	std::map<std::pair<int,int> , std::vector<int> >  LoopMultiClusterIDMap;
	int conflictNotInValid;
	std::vector<int>	vectorNewSplit,	vectorNewSplitCloser;
	LoopToClusterIDMap  loopToClusterIDMap, conflictClusterPairDistanceMap;
	std::set<std::pair<int, int> > conflictPairSet;
	std::map<std::pair<int, int>, std::set<std::pair<int, int> > > conflit_clsuter_pair_map_cause_loops;
	std::vector<std::pair<int,std::vector<std::pair<int, int> >  > > conflict_cause;
	std::vector<std::array<double,11>> OdoInf;
	std::vector<double> syntheticOdoAccumulate_dis;
	std::vector<std::pair<g2o::SE2, Matrix3d> > syntheticOdoAccumulate_trans_varMatrix;

	std::vector<std::pair<g2o::SE2, Matrix3d> > syntheticOdo10_trans_varMatrix, 
		syntheticOdo100_trans_varMatrix, syntheticOdo1000_trans_varMatrix, syntheticOdo10000_trans_varMatrix;
	std::vector<double > syntheticOdo10_dis, syntheticOdo100_dis, syntheticOdo1000_dis, syntheticOdo10000_dis;

	std::vector<std::pair<g2o::SE2, Matrix3d> > syntheticOdo10_trans_varMatrix_adjoint, 
		syntheticOdo100_trans_varMatrix_adjoint, syntheticOdo1000_trans_varMatrix_adjoint, syntheticOdo10000_trans_varMatrix_adjoint;
	std::vector<double > syntheticOdo10_dis_adjoint, syntheticOdo100_dis_adjoint, syntheticOdo1000_dis_adjoint, syntheticOdo10000_dis_adjoint;

	int NumOfRealLoops = -1;
	string     nameofclusterfile;
	std::vector<cluster> _clustersFound;
	std::set<std::pair<int,int> > set4GTL;
	std::vector<std::pair<int, std::vector<int> > > all_conflict_cluster;
	std::vector<std::pair<int, std::vector<int> > > all_uncertain_cluster;	
	std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	LP_Trans_Covar_Map;	
	find_element find_ele;
	std::vector<std::pair<int, int> > clusterSizeVector;

	static bool cmp(  const std::pair<int, int> p,  const std::pair<int, int> q)
	{
		return p.second > q.second;
	}

	static bool cmp_reverse(  const std::pair<int, int> p,  const std::pair<int, int> q)
	{
		return p.second < q.second;
	}
	// Assumes that in the vector the values are passed as (start_1,end_1), (start_2,end_2), ...
	int getClusterID(const IntPair& loop)//input a pair, lC, output the ID of the cluster to which the pair belongs
	{
		LoopToClusterIDMap::iterator it;
		it = loopToClusterIDMap.find(loop);
		if(it == loopToClusterIDMap.end())
		{
			cout<<"the loop "<<loop.first<<" "<<loop.second <<" can not be found in map"<<endl;
			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);	
		}
		// else
		// 	cout<<"check if get clUSter ID IS right "<<(*it).first.first<<" "<<(*it).first.second<<" id: "<<(*it).second <<endl;
		return loopToClusterIDMap[loop];
	}
	void getClusterID_multiExist(const IntPair& loop, std::vector<int> & mulCluster)//input a pair, lC, output the ID of the cluster to which the pair belongs
	{
		std::map<std::pair<int,int> , std::vector<int> >::iterator it;
		it = LoopMultiClusterIDMap.find(loop);
		if(it == LoopMultiClusterIDMap.end())
		{
			// cout<<"the loop can not be found in map"<<endl;
			// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			mulCluster.clear();
			return ;	
		}
		// else
		// 	cout<<"check if get clUSter ID IS right "<<(*it).first.first<<" "<<(*it).first.second<<" id: "<<(*it).second <<endl;
		mulCluster = LoopMultiClusterIDMap[loop];
		if(mulCluster.size() == 0)
		{
			cout<<"get a loop "<<loop.first<<" "<<loop.second <<" 's id , exist in multi cluster, but return an empty vector"<<endl;
			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		return ;
	}

	void updateConflictPairSet(std::vector<int> & conflict, int finalClusterID, std::vector<int> & uncertain, bool show)
	{
		std::pair<int, int> a, b;
		if(conflict.size() > 0)
		{
			if(finalClusterID >= 0)
			{
				for(int i = 0; i < conflict.size(); i++)
				{
					a.first = finalClusterID;
					a.second = conflict[i];
					b.first = conflict[i];
					b.second = finalClusterID;
					if(conflictPairSet.find(a) == conflictPairSet.end())
					{
						conflictPairSet.insert(a);
						conflictPairSet.insert(b);
						// cout<<"add conflict pair "<<finalClusterID<<" "<<conflict[i]<<endl;
					}

				}	
			}
			else
			{
				cout<<"clsuter id is  negative "<<endl;
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0) ;				
			}
		}
		else
		{
			// return;
		}
		if(uncertain.size()> 0)
		{
			if(finalClusterID >= 0)
			{
				for(int i = 0; i < uncertain.size(); i++)
				{
					a.first = finalClusterID;
					a.second = uncertain[i];
					b.first = uncertain[i];
					b.second = finalClusterID;
					if(conflictPairSet.find(a) != conflictPairSet.end())
					{
						conflictPairSet.erase(a);
						conflictPairSet.erase(b);
						// cout<<"delete from unvertain pair "<<finalClusterID<<" "<<uncertain[i]<<endl;
					}

				}	
			}
			else
			{
				cout<<"clsuter id is  negative "<<endl;
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0) ;				
			}
		}
		if(show)
		{
			if(conflict.size() > 0)
			{
				cout<<"conflict cluster pair:"<<endl;
				for(int i = 0; i < conflict.size(); i++)
				{
					cout<<finalClusterID<<" "<<conflict[i]<<endl;
				}
			}
			else
				cout<<"no conflict info"<<endl;
	
			if(uncertain.size() > 0)
			{
				cout<<"umcertain cluster pair:"<<endl;
				for(int i = 0; i < uncertain.size(); i++)
				{
					cout<<finalClusterID<<" "<<uncertain[i]<<endl;
				}
			}
			else
				cout<<"no uncertain info"<<endl;
				
		}
	}

	void updateConflictCause(int finalClusterID, 
		std::vector<std::pair<int,std::vector<std::pair<int, int> >  > >  & conflict_cause)
	{
		if(conflict_cause.size() > 0)
		{
			if(finalClusterID >= 0)
			{
				std::pair<int, int> a, b;
				for(int i = 0; i < conflict_cause.size(); i++)
				{
					a.first = finalClusterID;
					a.second = conflict_cause[i].first;
					b.first = conflict_cause[i].first;
					b.second = finalClusterID;
					for(int j = 0; j <conflict_cause[i].second.size();j++)
					{
						conflit_clsuter_pair_map_cause_loops[a].insert(conflict_cause[i].second[j]);
						conflit_clsuter_pair_map_cause_loops[b].insert(conflict_cause[i].second[j]);
					}

				}	
			}
			else
			{
				cout<<"clsuter id is  negative "<<endl;
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0) ;				
			}
		}
		else
		{
			return;
		}
	}

	void deleteClusterID_multiExist(const IntPair& loop, int clustgerID)//input a pair, lC, output the ID of the cluster to which the pair belongs
	{
		std::map<std::pair<int,int> , std::vector<int> >::iterator it;
		it = LoopMultiClusterIDMap.find(loop);
		if(it == LoopMultiClusterIDMap.end())
		{
			cout<<"the loop can not be found in multi Cluster Exist"<<endl;
			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0) ;	
		}
		int num = LoopMultiClusterIDMap[loop].size();
		if(num == 1 )
		{
			auto it = LoopMultiClusterIDMap.find(loop);
			LoopMultiClusterIDMap.erase(it);
		}
		else
			LoopMultiClusterIDMap[loop].erase(std::remove(LoopMultiClusterIDMap[loop].begin(), LoopMultiClusterIDMap[loop].end(), clustgerID), LoopMultiClusterIDMap[loop].end());
		if(num == 1)
		{
			if(LoopMultiClusterIDMap.find(loop) != LoopMultiClusterIDMap.end())
			{
				cout<<"multiexist delete problem"<<endl;
				cout<<"size:"<<LoopMultiClusterIDMap[loop].size()<<endl;
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0) ;					
			}
		}
		if(num == 2)
        {

        }
		return ;
	}

	bool get_cluster_node_scequence(std::vector<std::vector< std::vector<int> > > & node_sequence, const std::vector<cluster> & clustersFound)
	{

		std::vector<int> nodes_dout, nodes_din;
		std::vector<std::vector<int> > element;
		bool retu = 0;
		node_sequence.clear();
		for(int i = 0; i < clustersFound.size(); ++i)
		{
			nodes_dout.clear();
			nodes_din.clear();	
			element.clear();	
			for(int j = 0; j < clustersFound[i].positionserial.size(); ++j)
			{
				nodes_dout.push_back(clustersFound[i].positionserial[j][0]);
				nodes_din.push_back(clustersFound[i].positionserial[j][3]);
			}
			element.push_back(nodes_dout);
			element.push_back(nodes_din);	
			node_sequence.push_back(element);
		}
		// exit(0);
		retu = 1;
		return retu;
	}

	void setClusterID(const IntPair& loopClosure, int ID)
	{
		if(loopToClusterIDMap.find(loopClosure)!=loopToClusterIDMap.end())
		{
			int oldID = loopToClusterIDMap[loopClosure];

			clusterIDtoLoopsMap[oldID].erase(loopClosure);

			clusterIDtoLoopsMap[ID].insert(loopClosure);
			loopToClusterIDMap[loopClosure] = ID;
			int clsuterSIZE = _clustersFound.size();
			if(ID > clsuterSIZE)
			{
				if(ID >_clustersFound.size() )
				{
					cout<<"cluster ID too big: "<<ID<<" while the whole clsuter number is "<<_clustersFound.size()<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);	
				}
			}
		}
	}

	std::array<double,6> get_LC_Pos(int start, int end)
	{
		std::array<double,6> fullLoopInfo;
		bool findStartVertexInMatrix=false,findEndVertexInMatrix=false;
		//get the position of start and end point std::vector<std::array<double,4>> VertexInf
		for(std::vector<std::array<double,4>>::const_iterator itVertex = VertexInf.begin(), 
			lendVertex = VertexInf.end();itVertex!=lendVertex;itVertex++)
		{
			if(int((*itVertex)[0])==start)
			{
                   fullLoopInfo[1]=(*itVertex)[1]; fullLoopInfo[2]=(*itVertex)[2]; //fullLoopInfo[2]=(*itVertex)[3];
				findStartVertexInMatrix=true;
			}
			if(int((*itVertex)[0])==end)
			{
				fullLoopInfo[4]=(*itVertex)[1]; fullLoopInfo[5]=(*itVertex)[2];// fullLoopInfo[2]=(*itVertex)[3];
				findEndVertexInMatrix=true;
			}
			if(findStartVertexInMatrix and findEndVertexInMatrix)
				break;
		}
		if(!(findStartVertexInMatrix and findEndVertexInMatrix))
			cout<<"can't find the position of the poses of the loop"<<endl;
		if(!findStartVertexInMatrix)
		{
			cout<<"start vertex id: "<<start<<endl;
		}
		if(!findEndVertexInMatrix)
		{
			cout<<"end vertex id: "<<end<<endl;
		}

		fullLoopInfo[0]=start;//  fullLoopInfo[1]=startPosition[0];  fullLoopInfo[2]=startPosition[1];
		fullLoopInfo[3]=end;  //  fullLoopInfo[4]=endPosition[0];    fullLoopInfo[5]=endPosition[1];
		return fullLoopInfo;
	}

	void prepare_to_signle_loop_pair_check_invalid_exist(std::pair<int, int> & loop_node_pair, IntPairSet::const_iterator & LP_nodes,
										   std::array<std::pair<g2o::SE2, Matrix3d>, 4> & FullInfo,
										   std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	& LP_Trans_Covar_Map,
										   std::array<std::array<double, 5>, 5> & length, bool & fultileBit1,  bool & fultileBit2,
										   std::pair<int, int> & validVertexPoint)
	{
		//input variable: std::pair<int, int>  and full info of four segments.
		int start, end;
		bool fultileBit_;
		g2o::SE2  Trans1, Trans2, midTrans, midTrans1, Trans1_;
		std::pair<g2o::SE2, Matrix3d > result_synthetic_odo, result_synthetic_odo2;
		Matrix3d  cov1, cov2, midCov, cov1_;
		bool fultileBit = 0;
		std::array<double ,5> length_;

		length[0][0] = 0; length[0][1] = 0; length[0][2] = 0;
		length[1][0] = 0; length[1][1] = 0; length[1][2] = 0;
		length[2][0] = 0; length[2][1] = 0; length[2][2] = 0;
		length[3][0] = 0; length[3][1] = 0; length[3][2] = 0;

		start = loop_node_pair.first;
		end   = (*LP_nodes).first;
		if(validVertexPoint.first > 0 and (abs(start - end) > validVertexPoint.first))
		{
			fultileBit1 = 1;
		}
		else
		{
			if(start > end)
			{
				// synthesize_odo_edges( end,  start,  OdoInf, Trans1,  cov1, length[0], fultileBit);///
				accelerated_synthetic_odo( end, start, OdoInf, 
					result_synthetic_odo,  length[0], 
					fultileBit1);

				Trans1 = result_synthetic_odo.first.inverse();
				cov1   = result_synthetic_odo.second;
			}
			else
			{
				// synthesize_odo_edges( start,  end,  OdoInf, Trans1,  cov1, length[0], fultileBit);
				accelerated_synthetic_odo(start, end, OdoInf, 
					result_synthetic_odo,  length[0], 
					fultileBit1);

				Trans1 = result_synthetic_odo.first;
				cov1   = result_synthetic_odo.second;
			}	
		}
		//fultileBit == 1 means the variance is larger then the distance, so we think there is no need to do next calculation
		if(fultileBit1 == 1)
		{
			if(validVertexPoint.first < 0)
				validVertexPoint.first = abs(start - end);
			else if (abs(start - end) < validVertexPoint.first)
			{
				validVertexPoint.first = abs(start - end);
			}
			return;
		}

		//synthesize the end odometry edges
		start = (*LP_nodes).second ;
		end   = loop_node_pair.second;

		if(validVertexPoint.second > 0 and (abs(start - end) > validVertexPoint.second))
		{
			fultileBit2 = 1;
		}
		else
		{
			if(start > end)
			{
				accelerated_synthetic_odo( end, start, OdoInf, 
					result_synthetic_odo,  length[2], 
					fultileBit2);	

				Trans2 = result_synthetic_odo.first.inverse();
				cov2   = result_synthetic_odo.second;
			}
			else
			{
				// synthesize_odo_edges( start,  end,  OdoInf, Trans2,  cov2, length[2], fultileBit);
				accelerated_synthetic_odo(start, end, OdoInf, 
					result_synthetic_odo,  length[2], 
					fultileBit2);

				Trans2 = result_synthetic_odo.first;
				cov2   = result_synthetic_odo.second;
			}
		}
		if(fultileBit2 == 1)
		{
			if(validVertexPoint.second < 0)
				validVertexPoint.second = abs(start - end);
			else if (abs(start - end) < validVertexPoint.second)
			{
				validVertexPoint.second = abs(start - end);
			}
			return;
		}

		//assign value of 4 parts which from one cycle, the sequence is from part1{lastloop.start > test_loop.start} to part2{test_loop}
		// to part3 {test_loop.second > last_loop.second} to part4{the inverse of last_loop} 
		FullInfo[0] = std::pair<g2o::SE2, Matrix3d> (Trans1, cov1);
		FullInfo[1] = LP_Trans_Covar_Map[(*LP_nodes)];
		FullInfo[2] = std::pair<g2o::SE2, Matrix3d> (Trans2, cov2);
		midTrans    = LP_Trans_Covar_Map[loop_node_pair].first;
		midTrans    = midTrans.inverse();
		midCov      = LP_Trans_Covar_Map[loop_node_pair].second;
		FullInfo[3] = std::pair<g2o::SE2, Matrix3d> (midTrans, midCov);

		length[1][0] = sqrt(FullInfo[1].first[0]*FullInfo[1].first[0]+FullInfo[1].first[1]*FullInfo[1].first[1]);
		length[1][1] = abs(FullInfo[1].first[0]);
		length[1][2] = abs(FullInfo[1].first[1]);
		length[1][3] = FullInfo[1].second(0,0);
		length[1][4] = FullInfo[1].second(1,1);

		length[3][0] = sqrt(FullInfo[3].first[0]*FullInfo[3].first[0]+FullInfo[3].first[1]*FullInfo[3].first[1]);
		length[3][1] = abs(FullInfo[3].first[0]);
		length[3][2] = abs(FullInfo[3].first[1]);
		length[3][3] = FullInfo[3].second(0,0);
		length[3][4] = FullInfo[3].second(1,1);	

		length[4][0] = length[0][0]+length[1][0]+length[2][0]+length[3][0];
		length[4][1] = length[0][1]+length[1][1]+length[2][1]+length[3][1];
		length[4][2] = length[0][2]+length[1][2]+length[2][2]+length[3][2];
	}

	void prepare_to_signle_loop_pair_check_adjoint(std::pair<int, int> & loop_node_pair, const std::pair<int, int> & LP_nodes,
										   std::array<std::pair<g2o::SE2, Matrix3d>, 4> & FullInfo,
										   std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	& LP_Trans_Covar_Map,
										   std::array<std::array<double, 5>, 5> & length, bool & fultileBit1,  bool & fultileBit2)
	{
		//input variable: std::pair<int, int>  and full info of four segments.
		int start, end;
		bool fultileBit_;
		g2o::SE2  Trans1, Trans2, midTrans, midTrans1, Trans1_;
		std::pair<g2o::SE2, Matrix3d > result_synthetic_odo, result_synthetic_odo2;
		Matrix3d  cov1, cov2, midCov, cov1_;
		bool fultileBit = 0;
		std::array<double ,5> length_;
		length[0][0] = 0; length[0][1] = 0; length[0][2] = 0;
		length[1][0] = 0; length[1][1] = 0; length[1][2] = 0;
		length[2][0] = 0; length[2][1] = 0; length[2][2] = 0;
		length[3][0] = 0; length[3][1] = 0; length[3][2] = 0;

		start = loop_node_pair.first;
		end   = (LP_nodes).first;
		if(start > end)
		{
			// synthesize_odo_edges( end,  start,  OdoInf, Trans1,  cov1, length[0], fultileBit);///
			accelerated_synthetic_odo_adjoint( end, start, OdoInf, 
				result_synthetic_odo,  length[0], 
				fultileBit1);
			Trans1 = result_synthetic_odo.first.inverse();
			cov1   = result_synthetic_odo.second;
		}
		else
		{
			// synthesize_odo_edges( start,  end,  OdoInf, Trans1,  cov1, length[0], fultileBit);
			accelerated_synthetic_odo_adjoint(start, end, OdoInf, 
				result_synthetic_odo,  length[0], 
				fultileBit1);
			Trans1 = result_synthetic_odo.first;
			cov1   = result_synthetic_odo.second;
		}
		//synthesize the end odometry edges
		start = (LP_nodes).second ;
		end   = loop_node_pair.second;

		if(start > end)
		{
			// synthesize_odo_edges( end,  start,  OdoInf, Trans2,  cov2, length[2], fultileBit);
			// Trans2 = Trans2.inverse();
			accelerated_synthetic_odo_adjoint( end, start, OdoInf, 
				result_synthetic_odo,  length[2], 
				fultileBit2);
			Trans2 = result_synthetic_odo.first.inverse();
			cov2   = result_synthetic_odo.second;
		}
		else
		{
			// synthesize_odo_edges( start,  end,  OdoInf, Trans2,  cov2, length[2], fultileBit);
			accelerated_synthetic_odo_adjoint(start, end, OdoInf, 
				result_synthetic_odo,  length[2], 
				fultileBit2);
			Trans2 = result_synthetic_odo.first;
			cov2   = result_synthetic_odo.second;
		}
		//assign value of 4 parts which from one cycle, the sequence is from part1{lastloop.start > test_loop.start} to part2{test_loop}
		// to part3 {test_loop.second > last_loop.second} to part4{the inverse of last_loop} 
		FullInfo[0] = std::pair<g2o::SE2, Matrix3d> (Trans1, cov1);
		FullInfo[1] = LP_Trans_Covar_Map[(LP_nodes)];
		FullInfo[2] = std::pair<g2o::SE2, Matrix3d> (Trans2, cov2);
		midTrans    = LP_Trans_Covar_Map[loop_node_pair].first;
		midTrans    = midTrans.inverse();
		midCov      = LP_Trans_Covar_Map[loop_node_pair].second;
		FullInfo[3] = std::pair<g2o::SE2, Matrix3d> (midTrans, midCov);

		length[1][0] = sqrt(FullInfo[1].first[0]*FullInfo[1].first[0]+FullInfo[1].first[1]*FullInfo[1].first[1]);
		length[1][1] = abs(FullInfo[1].first[0]);
		length[1][2] = abs(FullInfo[1].first[1]);
		length[1][3] = FullInfo[1].second(0,0);
		length[1][4] = FullInfo[1].second(1,1);

		length[3][0] = sqrt(FullInfo[3].first[0]*FullInfo[3].first[0]+FullInfo[3].first[1]*FullInfo[3].first[1]);
		length[3][1] = abs(FullInfo[3].first[0]);
		length[3][2] = abs(FullInfo[3].first[1]);
		length[3][3] = FullInfo[3].second(0,0);
		length[3][4] = FullInfo[3].second(1,1);	

		length[4][0] = length[0][0]+length[1][0]+length[2][0]+length[3][0];
		length[4][1] = length[0][1]+length[1][1]+length[2][1]+length[3][1];
		length[4][2] = length[0][2]+length[1][2]+length[2][2]+length[3][2];

		// cout<<"cov1:"<<endl<<cov1<<endl<<"cov2: "<<cov2<<endl;
	}
	void prepare_to_signle_loop_pair_check(std::pair<int, int> & loop_node_pair, const std::pair<int, int> & LP_nodes,
										   std::array<std::pair<g2o::SE2, Matrix3d>, 4> & FullInfo,
										   std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	& LP_Trans_Covar_Map,
										   std::array<std::array<double, 5>, 5> & length, bool & fultileBit1,  bool & fultileBit2)
	{
		//input variable: std::pair<int, int>  and full info of four segments.
		int start, end;
		bool fultileBit_;
		g2o::SE2  Trans1, Trans2, midTrans, midTrans1, Trans1_;
		std::pair<g2o::SE2, Matrix3d > result_synthetic_odo, result_synthetic_odo2;
		Matrix3d  cov1, cov2, midCov, cov1_;
		bool fultileBit = 0;
		std::array<double ,5> length_;
		length[0][0] = 0; length[0][1] = 0; length[0][2] = 0;
		length[1][0] = 0; length[1][1] = 0; length[1][2] = 0;
		length[2][0] = 0; length[2][1] = 0; length[2][2] = 0;
		length[3][0] = 0; length[3][1] = 0; length[3][2] = 0;

		start = loop_node_pair.first;
		end   = (LP_nodes).first;
		LOG(INFO)<<"start ,end = "<<start<<" "<<end;
		if(start > end)
		{
			// synthesize_odo_edges( end,  start,  OdoInf, Trans1,  cov1, length[0], fultileBit);///
			accelerated_synthetic_odo( end, start, OdoInf, 
				result_synthetic_odo,  length[0], 
				fultileBit1);
			Trans1 = result_synthetic_odo.first.inverse();
			cov1   = result_synthetic_odo.second;
		}
		else
		{
			// synthesize_odo_edges( start,  end,  OdoInf, Trans1,  cov1, length[0], fultileBit);
			accelerated_synthetic_odo(start, end, OdoInf, 
				result_synthetic_odo,  length[0], 
				fultileBit1);
			Trans1 = result_synthetic_odo.first;
			cov1   = result_synthetic_odo.second;
		}
		//synthesize the end odometry edges
		LOG(INFO)<<"trans is "<<Trans1[0]<<" "<<Trans1[1]<<" "<<Trans1[2];
		// cout<<"start ,end = "<<end<<" "<<start<<endl;
		// cout<<"trans is "<<Trans1.inverse()[0]<<" "<<Trans1.inverse()[1]<<" "<<Trans1.inverse()[2]<<endl;
		LOG(INFO)<<" ";

		start = (LP_nodes).second ;
		end   = loop_node_pair.second;
		LOG(INFO)<<"start ,end = "<<start<<" "<<end;
		if(start > end)
		{
			// synthesize_odo_edges( end,  start,  OdoInf, Trans2,  cov2, length[2], fultileBit);
			// Trans2 = Trans2.inverse();
			accelerated_synthetic_odo( end, start, OdoInf, 
				result_synthetic_odo,  length[2], 
				fultileBit2);
			Trans2 = result_synthetic_odo.first.inverse();
			cov2   = result_synthetic_odo.second;
		}
		else
		{
			// synthesize_odo_edges( start,  end,  OdoInf, Trans2,  cov2, length[2], fultileBit);
			accelerated_synthetic_odo(start, end, OdoInf, 
				result_synthetic_odo,  length[2], 
				fultileBit2);
			Trans2 = result_synthetic_odo.first;
			cov2   = result_synthetic_odo.second;
		}
		LOG(INFO)<<"trans is "<<Trans2[0]<<" "<<Trans2[1]<<" "<<Trans2[2]<<endl;
		// cout<<"start ,end = "<<end<<" "<<start<<endl;
		// cout<<"trans is "<<Trans2.inverse()[0]<<" "<<Trans2.inverse()[1]<<" "<<Trans2.inverse()[2]<<endl;
		LOG(INFO)<<" "<<endl;
		//assign value of 4 parts which from one cycle, the sequence is from part1{lastloop.start > test_loop.start} to part2{test_loop}
		// to part3 {test_loop.second > last_loop.second} to part4{the inverse of last_loop} 
		FullInfo[0] = std::pair<g2o::SE2, Matrix3d> (Trans1, cov1);
		FullInfo[1] = LP_Trans_Covar_Map[(LP_nodes)];
		LOG(INFO)<<"loop "<<LP_nodes.first<<" "<<LP_nodes.second<<" : ";
		LOG(INFO)<<LP_Trans_Covar_Map[(LP_nodes)].first[0]<<" "<<LP_Trans_Covar_Map[(LP_nodes)].first[1]<<" "<<LP_Trans_Covar_Map[(LP_nodes)].first[2]<<" ";
		// cout<<LP_Trans_Covar_Map[(LP_nodes)].first.toIsometry().matrix()<<endl;
		// cout<<(LP_Trans_Covar_Map[(LP_nodes)].first.inverse()).toIsometry().matrix()<<endl;

		FullInfo[2] = std::pair<g2o::SE2, Matrix3d> (Trans2, cov2);
		midTrans    = LP_Trans_Covar_Map[loop_node_pair].first;
		LOG(INFO)<<"loop "<<loop_node_pair.first <<" "<<loop_node_pair.second<<" : "<<LP_Trans_Covar_Map[(LP_nodes)].first[0]<<
			" "<<LP_Trans_Covar_Map[(LP_nodes)].first[1]<<" "<<LP_Trans_Covar_Map[(LP_nodes)].first[2];
		LOG(INFO)<<midTrans.toIsometry().matrix()<<endl;
		midTrans    = midTrans.inverse();
		// cout<<midTrans.toIsometry().matrix()<<endl;
		
		midCov      = LP_Trans_Covar_Map[loop_node_pair].second;
		FullInfo[3] = std::pair<g2o::SE2, Matrix3d> (midTrans, midCov);

		length[1][0] = sqrt(FullInfo[1].first[0]*FullInfo[1].first[0]+FullInfo[1].first[1]*FullInfo[1].first[1]);
		length[1][1] = abs(FullInfo[1].first[0]);
		length[1][2] = abs(FullInfo[1].first[1]);
		length[1][3] = FullInfo[1].second(0,0);
		length[1][4] = FullInfo[1].second(1,1);

		length[3][0] = sqrt(FullInfo[3].first[0]*FullInfo[3].first[0]+FullInfo[3].first[1]*FullInfo[3].first[1]);
		length[3][1] = abs(FullInfo[3].first[0]);
		length[3][2] = abs(FullInfo[3].first[1]);
		length[3][3] = FullInfo[3].second(0,0);
		length[3][4] = FullInfo[3].second(1,1);	

		length[4][0] = length[0][0]+length[1][0]+length[2][0]+length[3][0];
		length[4][1] = length[0][1]+length[1][1]+length[2][1]+length[3][1];
		length[4][2] = length[0][2]+length[1][2]+length[2][2]+length[3][2];
	}


	void find_cons_cluster(IntPairSet::const_iterator & LP_nodes,  std::vector<cluster>& _clustersFound, std::vector<int> & cons_cluster_number,
		std::vector<int>  & conflict_cluster, std::vector<int>  & uncertain_cluster, std::map<std::pair<int, int>, 
		std::pair<g2o::SE2, Matrix3d> > & LP_Trans_Covar_Map, 
		std::vector<std::array<double,4>> & VertexInf, bool debug, std::vector<double> & chiStatis, bool & futileBit)
	{
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		std::array<std::array<double, 5>, 5>  length;
		int start, end;
		std::array<double,3> returndis;
		std::array<double,2> returnmid;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair;
		cons_cluster_number.clear();
		conflict_cluster.clear();
		uncertain_cluster.clear();
		chiStatis.clear();
		std::pair<bool, double> reV, reV_second_last;
		double  transX_residual, transY_residual, transA_residual;
		bool futileBit1, futileBit2;

		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
		for(int i=_clustersFound.size()-1; i >= 0; i--)
		{
			bool pass  = 0;
			int limit = 1;
			double bigChiErr = 0, passChi = 100.0;
			if(_clustersFound[i].positionserial.size() > 1)
				limit = 2;

			//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
			double transformDistance = 100000000000.0;
			loop_node_pair.first = _clustersFound[i].positionserial.back()[0];
			loop_node_pair.second = _clustersFound[i].positionserial.back()[3];
			prepare_to_signle_loop_pair_check(loop_node_pair, *LP_nodes, FullInfo, LP_Trans_Covar_Map, length, futileBit1, futileBit2);


			if((futileBit1 or futileBit2) != 1)
			{	
				reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
				transformDistance = reV.second;
				bigChiErr = transformDistance;
				// if((reV.first == 1) and (futileBit == 0))
				if(reV.first == 1)
				{
					pass  = 1;
					if(transformDistance < passChi)
						passChi = transformDistance;

				}

				if(_clustersFound[i].positionserial.size() > 2)
				{
 					int si = _clustersFound[i].positionserial.size();
					loop_node_pair.first  = _clustersFound[i].positionserial[si-3][0];
					loop_node_pair.second = _clustersFound[i].positionserial[si-3][3];
					prepare_to_signle_loop_pair_check(loop_node_pair, *LP_nodes, FullInfo, LP_Trans_Covar_Map, length, futileBit1, futileBit2);
					reV_second_last = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual,
						length[4][0], futileBit);
					if(reV_second_last.second < transformDistance)
						transformDistance = reV_second_last.second;
					if(reV_second_last.second > bigChiErr)
						bigChiErr = reV_second_last.second;
					// if((reV_second_last.first == 1) and (futileBit == 0))
					if(reV_second_last.first == 1)
					{
						if(bigChiErr < upperLimitChi)
						{
							pass  = 1;
							if(transformDistance < passChi)
								passChi = transformDistance;
						}
					}
				}
					cout<<" "<<endl;
					cout<<"fultileBit == 0"<<endl;
					cout<<"test      loop: "<<(*LP_nodes).first<<" "<<(*LP_nodes).second<<endl;
					cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< i<<endl;
					cout<<"transformDistance: "<<transformDistance<<endl;
					cout<<"length: "<<length[0][0]<<" "<<length[1][0]<<" "<<length[2][0]<<" "<<length[3][0]<<" "<<length[4][0]<<" "<<endl;
					cout<<"cov:"<<endl;
					cout<<displayCov<<endl;
			}
			else
			{
				reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
				transformDistance = reV.second;
				cout<<" "<<endl;
				cout<<"fultileBit == 1"<<endl;
				cout<<"test      loop: "<<(*LP_nodes).first<<" "<<(*LP_nodes).second<<endl;
				cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< i<<endl;
				cout<<"transformDistance: "<<transformDistance<<endl;
				cout<<"length: "<<length[0][0]<<" "<<length[1][0]<<" "<<length[2][0]<<" "<<length[3][0]<<" "<<length[4][0]<<" "<<endl;
							printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				cout<<"cov:"<<endl;
				cout<<displayCov<<endl;

				if(transformDistance > fiveupperLimitChi )//tenUpperLimit
				{
					cout<<"add to conflict set although the NSR is higher than threshold"<<endl;
				}
			}
			//if current test cluster variance too big add it to the uncertain set
			if( abs(transformDistance - 100000000000.0) <0.001)
			{
				// uncertain_cluster.push_back(i);
				// continue;
			}
			//if pass chi2 test add this cluster id to the consistent set
			if(pass == 1)
			{
				cons_cluster_number.push_back(i);
				chiStatis.push_back(passChi);			
			}
			//if fail to pass the chi2 test, handle it based on the transform distance // nineNineBelief twelveBelief
			else
			{
				if(transformDistance >  upperLimitChi )//tenUpperLimit
					conflict_cluster.push_back(i);
				// else if(transform_distance > 16.266)//%0.1
				// {
				// 	uncertain_cluster.insert(i);
				// }
				else 
				{
					// cons_cluster_number.push_back(i);
					// chiStatis.push_back(transform_distance);
					// reV.first = 1;
					// reV.second =transform_distance;
					// break;
					uncertain_cluster.push_back(i);
				}
			}

		}
	}
 	void find_cons_cluster_varying_belief_all_cluster(IntPairSet::const_iterator & LP_nodes,  std::vector<cluster>& _clustersFound, std::vector<int> & cons_cluster_number,
		std::vector<int>  & conflict_cluster, std::vector<int>  & uncertain_cluster, std::map<std::pair<int, int>, 
		std::pair<g2o::SE2, Matrix3d> > & LP_Trans_Covar_Map, 
		std::vector<std::array<double,4>> & VertexInf, bool debug, std::vector<double> & chiStatis, bool & futileBit, double & beli)
	{
		// std::vector<std::pair<int,std::vector<std::pair<int, int> >  > > 
		node_distance_cons.clear();
		std::vector<std::pair<int, int> > conflict_cause_pairs;
		std::pair<int,std::vector<std::pair<int, int> >  > conflict_cause_element;
		conflict_cause.clear();
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		std::array<std::array<double, 5>, 5>  length;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair;
		cons_cluster_number.clear();
		conflict_cluster.clear();
		uncertain_cluster.clear();
		chiStatis.clear();
		std::pair<bool, double> reV;
		double  transX_residual, transY_residual, transA_residual;
		bool futileBit1, futileBit2;
		std::vector<double> transDisWhole;
		std::vector<int> consisBitVector;

		// reV =  check_single_loop_odo(*LP_nodes, LP_Trans_Covar_Map, futileBit1, beli);
		// if(reV.second > upperLimitChi)
		// {
		// 	bad_loops_set_self_check.insert(*LP_nodes);
		// 	cout<<"bad loop cause the self dis is "<<reV.second<<endl;	
		// 	return;
		// }
					// if(reV.second < 80)
					
		// find the biggest set
		int biggestCluster = 0;
		for(int i=_clustersFound.size()-1; i >= 0; i--)
		{
			if(getClusterByID(i).size() > getClusterByID(biggestCluster).size() )
				biggestCluster = getClusterByID(i).size();
		}
		std:set<int> checkedClusterSet;
		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
		std::pair<int, int> Pair;
		for(int i=_clustersFound.size()-1; i >= 0; i--)
		{
			conflict_cause_pairs.clear();
			bool pass  = 0;
			double bigChiErr = 0, passChi = 100.0, smallChiErr = -1.0;
			double transformDistance;
			nodeDisVector.clear();
			// cout<<"test clsuter "<<i<<endl;

			for(int wholeDis = 0; wholeDis < _clustersFound[i].positionserial.size(); wholeDis++)
			{
				Pair.first = wholeDis;
				Pair.second = abs((*LP_nodes).first - _clustersFound[i].positionserial[wholeDis][0]) + 
							  abs((*LP_nodes).second - _clustersFound[i].positionserial[wholeDis][3]);

				// cout<<" "<<endl<<"first: "<<wholeDis<<" second: "<<Pair.second<<endl;
				nodeDisVector.push_back(Pair);
			}	
			// cout<<"nodeDisVector size "<<nodeDisVector.size()<<endl;	
			std::sort (nodeDisVector.begin(), nodeDisVector.end(), cmp_reverse);
			//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
			int loopIDNearest, findTime;
			if(nodeDisVector.size() > 10)
				findTime = 10;
			else
				findTime = nodeDisVector.size() ;

			int node_dis_cons = -1;
			for(int wholeDis = 0; wholeDis < findTime; wholeDis++)
			{
				//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
				transformDistance = 100000000000.0;
				loopIDNearest = nodeDisVector[wholeDis].first;
				// cout<<" "<<endl;
				// cout<<"node dis "<<nodeDisVector[wholeDis].second<<" in cluster "<<i<<endl;

				//if the node dis is bigger than it self, do self check, if self check fail, put into bad one and exit, if 
				//pass self check then do next staff.
				// cout<<" "<<endl<<"first: "<<nodeDisVector[wholeDis].first<<" second: "<<nodeDisVector[wholeDis].second<<endl;
				if(nodeDisVector[wholeDis].second > abs((*LP_nodes).first - (*LP_nodes).second))
				{
					reV =  check_single_loop_odo(*LP_nodes, LP_Trans_Covar_Map, futileBit1, beli);
					if(reV.second < 80)
						cout<<"the self dis is "<<reV.second<<endl;	
					// if((futileBit1) != 1)
					// {
						// if((reV.first <= upperLimitChi))
					if((reV.first <= 7.8))
						{
							// disVector.push_back(reV.second);
							// if((*LP_nodes).first == 79)
							// 	cin.get();
						}		
						else
						{
							// if(reV.second > utils::chi2_continuous(3,0.95))
							// {
								conflict_cause_element.second = conflict_cause_pairs;
								conflict_cause_element.first = i;
								conflict_cause.push_back(conflict_cause_element);
								conflict_cluster.push_back(i);
								bad_loops_set_self_check.insert(*LP_nodes);
								return;
											
						}	

						if((futileBit1) == 1)
						{
							break;
						}				
					// }
					// else
					// 	break;
				}

				loop_node_pair.first = _clustersFound[i].positionserial[loopIDNearest][0];
				loop_node_pair.second = _clustersFound[i].positionserial[loopIDNearest][3];

				prepare_to_signle_loop_pair_check(loop_node_pair, *LP_nodes, FullInfo, LP_Trans_Covar_Map, 
					length, futileBit1, futileBit2);



				if((futileBit1 or futileBit2) != 1)
				{	
					reV = check_single_loop_inter_varying_belief(FullInfo, covX, covY, displayCov, transX_residual, 
						transY_residual, transA_residual, length[4][0], futileBit, beli);
					transformDistance = reV.second;
					transDisWhole.push_back(transformDistance);

					// cout<<(*LP_nodes).first<<" "<<(*LP_nodes).second<<" to fixed loop "
					// 	<<loop_node_pair.first<<" "<<loop_node_pair.second<< " dis: "<<transformDistance<<endl;

					if(futileBit == 0)
					{
						if(bigChiErr < transformDistance)
							bigChiErr = transformDistance;
						if(smallChiErr < 0)
							smallChiErr = transformDistance;	
						else if(smallChiErr > transformDistance )
							smallChiErr = transformDistance;						
						if(reV.first == 1)			
						{
							pass  = 1;
							if(node_dis_cons == -1)
								node_dis_cons = nodeDisVector[wholeDis].second;


							consisBitVector.push_back(1);
							if(passChi > transformDistance)
								passChi = transformDistance;
						}
						else
						{
							if(transformDistance > upperLimitChi)
							{
								conflict_cause_pairs.push_back(loop_node_pair);
							}
							consisBitVector.push_back(0);
						}
					}
					else
					{
						cout<<"futile inner"<<endl;
						consisBitVector.push_back(0);
						break;
					}
				}
				else
				{
					cout<<"futile"<<endl;
					consisBitVector.push_back(0);
					transDisWhole.push_back(-2);
					break;
				}

				// if(i == 25 and 4085 == ( *LP_nodes).second)
				// 	cout<<"futile1 futile2 futile dis: "<<futileBit1<<" "
				// 		<<futileBit2<<" "<<futileBit<<" "<<transformDistance<<endl;

				if(transformDistance == 0)
				{
					cout<<"distance is zero"<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
			}


			//if pass chi2 test add this cluster id to the consistent set
			// if(pass == 1 and bigChiErr< upperLimitChi )//fourNineBelief
			if(pass == 1 and bigChiErr< utils::chi2_continuous(3,0.999999995) ) //9999999
			{
				cons_cluster_number.push_back(i);
				chiStatis.push_back(passChi);		
				node_distance_cons.push_back(std::pair<int, int>(i, node_dis_cons))	;
				// cout<<"add to consistent with biggest error "<<bigChiErr<<endl;
			}
			//if fail to pass the chi2 test, handle it based on the transform distance // nineNineBelief twelveBelief
			else
			{
				// cout<<"smallChiErr: "<<smallChiErr<<"  "<<endl;
				// if(bigChiErr >  upperLimitChi and pass!= 1)//tenUpperLimit
				if(smallChiErr > upperLimitChi )
				{
					conflict_cause_element.second = conflict_cause_pairs;
					conflict_cause_element.first = i;
					conflict_cause.push_back(conflict_cause_element);
					conflict_cluster.push_back(i);
				}
				// else if(transform_distance > 16.266)//%0.1
				// {
				// 	uncertain_cluster.insert(i);
				// }
				else if(smallChiErr < utils::chi2_continuous(3,0.95) and smallChiErr > 0)
				{
					// cons_cluster_number.push_back(i);
					// chiStatis.push_back(transform_distance);
					// reV.first = 1;
					// reV.second =transform_distance;
					// break;
					uncertain_cluster.push_back(i);
				}
			}

		}
	} 
 	void find_cons_cluster_varying_belief(IntPairSet::const_iterator & LP_nodes,  std::vector<cluster>& _clustersFound, std::vector<int> & cons_cluster_number,
		std::vector<int>  & conflict_cluster, std::vector<int>  & uncertain_cluster, std::map<std::pair<int, int>, 
		std::pair<g2o::SE2, Matrix3d> > & LP_Trans_Covar_Map, 
		std::vector<std::array<double,4>> & VertexInf, bool debug, std::vector<double> & chiStatis, bool & futileBit, double & beli)
	{
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		std::array<std::array<double, 5>, 5>  length;
		int start, end;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair;
		cons_cluster_number.clear();
		conflict_cluster.clear();
		uncertain_cluster.clear();
		chiStatis.clear();
		std::pair<bool, double> reV, reV_second_last;
		double  transX_residual, transY_residual, transA_residual;
		bool futileBit1, futileBit2;

		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
		for(int i=_clustersFound.size()-1; i >= 0; i--)
		{
			bool pass  = 0;
			int limit = 1;
			double bigChiErr = 0, passChi = 100.0;
			if(_clustersFound[i].positionserial.size() > 1)
				limit = 2;

			//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
			double transformDistance = 100000000000.0;
			loop_node_pair.first = _clustersFound[i].positionserial.back()[0];
			loop_node_pair.second = _clustersFound[i].positionserial.back()[3];
			prepare_to_signle_loop_pair_check(loop_node_pair, *LP_nodes, FullInfo, LP_Trans_Covar_Map, 
				length, futileBit1, futileBit2);

			if((futileBit1 or futileBit2) != 1)
			{	
				reV = check_single_loop_inter_varying_belief(FullInfo, covX, covY, displayCov, transX_residual, 
					transY_residual, transA_residual, length[4][0], futileBit, beli);
				transformDistance = reV.second;
				bigChiErr = transformDistance;
				// if((reV.first == 1) and (futileBit == 0))
				if(reV.first == 1)
				{
					pass  = 1;
					if(transformDistance < passChi)
						passChi = transformDistance;

				}

				if(_clustersFound[i].positionserial.size() > 2)
				{
 					int si = _clustersFound[i].positionserial.size();
					loop_node_pair.first  = _clustersFound[i].positionserial[si-3][0];
					loop_node_pair.second = _clustersFound[i].positionserial[si-3][3];
					prepare_to_signle_loop_pair_check(loop_node_pair, *LP_nodes, FullInfo, LP_Trans_Covar_Map, length, futileBit1, futileBit2);
					if((futileBit1 or futileBit2) != 1)
					{
						reV_second_last = check_single_loop_inter_varying_belief(FullInfo, covX, covY, displayCov, transX_residual, 
							transY_residual, transA_residual, length[4][0], futileBit, beli);
						if(reV_second_last.second < transformDistance)
							transformDistance = reV_second_last.second;
						if(reV_second_last.second > bigChiErr)
							bigChiErr = reV_second_last.second;
						// if((reV_second_last.first == 1) and (futileBit == 0))
						if(reV_second_last.first == 1)
						{
							if(bigChiErr < upperLimitChi)
							{
								pass  = 1;
								if(transformDistance < passChi)
									passChi = transformDistance;
							}
						}
					}

				}
					cout<<" "<<endl;
					cout<<"fultileBit == 0"<<endl;
					cout<<"test      loop: "<<(*LP_nodes).first<<" "<<(*LP_nodes).second<<endl;
					cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< i<<endl;
					cout<<"transformDistance: "<<transformDistance<<endl;
					cout<<"length: "<<length[0][0]<<" "<<length[1][0]<<" "<<length[2][0]<<" "<<length[3][0]<<" "<<length[4][0]<<" "<<endl;
					cout<<"cov:"<<endl;
					cout<<displayCov<<endl;
			}
			else
			{
				reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
				transformDistance = reV.second;
				cout<<" "<<endl;
				cout<<"fultileBit == 1"<<endl;
				cout<<"test      loop: "<<(*LP_nodes).first<<" "<<(*LP_nodes).second<<endl;
				cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< i<<endl;
				cout<<"transformDistance: "<<transformDistance<<endl;
				cout<<"length: "<<length[0][0]<<" "<<length[1][0]<<" "<<length[2][0]<<" "<<length[3][0]<<" "<<length[4][0]<<" "<<endl;
							printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				cout<<"cov:"<<endl;
				cout<<displayCov<<endl;

				if(transformDistance > fiveupperLimitChi )//tenUpperLimit
				{
					cout<<"add to conflict set although the NSR is higher than threshold"<<endl;
				}
			}
			//if current test cluster variance too big add it to the uncertain set
			if( abs(transformDistance - 100000000000.0) <0.001)
			{
				// uncertain_cluster.push_back(i);
				// continue;
			}
			//if pass chi2 test add this cluster id to the consistent set
			if(pass == 1)
			{
				cons_cluster_number.push_back(i);
				chiStatis.push_back(passChi);			
			}
			//if fail to pass the chi2 test, handle it based on the transform distance // nineNineBelief twelveBelief
			else
			{
				if(transformDistance >  upperLimitChi )//tenUpperLimit
					conflict_cluster.push_back(i);
				// else if(transform_distance > 16.266)//%0.1
				// {
				// 	uncertain_cluster.insert(i);
				// }
				else 
				{
					// cons_cluster_number.push_back(i);
					// chiStatis.push_back(transform_distance);
					// reV.first = 1;
					// reV.second =transform_distance;
					// break;
					uncertain_cluster.push_back(i);
				}
			}

		}
	} 

 // 	void find_cons_cluster_varying_belief_No_Conflict(IntPairSet::const_iterator & LP_nodes,  std::vector<cluster>& _clustersFound, 
 // 		std::vector<int> & cons_cluster_number,
	// 	std::vector<int>  & conflict_cluster, std::vector<int>  & partConsis_cluster, std::map<std::pair<int, int>, 
	// 	std::pair<g2o::SE2, Matrix3d> > & LP_Trans_Covar_Map, 
	// 	std::vector<std::array<double,4>> & VertexInf, bool debug, std::vector<double> & chiStatis, bool & futileBit, double & beli,
	// 	std::vector<std::vector<int>> &partConsisVV)
	// {
	// 	double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
	// 	std::array<std::array<double, 5>, 5>  length;
	// 	int start, end;
	// 	std::array<double,3> returndis;
	// 	std::array<double,2> returnmid;
	// 	std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
	// 	g2o::SE2  Trans1, Trans2, midTrans;
	// 	Matrix3d  cov1, cov2, midCov;
	// 	std::pair<int, int> loop_node_pair;
	// 	cons_cluster_number.clear();
	// 	conflict_cluster.clear();
	// 	partConsisVV.clear();
	// 	partConsis_cluster.clear();
	// 	chiStatis.clear();
	// 	std::pair<bool, double> reV, reV_second_last;
	// 	double  transX_residual, transY_residual, transA_residual;
	// 	bool futileBit1, futileBit2;
	// 	std::vector<int> transDisWhole, consisBitVector;


	// 	//iterate the elements in clusters, if find one consistent cluster, jump out loop.
	// 	for(int i=_clustersFound.size()-1; i >= 0; i--)
	// 	{
	// 		transDisWhole.clear();
	// 		consisBitVector.clear();
	// 		//two menbers should be equal
	// 		if(_clustersFound[i].positionserial.size() != getClusterByID(i).size())
	// 		{
	// 			for(int deNotE = 0; deNotE < _clustersFound[i].positionserial.size(); deNotE++)
	// 			{
	// 				cout<<"position serial loop "<<_clustersFound[i].positionserial[deNotE][0]<<" "<<_clustersFound[i].positionserial[deNotE][3]<<endl;
	// 			}
	// 			std::set<std::pair<int,int > > loopSet = getClusterByID(i);
	// 			for(auto it = loopSet.begin(); it != loopSet.end(); it++)
	// 			{
	// 				cout<<"getby id serial loop "<<(*it).first<<" "<<(*it).second<<endl;
	// 			}
	// 			cout<<_clustersFound[i].positionserial.size()<<" "<<getClusterByID(i).size()<<endl;
	// 			cout<<"loop set not equal to positionserial, in cluster "<<i<<endl;
	// 			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
	// 			exit(0);
	// 		}
	// 		bool pass  = 0;
	// 		int limit = 1;
	// 		double bigChiErr = 0, passChi = 100.0, minChi =200.0;
	// 		if(_clustersFound[i].positionserial.size() > 1)
	// 			limit = 2;

	// 		//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop

	// 		for(int wholeDis = 0; wholeDis < _clustersFound[i].positionserial.size(); wholeDis++)
	// 		{
	// 			//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
	// 			double transformDistance = 100000000000.0;
	// 			loop_node_pair.first = _clustersFound[i].positionserial[wholeDis][0];
	// 			loop_node_pair.second = _clustersFound[i].positionserial[wholeDis][3];
	// 			prepare_to_signle_loop_pair_check(loop_node_pair, *LP_nodes, FullInfo, LP_Trans_Covar_Map, 
	// 				length, futileBit1, futileBit2);

	// 			if((futileBit1 or futileBit2) != 1)
	// 			{	
	// 				reV = check_single_loop_inter_varying_belief(FullInfo, covX, covY, displayCov, transX_residual, 
	// 					transY_residual, transA_residual, length[4][0], futileBit, beli);
	// 				transformDistance = reV.second;
	// 				transDisWhole.push_back(transformDistance);
	// 				if(bigChiErr < transformDistance)
	// 					bigChiErr = transformDistance;

	// 				// if((reV.first == 1) and (futileBit == 0))
	// 				// cout<<" "<<endl;
	// 				// cout<<"fultileBit == 0"<<endl;
					
	// 				if(reV.first == 1)			
	// 				{
	// 					pass  = 1;
	// 					consisBitVector.push_back(1);
	// 					// if(minChi > transformDistance)
	// 					// 	minChi = transformDistance;
	// 				}
	// 				else
	// 					consisBitVector.push_back(0);
	// 			}
	// 			else
	// 			{
	// 				consisBitVector.push_back(0);
	// 				transDisWhole.push_back(-2);
	// 			}

	// 			if(transformDistance == 0)
	// 			{
	// 				cout<<"distance is zero"<<endl;
	// 				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
	// 				exit(0);
	// 			}
	// 		}

	// 		if(consisBitVector.size() != _clustersFound[i].positionserial.size())
	// 		{
	// 			cout<<"loop set not equal to positionserial"<<endl;
	// 			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
	// 			exit(0);
	// 		}

	// 		int allConsNum = 0;
	// 		// cout<<"consBit: ";
	// 		for(int SUM = 0; SUM < consisBitVector.size(); SUM++)
	// 		{
	// 			// cout<<consisBitVector[SUM]<<" ";
	// 			allConsNum = allConsNum + consisBitVector[SUM];
	// 		}
	// 		// cout<<endl;
	// 		// allConsNum = 0;
	// 		// cout<<"each allConsNum: ";
	// 		// cout<<allConsNum<<" ";
	// 		// for(int SUM = 0; SUM < consisBitVector.size(); SUM++)
	// 		// {
				
	// 		// 	allConsNum = allConsNum + consisBitVector[SUM];
	// 		// 	cout<<allConsNum<<" ";
	// 		// }
	// 		// cout<<endl;
	// 		int sum_of_elems = 0;
	// 		std::for_each(consisBitVector.begin(), consisBitVector.end(), [&] (int n) {
	// 		    sum_of_elems += n;
	// 		});
	// 		// cout<<"sum bit using for each: "<<sum_of_elems<<endl;
	// 		if(allConsNum != sum_of_elems)
	// 		{
	// 			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
	// 			exit(0);
	// 		}

	// 		if(allConsNum == _clustersFound[i].positionserial.size())//all member loop is consistent with the test loop
	// 		{
	// 			// cout<<"sum all 1 bit: "<<allConsNum<<endl;
	// 			// cout<<"its position serial number : "<<_clustersFound[i].positionserial.size()<<endl;
	// 			cout<<"add cluster "<<i<<" to cons"<<endl;
	// 			cons_cluster_number.push_back(i);
	// 		}
	// 		else if(allConsNum > 0 and bigChiErr > threeNineBelief)
	// 		{
	// 			//the max error > %99.9
	// 			cout<<"add cluster "<<i<<" to uncertain"<<endl;
	// 			partConsis_cluster.push_back(i);
	// 			partConsisVV.push_back(consisBitVector);
	// 		}
	// 		else
	// 			cout<<"no cons or uncertain"<<endl;

	// 	}
	// 	cin.get();
	// } 
 
	void calTrandformDis_PairSet_Pair(std::set<std::pair<int,int> > & goodSET, std::pair<int, int> doubtLoop, std::vector<double> & disVector)
	{

	    std::array<std::pair<g2o::SE2, Matrix3d>, 4>  FullInfo;
	    std::array<std::array<double, 5>, 5>  length;
	    std::pair<int,int> eleInSet;
	    std::pair<bool, double> reV;
	    bool  futileBit1, futileBit2, futileBit;
	    double  deltaX1, deltaY1, deltaX2, deltaY2, covX, covY , transX_residual, transY_residual, transA_residual;
	    disVector.clear();

		std::set<std::pair<int,int> >::iterator it = goodSET.begin(), end = goodSET.end();
		for(; it !=end; it++ )
		{
			eleInSet.first = (*it).first;
			eleInSet.second = (*it).second;
			prepare_to_signle_loop_pair_check(doubtLoop, eleInSet, FullInfo, LP_Trans_Covar_Map, length, futileBit1, futileBit2);
			reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
			disVector.push_back(reV.second);
		}
		cout<<"good set size : "<<goodSET.size();
		if(disVector.size() != goodSET.size())
		{
			cout<<"two vector size is not equal"<<endl;
			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
	}
	double calTransformDis_loop_loop( std::pair<int, int> Loop1, std::pair<int, int> Loop2)
	{

	    std::array<std::pair<g2o::SE2, Matrix3d>, 4>  FullInfo;
	    std::array<std::array<double, 5>, 5>  length;
	    std::pair<int,int> eleInSet;
	    std::pair<bool, double> reV;
	    bool  futileBit1, futileBit2, futileBit;
	    double  deltaX1, deltaY1, deltaX2, deltaY2, covX, covY , transX_residual, transY_residual, transA_residual;

		prepare_to_signle_loop_pair_check(Loop1, Loop2, FullInfo, LP_Trans_Covar_Map, length, futileBit1, futileBit2);
		reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
		if(futileBit == 1)
		{
			cout<<"futileBit = 1"<<"  "<<"transform dis is "<<reV.second <<endl;
			return reV.second*(-1);
		}

		return reV.second;
	}

	void find_cons_cluster_second(IntPairSet::const_iterator & LP_nodes,  std::vector<cluster>& _clustersFound, std::vector<int> & cons_cluster_number,
		std::vector<int>  & conflict_cluster, std::vector<int>  & uncertain_cluster, std::map<std::pair<int, int>, 
		std::pair<g2o::SE2, Matrix3d> > & LP_Trans_Covar_Map, 
		std::vector<std::array<double,4>> & VertexInf, bool debug, std::vector<double> & chiStatis, bool & futileBit,
		std::vector<std::pair<int,std::vector<int> > > & split, std::vector<std::pair<int,std::vector<double> > > & splitDis,
		std::map<int, double> & conflictClusterDistanceMap) // std::vector<std::vector<int> > & disForConsistent
	{
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		std::array<std::array<double, 5>, 5>  length;
		int start, end;
		std::array<double,3> returndis;
		std::array<double,2> returnmid;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair;
		cons_cluster_number.clear();
		conflict_cluster.clear();
		uncertain_cluster.clear();
		chiStatis.clear();
		std::pair<bool, double> reV, reV_second_last;
		double  transX_residual, transY_residual, transA_residual, mind, maxd, maxdValid;
		bool futileBit1, futileBit2;
		std::vector<double> transDisWhole, distanceWholeValid;
		split.clear();
		splitDis.clear();
		std::pair<int,std::vector<int> > eleToSplit;
		std::pair<int, std::vector<double> > eleSplitDis;
		conflictClusterDistanceMap.clear();
		std::pair<int, int> validVertexPoint(-1,-1);

		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
		for(int i=_clustersFound.size()-1; i >= 0; i--)
		{
			bool pass  = 0;
			int limit = 1;

			transDisWhole.clear();
			distanceWholeValid.clear();
			eleToSplit.second.clear();
			eleSplitDis.second.clear();
			maxdValid = 0;

			for(int wholeDis = 0; wholeDis < _clustersFound[i].positionserial.size(); wholeDis++)
			{
				//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
				double transformDistance = 100000000000.0;
				loop_node_pair.first = _clustersFound[i].positionserial[wholeDis][0];
				loop_node_pair.second = _clustersFound[i].positionserial[wholeDis][3];
				prepare_to_signle_loop_pair_check(loop_node_pair, *LP_nodes, FullInfo, LP_Trans_Covar_Map, 
					length, futileBit1, futileBit2);

				if((futileBit1 or futileBit2) != 1)
				{	
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transDisWhole.push_back(transformDistance);

					// if((reV.first == 1) and (futileBit == 0))
					cout<<" "<<endl;
					cout<<"fultileBit == 0"<<endl;

					if(futileBit == 0)
                    {
                        distanceWholeValid.push_back(transformDistance);
                        if(transformDistance > maxdValid)
                            maxdValid = transformDistance;

                        if(reV.first == 1)
                        {
                            pass  = 1;
                        }
                    }
				}
				else
				{
					cout<<" "<<endl;
					cout<<"fultileBit == 1"<<endl;
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transDisWhole.push_back(transformDistance);
					if(transformDistance > upperLimitChi )//tenUpperLimit
						cout<<"add to conflict set although the NSR is higher than threshold"<<endl;				
				}
				if(transformDistance == 0)
				{
					cout<<"distance is zero"<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
			}
			//find the max and min distacne in whole distace vector and their corresponding index
			mind = transDisWhole[0]; maxd = 0;
			int minSerial = 0, maxSerial = 0;
			for(int findMax = 0; findMax < transDisWhole.size(); findMax++)
			{
				if(transDisWhole[findMax] > maxd)
				{
					maxd = transDisWhole[findMax];
					maxSerial = findMax;
				}
				if(transDisWhole[findMax] < transDisWhole[minSerial])
				{
					mind = transDisWhole[findMax];
					minSerial = findMax;
				}
			}

			//display the distacnes stored in whole vector and that in valid vector
			cout<<"cluster "<<i<<endl;
			cout<<"mind: "<<mind<<" "<<"maxd: "<<maxd<<endl;	
			cout<<"distanceWholeValid  serial: ";
			for(int toPrint = 0; toPrint < distanceWholeValid.size(); toPrint++)
			{
				cout<<distanceWholeValid[toPrint]<<" ";
				if(distanceWholeValid[toPrint] > maxdValid)
					maxdValid = distanceWholeValid[toPrint];
			}
			cout<<endl;
			cout<<"distance serial: ";
			for(int toPrint = 0; toPrint < transDisWhole.size(); toPrint++)
			{
				cout<<transDisWhole[toPrint]<<" ";
			}
			cout<<endl;

			//check if the number in wholeDistance equal to element in the cluster 
			if(transDisWhole.size() != _clustersFound[i].positionserial.size())
			{
				cout<<"the size of element in cluster and whole valid distance is not equal"<<endl;
				cout<<"whole element: "<<_clustersFound[i].positionserial.size()<<
					" distance whole: "<<transDisWhole.size()<<endl;
				std::cin.get();
			}

			//check if the condition is met, which valid number should not bigger than whole elemnt
			if(transDisWhole.size() < distanceWholeValid.size())
			{
				cout<<"the size of valid distance and whole valid distance is not equal"<<endl;
				cout<<"the size of whole should be "<<_clustersFound[i].positionserial.size()<<endl;
				cout<<"valid: "<<distanceWholeValid.size()<<" whole: "<<transDisWhole.size()<<endl;
				std::cin.get();
			}

			//if pass chi2 test add this cluster id to the consistent set
			bool  futileBit69;
			std::vector<double> transformDistanceCluster, transformDistanceCluster2;
			if(pass == 1)
			{
				/*if pass but there is also a loop exist in the cluster that has a larger distacne,
					we need to propose it*/
				// if(maxdValid > Belief95)//fiveNineBelief
				if(maxdValid > threeNineBelief)
				{
					cout<<"clusterID: "<<i<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					// exit(0);

					cal_whole_cluster_distance_calculation(maxSerial,  _clustersFound,  
					LP_Trans_Covar_Map,  0, futileBit69, transformDistanceCluster, i);
					// double reBiggestDistance = *(std::max_element(transformDistanceCluster.begin(), transformDistanceCluster.end()));
					double reBiggestDistance = transformDistanceCluster[minSerial];
					double maxMan = *(std::max_element(transformDistanceCluster.begin(), transformDistanceCluster.end()));
					cout<<"distance to the biggest error element: ";
					for(int toPrint = 0; toPrint < transformDistanceCluster.size(); toPrint++)
					{
						cout<<transformDistanceCluster[toPrint]<<" ";
						if(transformDistanceCluster[toPrint] > Belief95)
						{
							// cout<<"this distacne "<<transformDistanceCluster[toPrint]<<" is bigger than %95 limit"<<endl;
							// std::cin.get();
						}
					}
					cout<<endl;
					//it means the test loop shold not be added to the cluster and the cluster should not be splited
					if(reBiggestDistance < mind and reBiggestDistance < Belief95 and maxMan < Belief95)
					{
						cout<<"reBiggestDistance is "<<reBiggestDistance<<endl;
						uncertain_cluster.push_back(i);
						cout<<"bu zuo split because origianl distance is smaller"<<endl;
						// std::cin.get();
						continue;
					}

					// //the biggest distance of the potential element with 
					// double BDP = *(std::max_element(transformDistanceCluster.begin(), transformDistanceCluster.end()));
					// //the second biggest distance 
					// int find_SecondBiggest(std::vector<double > & vector, double & obj)

					// transformDistanceCluster.clear();
					cal_whole_cluster_distance_calculation(minSerial,  _clustersFound,  
					LP_Trans_Covar_Map,  0, futileBit69, transformDistanceCluster2, i);
					cout<<"distance to the smallest error element: ";
					for(int toPrint = 0; toPrint < transformDistanceCluster2.size(); toPrint++)
					{
						cout<<transformDistanceCluster2[toPrint]<<" ";
					}
					cout<<endl;	


					// if the distane of min < dis bewteen min max,split
					if(transformDistanceCluster[minSerial] > mind)
					{
						//do split 
						eleToSplit.first = i;
						eleSplitDis.first = i;

						if(transformDistanceCluster.size() != transDisWhole.size())
						{
							printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
							exit(0);
						}
						for(int toPrint = 0; toPrint < transformDistanceCluster.size(); toPrint++)
						{
							//transformDistanceCluster stores the distacne bewteen maxSerial and other loops in the clsuter 
							if((transformDistanceCluster[toPrint] < transformDistanceCluster2[toPrint]) and 
							   (transformDistanceCluster[toPrint] < transDisWhole[toPrint]))
							{
								// //the sistacne of loops that need to be splicted out is bigger than %95
								// if((transformDistanceCluster[toPrint] > Belief95))
								// {
								// 	cout<<"transformDistanceCluster[toPrint]: "<<transformDistanceCluster[toPrint]<<endl;
								// 	printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
								// 	// exit(0);	
								// }
								// else
								{
									eleSplitDis.second.push_back(transformDistanceCluster[toPrint]);
									eleToSplit.second.push_back(toPrint);					
								}
							}
						}	

						split.push_back(eleToSplit);
						splitDis.push_back(eleSplitDis);
						cout<<"min distance ok but max diatance bigger than %95"<<endl;
						// std::cin.get();
						//add to split info, split						
						cons_cluster_number.push_back(i);
						chiStatis.push_back(mind);
						cout<<"add to split and consistent"<<endl;
						printf("This is in %s on line %d\n",  __FILE__, __LINE__);

					}
					else
					{
						uncertain_cluster.push_back(i);
						cout<<"dis biggest is smaller than that to test loop, so no split"<<endl;
						cout<<"add to unvertain"<<endl;
						printf("This is in %s on line %d\n",  __FILE__, __LINE__);
						
					}

				}
				else
				{
					/*	if pass is 1, it means this cluster is consistent to the test loop, despite there may some loops 
						also exitst in this cluster have a larger distance with the test loop
					*/
					cout<<"add to consistent"<<endl;
					// printf("This is in %s on line %d\n",  __FILE__, __LINE__);
					cons_cluster_number.push_back(i);
					chiStatis.push_back(mind);		
				}
	
			}
			//if fail to pass the chi2 test, handle it based on the transform distance // nineNineBelief
			else
			{
				// if(transformDistance >  fiveupperLimitChi )//tenUpperLimit  nineNineBelief   upperLimitChi twelveBelief
				// if(maxd > upperLimitChi)
					
				if(maxdValid > upperLimitChi)
				{
					if(maxdValid < upperLimitChi and distanceWholeValid.size() != 0 and maxdValid < 15)
					{
						conflictNotInValid++;
						cout<<"maxd("<<maxd<<")"<<" > "<<upperLimitChi<<" but maxdValid("<<maxdValid<<")" <<" < "<<upperLimitChi<<endl;
						cout<<"current num in conflictNotInValid is "<<conflictNotInValid<<endl;
						// conflictPairOnlyInWhole.push_back()
						std::cin.get();
					}
					cout<<"add to conflict"<<endl;
					// printf("This is in %s on line %d\n",  __FILE__, __LINE__);					
					conflict_cluster.push_back(i);
					conflictClusterDistanceMap[i] = maxd;
				}
				else 
				{
					cout<<"add to unvertain"<<endl;
					// printf("This is in %s on line %d\n",  __FILE__, __LINE__);
					uncertain_cluster.push_back(i);
				}
			}
			// std::cin.get();

		}
		// wholeLoopN++;
	}


	void cal_whole_cluster_distance_calculation(
		int serial,  std::vector<cluster>& _clustersFound,  
		std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> > & LP_Trans_Covar_Map,  
		bool debug, bool & futileBit, 
		std::vector<double> & transformDistanceCluster, int & toTestClusterID)
	{
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		std::array<std::array<double, 5>, 5>  length;
		int start, end, eleSize = _clustersFound[toTestClusterID].positionserial.size();
		std::array<double,3> returndis;
		std::array<double,2> returnmid;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair, loop_node_pair_base; 
		std::pair<bool, double> reV, reV_second_last;
		double  transX_residual, transY_residual, transA_residual;
		bool futileBit1, futileBit2;
		
		transformDistanceCluster.clear();

		loop_node_pair_base.first = _clustersFound[toTestClusterID].positionserial[serial][0];
		loop_node_pair_base.second = _clustersFound[toTestClusterID].positionserial[serial][3];


		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
		// for(int i=_clustersFound.size()-1; i >= 0; i--)

		{
			bool pass  = 0;
			for(int wholeClusterTransformDistance = 0; 
				wholeClusterTransformDistance < eleSize; wholeClusterTransformDistance++)
			{
				//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
				double transformDistance = 100000000000.0;
				loop_node_pair.first = _clustersFound[toTestClusterID].positionserial[wholeClusterTransformDistance][0];
				loop_node_pair.second = _clustersFound[toTestClusterID].positionserial[wholeClusterTransformDistance][3];
				prepare_to_signle_loop_pair_check(loop_node_pair, loop_node_pair_base, FullInfo, LP_Trans_Covar_Map, length, futileBit1, futileBit2);

				if((futileBit1 or futileBit2) != 1)
				{	
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transformDistanceCluster.push_back(transformDistance);
					if(reV.first == 1 and (futileBit == 0))
					// if(reV.first == 1 )
						pass  = 1;
					if(debug)
					{
						cout<<" "<<endl;
						cout<<"fultileBit == 0"<<endl;
						cout<<"test      loop: "<<(loop_node_pair_base).first<<" "<<(loop_node_pair_base).second<<endl;
						cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< toTestClusterID<<endl;
						cout<<"transformDistance: "<<transformDistance<<endl;
						cout<<"cov:"<<endl;
						cout<<displayCov<<endl;
					}
				}
				else
				{
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transformDistanceCluster.push_back(transformDistance);
					if(debug)
					{
						cout<<" "<<endl;
						cout<<"fultileBit == 1"<<endl;
						cout<<"test      loop: "<<(loop_node_pair_base).first<<" "<<(loop_node_pair_base).second<<endl;
						cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< toTestClusterID<<endl;
						cout<<"transformDistance: "<<transformDistance<<endl;
						cout<<"cov:"<<endl;
						cout<<displayCov<<endl;
					}
				}

			}
		}
	}
	void cal_whole_cluster_distance_calculation_4intra(std::pair<int,int > loopToTest,    bool debug, 
		std::vector<double> & transformDistanceCluster, int & toTestClusterID, std::pair<int,int > & confPair, double & disMax)
	{
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		std::array<std::array<double, 5>, 5>  length;
		int start, end, eleSize = _clustersFound[toTestClusterID].positionserial.size();
		std::array<double,3> returndis;
		std::array<double,2> returnmid;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair, loop_node_pair_base; 
		std::pair<bool, double> reV, reV_second_last;
		double  transX_residual, transY_residual, transA_residual;
		bool futileBit1, futileBit2, futileBit;
		
		transformDistanceCluster.clear();

		loop_node_pair_base = loopToTest;
		cout<<"test      loop: "<<(loop_node_pair_base).first<<" "<<(loop_node_pair_base).second<<endl;

		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
		// for(int i=_clustersFound.size()-1; i >= 0; i--)
		bool pass  = 0;
		disMax = 0.1;
		for(int wholeClusterTransformDistance = 0; 
			wholeClusterTransformDistance < eleSize; wholeClusterTransformDistance++)
		{
			//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
			double transformDistance = 100000000000.0;
			loop_node_pair.first = _clustersFound[toTestClusterID].positionserial[wholeClusterTransformDistance][0];
			loop_node_pair.second = _clustersFound[toTestClusterID].positionserial[wholeClusterTransformDistance][3];
			prepare_to_signle_loop_pair_check(loop_node_pair, loop_node_pair_base, FullInfo, LP_Trans_Covar_Map, length, futileBit1, futileBit2);

			if((futileBit1 or futileBit2) != 1)
			{	
				reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
				transformDistance = reV.second;
				transformDistanceCluster.push_back(transformDistance);
				if(transformDistance > disMax)
				{
					confPair = loop_node_pair;
					disMax = transformDistance;
				}

				if(reV.first == 1 and (futileBit == 0))
				// if(reV.first == 1 )
					pass  = 1;
				if(debug)
				{
					cout<<" "<<endl;
					if(transformDistance > utils::chi2(3))
						cout<<"fultileBit == 0 ** ** ** ** ** ** ** ** ** ** ** ** ** ** **"<<endl;
					else
						cout<<"fultileBit == 0"<<endl;
					cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< toTestClusterID<<endl;
					cout<<"transformDistance: "<<transformDistance<<endl;
					// cout<<"cov:"<<endl;
					// cout<<displayCov<<endl;
				}
			}
			else
			{
				reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
				transformDistance = reV.second;
				if(transformDistance > disMax)
				{
					confPair = loop_node_pair;
					disMax = transformDistance;
				}

				transformDistanceCluster.push_back(transformDistance);
				if(debug)
				{
					cout<<" "<<endl;
					if(transformDistance > utils::chi2(3))
						cout<<"fultileBit == 1 ** ** ** ** ** ** ** ** ** ** ** ** ** ** **"<<endl;
					else
						cout<<"fultileBit == 1"<<endl;
					cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< toTestClusterID<<endl;
					cout<<"transformDistance: "<<transformDistance<<endl;
					// cout<<"cov:"<<endl;
					// cout<<displayCov<<endl;
				}
			}
		}
		if(debug == 1)
		{
			cout<<"trans serial: ";
			for(int i =0; i < transformDistanceCluster.size(); i++)
			{
				cout<<transformDistanceCluster[i]<<" ";
			}
			cout<<endl;
		}
	}
	bool judgeVarianceAndDistance(std::array<std::array<double, 5>, 5> & length)
	{
		bool accept = 1;
		for(int i =0; i < 4; i++)
		{
			for(int j =3; j<5; j++)
			{
				if(length[i][j] > length[i][0])
				{
					accept = 0;
					break;
				}	
			}
		}
		return accept;
	}

	bool judgeVarianceAndDistance_Small(std::array<double, 5> & length, Matrix3d & cov)
	{
		bool futileBitLocal = 0;

		// if(sqrt(sqrt(length[3]*length[3] + length[4]*length[4]) > 0.3*length[0]))
		// if(sqrt(length[3]*length[3] + length[4]*length[4]) > 13.25*length[0])
		// double SNR = (length[0])/sqrt(length[3]*length[3] + length[4]*length[4]);
		// double correlation = 
		if(length[0] < 0.1)
			return 0;
		// index_dispersion = sqrt(cov(0,0)*cov(0,0) + cov(1,1)*cov(1,1))/(length[0]);
		index_dispersion = sqrt(cov(0,0)*cov(0,0) + cov(0,1)*cov(0,1) + cov(1,0)*cov(1,0) + cov(1,1)*cov(1,1))/(length[0]);

		// cout<<"length[0]: "<<length[0]<<" length[3]: "<<length[3]<<" length[4]: "<<length[4]<<endl;
		// cout<<"sqrt(length[3]*length[3] + length[4]*length[4]): "<<sqrt(length[3]*length[3] + length[4]*length[4])<<endl;
		if(index_dispersion > snrThres)//3.16 corresponding to 5 dB
		// if(SNR < 1)
		// if(SNR < 10)
		{
			if(sqrt(length[3]*length[3] + length[4]*length[4]) > 0.1)
			{
				futileBitLocal = 1;
				// cout<<"futile bit  = 1"<<endl;
				// cin.get();
			}
			else
			{
				cout<<"something abmormal about futile check"<<endl;
				cout<<"length[0] length[3] length[4]: "<<length[0]<<" "<<length[3]<<" "<<length[4]<<endl;
				cout<<"cov:"<<endl;
				cout<<cov<<endl;
				cin.get();
			}				
		}
		else
		{
			// cout<<"futile bit = "
		}
		return futileBitLocal;
	}
	bool judgeVarianceAndDistance_Small(double length, Matrix3d & cov)
	{
		bool futileBitLocal = 0;

		// if(sqrt(sqrt(length[3]*length[3] + length[4]*length[4]) > 0.3*length[0]))
		// if(sqrt(length[3]*length[3] + length[4]*length[4]) > 13.25*length[0])

		// double SNR = (length[0])/sqrt(length[3]*length[3] + length[4]*length[4]);
		// double correlation = 
		if(length < 0.1)
			return 0;		
		// index_dispersion = sqrt(cov(0,0)*cov(0,0) + cov(1,1)*cov(1,1))/(length);
		index_dispersion = sqrt(cov(0,0)*cov(0,0) + cov(0,1)*cov(0,1) + cov(1,0)*cov(1,0) + cov(1,1)*cov(1,1))/(length);
		// cout<<"length[0]: "<<length[0]<<" length[3]: "<<length[3]<<" length[4]: "<<length[4]<<endl;
		// cout<<"sqrt(length[3]*length[3] + length[4]*length[4]): "<<sqrt(length[3]*length[3] + length[4]*length[4])<<endl;
		if(index_dispersion > snrThres)//3.16 corresponding to 5 dB
		// if(SNR < 1)
		// if(SNR < 10)
		{
			if(length > 0.1)
			{
				futileBitLocal = 1;
				// cout<<"futile bit  = 1"<<endl;
				// cin.get();
			}
			else
			{
				cout<<"something abnormal about futile check"<<endl;
				cin.get();
			}
		}
		else
		{
			// cout<<"futile bit = "
		}

		return futileBitLocal;
	}

	std::array<double,2> chi2_test(const std::vector<double> & dis_clu, const std::vector<double> & dis_clu_real)
	{
		std::array<double,2> p_value = {0, 0};
		double dis_backup;

		double sum = std::accumulate(std::begin(dis_clu), std::end(dis_clu), 0.0);  
		double mean =  sum / dis_clu.size(); //均值  
						  
		double accum  = 0.0;  
		std::for_each (std::begin(dis_clu), std::end(dis_clu), [&](const double d) {  
			accum  += (d-mean)*(d-mean);  
			});  
						  
		double chi_statis = (accum/mean); //卡方统计量
		p_value[0] = 1-chi2_p(dis_clu.size()-1, chi_statis);	
		if(dis_clu.size()>2)
		{
			if(dis_clu_real.size() != dis_clu.size()-1)
			{
				cout<<"dis_clu_real.size: "<<dis_clu_real.size()<<endl;
				cout<<"dis_clu.size: "<<dis_clu.size()<<endl;
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			
			sum = std::accumulate(std::begin(dis_clu_real), std::end(dis_clu_real), 0.0);  
			mean =  sum / dis_clu_real.size(); //均值  
			accum  = 0.0;  
			std::for_each (std::begin(dis_clu_real), std::end(dis_clu_real), [&](const double d) {  
			accum  += (d-mean)*(d-mean);  
			});  
						  
			 chi_statis = (accum/mean); //卡方统计量
			p_value[1] = 1-chi2_p(dis_clu_real.size()-1, chi_statis);	
		}
		else
		{
			p_value[1] = p_value[0];
		}
		// cout<<"p_value:"<<p_value[0]<<"  pre_p_value:"<<p_value[1]<<"  dof:"<<(dis_clu.size()-1)<<"  chi_statis:"<<
		// 	chi_statis<<endl;
		return p_value;
	}

	double chi2_test(const std::vector<double> & dis_clu)
	{
		double p_value ;
		double dis_backup;

		double sum = std::accumulate(std::begin(dis_clu), std::end(dis_clu), 0.0);  
		double mean =  sum / dis_clu.size(); //均值  
						  
		double accum  = 0.0;  
		std::for_each (std::begin(dis_clu), std::end(dis_clu), [&](const double d) {  
			accum  += (d-mean)*(d-mean);  
			});  
						  
		double chi_statis = (accum/mean); //卡方统计量
		p_value = 1-chi2_p(dis_clu.size()-1, chi_statis);	

		// cout<<"p_value:"<<p_value[0]<<"  pre_p_value:"<<p_value[1]<<"  dof:"<<(dis_clu.size()-1)<<"  chi_statis:"<<
		// 	chi_statis<<endl;
		return p_value;
	}


	void clusterize_zihao( const IntPairSet& loops, const char* filename)//, std::vector<int>& membership, int& clusterCount)
	{
		// std::map<int, std::set<int>> conflict, uncertain;
		// seperate good loops from bad loops. then calculate UCD and transform distance

		std::map<int, double>  conflictClusterDistanceMap;
		double disSTART, disEND,  pre_p_value = 0.5;
		int num_loop=0, check = 0, consCluster_serialNum = -1;
		bool findStartVertexInMatrix=false,findEndVertexInMatrix=false, ready2chiTest=false, futileBit = 0, fultileBit = 0;
		std::array<double,3> startPosition, endPosition, nearest_cluster;
		std::array<double,2> tem_dis, p_value;
		std::array<double,6> fullLoopInfo;//six elements:start ID,X,Y,end ID,X,Y
		std::pair<double, std::pair<std::pair<int, int>, int> > ele_chiInfo;
		std::set<int> splitSET, split2SET, spCONset, sp2CONset;

		std::vector<std::set<int> > conflict_cluster_set_vector;
		std::vector<int> conflict_cluster_set, uncertain_set,conflict_cluster_set_, uncertain_set_, 
			cons_cluster_number, cons_cluster_number2, cons_cluster_number_;
		std::pair<int,std::vector<int> > uncertain_conflict_ele;

		std::vector<double>  chiStatis, chiStatis2;
		std::vector<std::pair<int,std::vector<int> > >  split, split2 ;
		std::vector<std::pair<int,std::vector<double> > >  splitDis, splitDis2;

		//LC_Inf is a map, from nodes pair of loop closure to the six element array of 
		collect_vertexAndLP(filename, LC_Inf, LP_Trans_Covar_Map);

		prepare_accelerateBySynthetic(OdoInf, 
				syntheticOdo10_trans_varMatrix,
				syntheticOdo100_trans_varMatrix,
				syntheticOdo1000_trans_varMatrix,
				syntheticOdo10000_trans_varMatrix,
				syntheticOdo10_dis, 
				syntheticOdo100_dis,
				syntheticOdo1000_dis, 
				syntheticOdo10000_dis);


		prepare_accelerateBySynthetic_adjoint(OdoInf, 
			syntheticOdo10_trans_varMatrix_adjoint,
			syntheticOdo100_trans_varMatrix_adjoint,
			syntheticOdo1000_trans_varMatrix_adjoint,
			syntheticOdo10000_trans_varMatrix_adjoint,
			syntheticOdo10_dis_adjoint, syntheticOdo100_dis_adjoint,
			syntheticOdo1000_dis_adjoint, syntheticOdo10000_dis_adjoint);


		ofstream fileStreamr; 
		
		for(IntPairSet::const_iterator it = loops.begin(), lend = loops.end();it!=lend;it++)
		{
			if(set4GTL.find(*it) != set4GTL.end())
				good_loops_vector.push_back(*it);
			else
				bad_loops_vector.push_back(*it);
		}		

		// prepare_accelerateBySynthetic(OdoInf, 
		// 		syntheticOdo10_trans_varMatrix_adjoint,
		// 		syntheticOdo100_trans_varMatrix_adjoint,
		// 		syntheticOdo1000_trans_varMatrix_adjoint,
		// 		syntheticOdo10000_trans_varMatrix_adjoint,
		// 		syntheticOdo10_dis_adjoint, 
		// 		syntheticOdo100_dis_adjoint,
		// 		syntheticOdo1000_dis_adjoint, 
		// 		syntheticOdo10000_dis_adjoint);


		// fileStreamr.open("loop_serial.txt",ios::trunc);
		// for(int i = 0; i < good_loops_vector.size(); i++)
		// {
		// 	fileStreamr<<i<<" "<<good_loops_vector[i].first<<" "<<good_loops_vector[i].second<<" good"<<"\n";
		// }
		// for(int i = 0; i < bad_loops_vector.size(); i++)
		// {
		// 	fileStreamr<<i+good_loops_vector.size()<<" "<<bad_loops_vector[i].first<<" "<<bad_loops_vector[i].second<<" bad"<<"\n";
		// }
		// fileStreamr.close();

		// fileStreamr.open("index_dispersion_transfrom_distance.txt",ios::trunc);
		// fileStreamr<<"good_loops_size "<<good_loops_vector.size()<<" "<<"bad_loops_size "<<bad_loops_vector.size()<<"\n";
		// int goodLoopNum = good_loops_vector.size();

		// good_loops_vector.insert( good_loops_vector.end(), bad_loops_vector.begin(), bad_loops_vector.end() );

		// double covX, covY, ucd_two_odo, ucd_two_odo_1, ucd_two_odo_2, ucd_whole, beli = 0.95;
		// int nodeDis1, nodeDis2, nodeDis3;
		// for(int i = 0; i < good_loops_vector.size()-1; i++)
		// {
		// 	for(int j = i+1; j < good_loops_vector.size(); j++)
		// 	{
		// 		std::array<std::array<double, 5>, 5>  length;
		// 		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		// 		std::pair<int, int> loop_node_pair;
		// 		std::pair<bool, double> reV, reV_adjoint;
		// 		double  transX_residual, transY_residual, transA_residual, transformDistance, x_cov, y_cov;
		// 		bool futileBit1, futileBit2, futileBit;
		// 		Matrix3d backup;

		// 		nodeDis1 = abs(good_loops_vector[i].first - good_loops_vector[j].first);
		// 		nodeDis2 = abs(good_loops_vector[i].second - good_loops_vector[j].second);
		// 		nodeDis3 = nodeDis2 + nodeDis1;
		// 		// cout<<"dis 1 2 3"<<nodeDis1<<" "<<nodeDis2<<" "<<nodeDis3<<endl;

		// 		prepare_to_signle_loop_pair_check(good_loops_vector[i], good_loops_vector[j], FullInfo, LP_Trans_Covar_Map, 
		// 			length, futileBit1, futileBit2);
		// 		// calculate the index of dispersion for systhsized segment form two odo segements 1
		// 		ucd_two_odo   = sqrt(Cov_two_odo(0, 0)*Cov_two_odo(0, 0) + Cov_two_odo(1, 1)*Cov_two_odo(1, 1)
		// 			+ Cov_two_odo(0,1)*Cov_two_odo(0,1) + Cov_two_odo(1,0)*Cov_two_odo(1,0))/(length[0][0]+length[2][0]);

		// 		// calculate the index of dispersion for odo segement 1
		// 		x_cov = FullInfo[0].second(0,0);
		// 		y_cov = FullInfo[0].second(1,1);
		// 		ucd_two_odo_1 = sqrt(x_cov*x_cov + FullInfo[0].second(0,1)*FullInfo[0].second(0,1) 
		// 			+ FullInfo[0].second(1,0)*FullInfo[0].second(1,0) + y_cov*y_cov)/(length[0][0]);

		// 		if(length[0][0] == 0)
		// 		{
		// 			ucd_two_odo_1 = 0.01;
		// 		}

		// 		// calculate the index of dispersion for odo segement 2
		// 		x_cov = FullInfo[2].second(0,0);
		// 		y_cov = FullInfo[2].second(1,1);
		// 		ucd_two_odo_2 = sqrt(x_cov*x_cov + FullInfo[0].second(0,1)*FullInfo[0].second(0,1) 
		// 			+ FullInfo[0].second(1,0)*FullInfo[0].second(1,0) + y_cov*y_cov)/(length[2][0]);

		// 		if(length[2][0] == 0)
		// 		{
		// 			ucd_two_odo_2 = 0.01;
		// 		}

		// 		// cout<<"start loop "<<good_loops_vector[i].first<<" "<<good_loops_vector[i].second<<endl;
		// 		// cout<<"end   loop "<<good_loops_vector[j].first<<" "<<good_loops_vector[j].second<<endl;
		// 		if(std::isnan(ucd_two_odo_2))
		// 		{

		// 			cout<<"nan, the length is "<<length[2][0]<<endl;\
		// 			cin.get();
		// 		}

		// 		reV = check_single_loop_inter_varying_belief(FullInfo, covX, covY, displayCov, transX_residual, 
		// 			transY_residual, transA_residual, length[4][0], futileBit, beli);
		// 		backup = displayCov;
		// 		// calculate the index of dispersion for the whole loop
		// 		ucd_whole = index_dispersion;

		// 		prepare_to_signle_loop_pair_check_adjoint(good_loops_vector[i], good_loops_vector[j], FullInfo, LP_Trans_Covar_Map, 
		// 			length, futileBit1, futileBit2);

		// 		reV_adjoint = check_single_loop_inter_varying_belief_adjoint(FullInfo, covX, covY, displayCov, transX_residual, 
		// 			transY_residual, transA_residual, length[4][0], futileBit, beli);

		// 		// if(i < goodLoopNum and j >= goodLoopNum and (reV_adjoint.second < 8 ))  //or reV_adjoint.second <= 8
		// 		// // if(reV.second > reV_adjoint.second and (reV.second >8 and reV_adjoint.second <= 8))
		// 		// {
		// 		// 	cout<<"i: "<<i<<" j: "<<j<<endl;
		// 		// 	cout<<"original dis is smaller, "<<reV.second<<" "<<reV_adjoint.second<<endl;
		// 		// 	cout<<"cov_original"<<endl<<backup<<endl<<"adjoint: "<<endl<<displayCov<<endl;
		// 		// 	cin.get();
		// 		// }

		// 		// fileStreamr<<" "<< ucd_two_odo_1<<" "<<ucd_two_odo_2<<" "
		// 		// 	<<ucd_two_odo<<" "<<ucd_whole<<" "<<reV.second <<" "<<reV_adjoint.second<<" "<<nodeDis1<<" "<<
		// 		// 	nodeDis2<<" "<<nodeDis3<<"\n";

		// 		fileStreamr<<"i "<<i<<" j "<<j<<" "<< ucd_two_odo_1<<" "<<ucd_two_odo_2<<" "
		// 			<<ucd_two_odo<<" "<<ucd_whole<<" "<<reV.second <<" "<<reV_adjoint.second<<" "<<nodeDis1<<" "<<
		// 			nodeDis2<<" "<<nodeDis3<<"\n";
		// 	}
		// }
		// fileStreamr.close();

		// // // nodeID   trajectory_distace   tranform_distance   variance_x variance_y
		// // fileStreamr.open("indexVStrajectory.txt",ios::trunc);
		// // for(int i =0; i <= OdoInf.size(); i++)
		// // // for(int i = OdoInf.size(); i <= OdoInf.size(); i++)
		// // {
		// // 	int start, end;
		// // 	std::pair<g2o::SE2, Matrix3d> result_synthetic_odo;
		// // 	std::array<double, 5>  length;
		// // 	bool  fultileBit;
		// // 	double trans_X, trans_Y, tran_XY, varX, varY, varXY,eigen1, eigen2, eigen_XY, covxy_1, covxy_2, var_cov;
		// // 	// for(int i = 0; i*10+10 <= OdoInf.size(); i++ )
		// // 	{
		// // 		end   = i;
		// // 		synthesize_odo_edges( 0, end, OdoInf, result_synthetic_odo.first, 
		// // 			result_synthetic_odo.second, length, fultileBit);
		// // 		trans_X = result_synthetic_odo.first[0];
		// // 		trans_Y = result_synthetic_odo.first[1];
		// // 		tran_XY = sqrt(trans_X*trans_X + trans_Y*trans_Y);
		// // 		varX    = result_synthetic_odo.second(0,0);
		// // 		varY    = result_synthetic_odo.second(1,1);
		// // 		varXY   = sqrt(varX*varX + varY*varY);
		// // 		covxy_1 = result_synthetic_odo.second(0,1);
		// // 		covxy_2 = result_synthetic_odo.second(1,0);
		// // 		var_cov   = sqrt(varX*varX + varY*varY + covxy_1*covxy_1 + covxy_2*covxy_2);
		// // 		double determinent = result_synthetic_odo.second.determinant();

		// // 		EigenSolver<MatrixXd> es(result_synthetic_odo.second.block(0,0,2,2));
		// // 		// cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
		// // 		eigen1 = es.eigenvalues()[0].real();
		// // 		eigen2 = es.eigenvalues()[1].real();
		// // 		eigen_XY = sqrt(eigen1*eigen1 + eigen2*eigen2);
		// // 		// eigen_XY = 0;

		// // 		cout<<endl<<i<<" var matrix "<<endl<<result_synthetic_odo.second<<endl;
				
		// // 		if(fultileBit == 1)
		// // 		{
		// // 			futilePairVector.push_back(std::pair<int,int>(start,end));
		// // 			// cout<<"length:"<<length[0]<<" "<<length[1]<<" "<<length[2]<<" "<<length[3]<<" "<<length[4]<<" "<<length[5]<<" "<<length[6]<<" "<<endl;
		// // 			// cout<<"futile from "<<start<<" to "<<end<<endl;
		// // 		}		
		// // 		fileStreamr<<"vertexID "<<i<<" "<<length[0]<<" "<<tran_XY<<" "<<varXY<<" "<<var_cov<<" "<<eigen_XY<<" "<<determinent<<"\n";	
		// // 	}			
		// // }
		// // cout<<"please check the covariance of each node sequence"<<endl;
		// // cin.get();
		// // fileStreamr.close();

		// // nodeID   trajectory_distace   tranform_distance   variance_x variance_y
		// fileStreamr.open("indexVStrajectory_accelerated_bovisa04.txt",ios::trunc);
		// for(int i =0; i <= OdoInf.size(); i++)
		// {
		// 	int start, end;
		// 	std::pair<g2o::SE2, Matrix3d> result_synthetic_odo;
		// 	std::array<double, 5>  length;
		// 	bool  fultileBit;
		// 	double trans_X, trans_Y, tran_XY, varX, varY, varXY, covxy_1, covxy_2, var_cov;
		// 	// for(int i = 0; i*10+10 <= OdoInf.size(); i++ )
		// 	{
		// 		end   = i;
		// 		accelerated_synthetic_odo(0,  i, OdoInf, 
		// 				result_synthetic_odo,  length, fultileBit);
		// 		// cout<<"length[0]:"<<length[0]<<endl;
		// 		// cin.get();

		// 		trans_X = result_synthetic_odo.first[0];
		// 		trans_Y = result_synthetic_odo.first[1];
		// 		tran_XY = sqrt(trans_X*trans_X + trans_Y*trans_Y);
		// 		varX    = result_synthetic_odo.second(0,0);
		// 		varY    = result_synthetic_odo.second(1,1);
		// 		varXY   = sqrt(varX*varX + varY*varY);
		// 		covxy_1 = result_synthetic_odo.second(0,1);
		// 		covxy_2 = result_synthetic_odo.second(1,0);
		// 		var_cov   = sqrt(varX*varX + varY*varY + covxy_1*covxy_1 + covxy_2*covxy_2);				
		// 		double determinent = result_synthetic_odo.second.determinant();
		// 		if(fultileBit == 1)
		// 		{
		// 			futilePairVector.push_back(std::pair<int,int>(start,end));
		// 			// cout<<"length:"<<length[0]<<" "<<length[1]<<" "<<length[2]<<" "<<length[3]<<" "<<length[4]<<" "<<length[5]<<" "<<length[6]<<" "<<endl;
		// 			// cout<<"futile from "<<start<<" to "<<end<<endl;
		// 		}		
		// 		fileStreamr<<"vertexID "<<i<<" "<<length[0]<<" "<<tran_XY<<" "<<varXY<<" "<<var_cov<<" "<<determinent<<"\n";	
		// 		// cout<<"length[0]:"<<length[0]<<endl;
		// 		// cin.get();
		// 	}			
		// }
		// fileStreamr.close();	

		
		// // nodeID   trajectory_distace   tranform_distance   variance_x variance_y
		// fileStreamr.open("indexVStrajectory_adjoint.txt",ios::trunc);
		// for(int i =0; i <= OdoInf.size(); i++)
		// // for(int i = OdoInf.size(); i <= OdoInf.size(); i++)
		// {
		// 	int start, end;
		// 	std::pair<g2o::SE2, Matrix3d> result_synthetic_odo;
		// 	std::array<double, 5>  length;
		// 	bool  fultileBit;
		// 	double trans_X, trans_Y, tran_XY, varX, varY, varXY,eigen1, eigen2, eigen_XY, covxy_1, covxy_2, var_cov;
		// 	// for(int i = 0; i*10+10 <= OdoInf.size(); i++ )
		// 	{
		// 		end   = i;
		// 		synthesize_odo_edges_adjoint( 0, end, OdoInf, result_synthetic_odo.first, 
		// 			result_synthetic_odo.second, length, fultileBit);
		// 		trans_X = result_synthetic_odo.first[0];
		// 		trans_Y = result_synthetic_odo.first[1];
		// 		tran_XY = sqrt(trans_X*trans_X + trans_Y*trans_Y);
		// 		varX    = result_synthetic_odo.second(0,0);
		// 		varY    = result_synthetic_odo.second(1,1);
		// 		varXY   = sqrt(varX*varX + varY*varY);
		// 		covxy_1 = result_synthetic_odo.second(0,1);
		// 		covxy_2 = result_synthetic_odo.second(1,0);
		// 		var_cov   = sqrt(varX*varX + varY*varY + covxy_1*covxy_1 + covxy_2*covxy_2);
		// 		double determinent = result_synthetic_odo.second.determinant();

		// 		EigenSolver<MatrixXd> es(result_synthetic_odo.second.block(0,0,2,2));
		// 		// cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
		// 		eigen1 = es.eigenvalues()[0].real();
		// 		eigen2 = es.eigenvalues()[1].real();
		// 		eigen_XY = sqrt(eigen1*eigen1 + eigen2*eigen2);
		// 		// eigen_XY = 0;
		// 		cout<<endl<<i<<" var adjoint matrix "<<endl<<result_synthetic_odo.second<<endl;
				
		// 		if(fultileBit == 1)
		// 		{
		// 			futilePairVector.push_back(std::pair<int,int>(start,end));
		// 			// cout<<"length:"<<length[0]<<" "<<length[1]<<" "<<length[2]<<" "<<length[3]<<" "<<length[4]<<" "<<length[5]<<" "<<length[6]<<" "<<endl;
		// 			// cout<<"futile from "<<start<<" to "<<end<<endl;
		// 		}		
		// 		fileStreamr<<"vertexID "<<i<<" "<<length[0]<<" "<<tran_XY<<" "<<varXY<<" "<<var_cov<<" "<<eigen_XY<<" "<<determinent<<"\n";	
		// 	}			
		// }
		// cout<<"please check the covariance matrix produced by adjoint method "<<endl;
		// cin.get();
		// fileStreamr.close();



		// second_prepare_accelerated_synthetic_odo(OdoInf, syntheticOdoAccumulate_trans_varMatrix,
		// syntheticOdoAccumulate_dis);

		if(1)
		{
			double covX, covY, ucd_two_odo, ucd_two_odo_1, ucd_two_odo_2, ucd_whole, beli = 0.95;
			std::array<std::array<double, 5>, 5>  length;
			std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
			std::pair<int, int> loop_node_pair;
			std::pair<bool, double> reV, reV_adjoint;
			double  transX_residual, transY_residual, transA_residual, transformDistance, x_cov, y_cov;
			bool futileBit1, futileBit2, futileBit;
			Matrix3d backup;

			std::pair<int,int> loop1(10,975), loop2(35,667), loop3(35,667);
// loop 4 980 in cluster 0 has error 0.0658927
// loop 10 975 in cluster 0 has error 0.297352
// loop 35 667 in cluster 1 has error 0.978617


			reV =  check_single_loop_odo(loop3, LP_Trans_Covar_Map, futileBit1, beli);
			// if(reV.second < 80)
				cout<<"the self dis is "<<reV.second<<endl;	

			prepare_to_signle_loop_pair_check(loop1, loop2, FullInfo, LP_Trans_Covar_Map, 
				length, futileBit1, futileBit2);
			cin.get();
			reV = check_single_loop_inter_varying_belief(FullInfo, covX, covY, displayCov, transX_residual, 
				transY_residual, transA_residual, length[4][0], futileBit, beli);

							prepare_to_signle_loop_pair_check_adjoint(loop1, loop2, FullInfo, LP_Trans_Covar_Map, 
					length, futileBit1, futileBit2);

				reV_adjoint = check_single_loop_inter_varying_belief_adjoint(FullInfo, covX, covY, displayCov, transX_residual, 
					transY_residual, transA_residual, length[4][0], futileBit, beli);
				cout<<reV_adjoint.second<<endl;
		}
		cin.get();


		if(loops.empty())
		{
			std::cerr<<"clusterize(): "<<__LINE__<<" no loops to make clusters"<<std::endl;
			return;
		}
		_clustersFound.clear();


		// }
		for(IntPairSet::const_iterator it = loops.begin(), lend = loops.end();it!=lend;it++)
		{

			int start 	= it->first;
			int end 	= it->second;

			consCluster_serialNum = -2;
			// if(start<end)
			// {
			// 	printf("This error about start node and end node is in %s on line %d\n",  __FILE__, __LINE__);
			// 	exit(0);
			// }
			//get loop closure vextexes position
			fullLoopInfo = get_LC_Pos( start,  end);
			//print loop number and cluster id of the loop
			num_loop++;
			// cout<<" "<<endl<<"loop "<<num_loop<<" "<<start<<" "<<end<<endl;
			// cout<<"current has "<<_clustersFound.size()<<" clusters"<<endl;

			if(_clustersFound.empty())
			{
				cluster s(start,end,_clustersFound.size());
				s.positionserial.push_back(fullLoopInfo);
				
				ele_chiInfo.first = 0;
				ele_chiInfo.second.first = *it;
				ele_chiInfo.second.second = s.positionserial.size()-1;
				s.chiStatisWhenAdd2Cluster.push_back( ele_chiInfo);//
				// s.nodesInCluster.push_back(*it);

				_clustersFound.push_back(s);

						cout<<"fullLoopInfo: "<<s.positionserial[0][0]<<" "<<s.positionserial[0][1]<<" "
						<<s.positionserial[0][2]<<" "<<s.positionserial[0][3]<<
			" "<<s.positionserial[0][4]<<" "<<s.positionserial[0][5]<<endl;

				clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
				loopToClusterIDMap[*it] = _clustersFound.size()-1;

			}
			else
			{
				// if(it->first == 1441)
				// 	exit(0);
				//search for the nearest cluster to the loop
				cons_cluster_number.clear();
				chiStatis.clear();

				double strictBelief = 0.95;
				find_cons_cluster_varying_belief_all_cluster( it,  _clustersFound,  cons_cluster_number, conflict_cluster_set,
						uncertain_set, LP_Trans_Covar_Map, VertexInf, 1 , chiStatis,  futileBit, strictBelief);

				bool bad  = 0;
				for(int i = 0; i <bad_loops_vector.size(); i++)
				{
					if(bad_loops_vector[i].first == start and bad_loops_vector[i].second == end)
					{
						bad = 1;
						DLOG(INFO)<<"loop "<<start<<" "<<end<<" is bad "<<endl;
						break;
					}
				}
				if(bad == 0)
				{
					cout<<"loop "<<start<<" "<<end<<" is good ***********************8 "<<endl;
				}

				if(bad_loops_set_self_check.find(*it) != bad_loops_set_self_check.end())
				{
					loopToClusterIDMap[*it] = -3;
					continue;
				}
				// cons_cluster_number.clear();
				// std::cin.get();
				int consClusterSerial = -1, uncertainClusterSerial = -1;
				bool exitBit = 0;
				std::vector<int > signToCalWholeCluster;

				//size equal to 0, then it means find no constent cluster, so construct a new cluster
				if(cons_cluster_number.size() == 0 )
				{
					// cout<<"no consistent cluster found"<<endl;
					cluster s(start,end,_clustersFound.size());
					s.positionserial.push_back(fullLoopInfo);

					// s.chiStatisWhenAdd2Cluster.push_back( 0);
					// s.nodesInCluster.push_back(*it);
					ele_chiInfo.first = 0;
					ele_chiInfo.second.first = *it;
					ele_chiInfo.second.second = s.positionserial.size()-1;
					s.chiStatisWhenAdd2Cluster.push_back( ele_chiInfo);//

					_clustersFound.push_back(s);

					clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
					loopToClusterIDMap[*it] = _clustersFound.size()-1;
	

					consCluster_serialNum = _clustersFound.size()-1;

					//add element to the cluster pair adn distance mapp
					if(conflict_cluster_set.size() + conflictClusterDistanceMap.size() > 0)
					{
						if(conflict_cluster_set.size() != conflictClusterDistanceMap.size())
						{
							// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
							// exit(0);		
						}

						for(int addMap = 0; addMap < conflict_cluster_set.size(); addMap++)
						{
							// cout<<consCluster_serialNum<<" conflict to clsuter "<<conflict_cluster_set[addMap]<<endl;
							conflictClusterPairDistanceMap[std::pair<int,int> (conflict_cluster_set[addMap], consCluster_serialNum)] =
							 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];
							conflictClusterPairDistanceMap[std::pair<int,int> (consCluster_serialNum, conflict_cluster_set[addMap])] =
							 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];  
						}	
					}					
				}
				//find one constent cluster, add the loop to it
				else if (cons_cluster_number.size() == 1)
				{
					// cout<<"only one consistent cluster: "<<cons_cluster_number[0]<<endl;

					consCluster_serialNum = cons_cluster_number[0];
					_clustersFound[consCluster_serialNum].positionserial.push_back(fullLoopInfo);

					clusterIDtoLoopsMap[consCluster_serialNum].insert(*it);
					loopToClusterIDMap[*it] = consCluster_serialNum;
					if(consCluster_serialNum > _clustersFound.size())
					{
						cout<<"find the point that cluster ID is too big: "<<consCluster_serialNum<<endl;
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);
					}
		
					//save the chi statis when it pass the chi2 test and being add to this cluster, we save this for further check in intra check

					ele_chiInfo.first = chiStatis[0];
					ele_chiInfo.second.first = *it;
					ele_chiInfo.second.second = _clustersFound[consCluster_serialNum].positionserial.size()-1;
					_clustersFound[consCluster_serialNum].chiStatisWhenAdd2Cluster.push_back( ele_chiInfo);//
					if(_clustersFound[consCluster_serialNum].positionserial.size() == 2)
						_clustersFound[consCluster_serialNum].chiStatisWhenAdd2Cluster[0].first = chiStatis[0];

					if(chiStatis.size() != 1)
					{
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);	
					}

					//add element to the cluster pair adn distance mapp
					if(conflict_cluster_set.size() + conflictClusterDistanceMap.size() > 0)
					{
						if(conflict_cluster_set.size() != conflictClusterDistanceMap.size())
						{
							// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
							// exit(0);		
						}

						for(int addMap = 0; addMap < conflict_cluster_set.size(); addMap++)
						{
							// cout<<consCluster_serialNum<<" conflict to clsuter "<<conflict_cluster_set[addMap]<<endl;
							conflictClusterPairDistanceMap[std::pair<int,int> (conflict_cluster_set[addMap], consCluster_serialNum)] =
							 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];
							conflictClusterPairDistanceMap[std::pair<int,int> (consCluster_serialNum, conflict_cluster_set[addMap])] =
							 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];  
						}	
					}
				}
				//find more than one consistent cluster
				else 
				{
					// cout<<"find "<<cons_cluster_number.size()<<" consistent clusters"<<endl;
					
					std::sort (node_distance_cons.begin(), node_distance_cons.end(), cmp_reverse); 
					int best_cons = 0;
					if(chiStatis.size() != cons_cluster_number.size())
					{
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);	
					}

					// cout<<" consistent cluster: "<<cons_cluster_number[0]<<endl;

					for(int selectFromMultiConsCluster = 1; selectFromMultiConsCluster < chiStatis.size(); selectFromMultiConsCluster++)
					{
						// cout<<" consistent cluster: "<<cons_cluster_number[selectFromMultiConsCluster]<<endl;

						if(chiStatis[selectFromMultiConsCluster] < chiStatis[best_cons])
							best_cons = selectFromMultiConsCluster;
						if(chiStatis[selectFromMultiConsCluster] > Belief95)
						{
							cout<<"consistent cluster has a distacne as "<<chiStatis[selectFromMultiConsCluster]<<endl;
							printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
							// exit(0);	
						}
					}

					// cout<<"debug dumpted 1"<<endl;
					// find the best one from consisent candidates
					// consCluster_serialNum = cons_cluster_number[best_cons];
					consCluster_serialNum = node_distance_cons[0].first;

					clusterIDtoLoopsMap[consCluster_serialNum].insert(*it);
					loopToClusterIDMap[*it] = consCluster_serialNum;

					// cout<<"_clustersFound[consCluster_serialNum].size: "<<_clustersFound[consCluster_serialNum].positionserial.size()<<endl;
					

					_clustersFound[consCluster_serialNum].positionserial.push_back(fullLoopInfo);
					cout<<"the final cluster selected is "<<consCluster_serialNum<<endl;

					if(consCluster_serialNum > _clustersFound.size())
					{
						cout<<"find the point that cluster ID is too big: "<<consCluster_serialNum<<endl;
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);
					}

					//save the chi statis when it pass the chi2 test and being add to this cluster, we save this for further check in intra check
					ele_chiInfo.first = chiStatis[best_cons];
					ele_chiInfo.second.first = *it;
					ele_chiInfo.second.second = _clustersFound[consCluster_serialNum].positionserial.size()-1;

					_clustersFound[consCluster_serialNum].chiStatisWhenAdd2Cluster.push_back( ele_chiInfo);//
					// cout<<"debug dumpted 4"<<endl;

					if(_clustersFound[consCluster_serialNum].positionserial.size() == 2)
						_clustersFound[consCluster_serialNum].chiStatisWhenAdd2Cluster[0].first = chiStatis[best_cons];

					//put other consistent cluster into uncertain set
					for(int selectFromMultiConsCluster = 0; selectFromMultiConsCluster < chiStatis.size(); selectFromMultiConsCluster++)
					{
						if(selectFromMultiConsCluster == best_cons)
							continue;
						uncertain_set.push_back(cons_cluster_number[selectFromMultiConsCluster]);
					}

					//add element to the cluster pair adn distance mapp
					if(conflict_cluster_set.size() + conflictClusterDistanceMap.size() > 0)
					{
						if(conflict_cluster_set.size() != conflictClusterDistanceMap.size())
						{
							// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
							// exit(0);		
						}

						for(int addMap = 0; addMap < conflict_cluster_set.size(); addMap++)
						{
							// cout<<consCluster_serialNum<<" conflict to clsuter "<<conflict_cluster_set[addMap]<<endl;
							conflictClusterPairDistanceMap[std::pair<int,int> (conflict_cluster_set[addMap], consCluster_serialNum)] =
							 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];
							conflictClusterPairDistanceMap[std::pair<int,int> (consCluster_serialNum, conflict_cluster_set[addMap])] =
							 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];  
						}
					}
				}
				// update the conflict cluster pair
				updateConflictPairSet(conflict_cluster_set,  consCluster_serialNum, uncertain_set, 0);
				updateConflictCause(consCluster_serialNum, conflict_cause)	;
						
			} 
 
			// int sf = 387, ef = 515;
			// if((start == sf  and end == ef) or (end == sf  and  start == ef))
			// {
			// 		if(node_distance_cons.size() > 1)
			// 			cout<<"candidate cluster "<<node_distance_cons[0].first<<" has smallest node dis "<<node_distance_cons[0].second<<endl;
			// 		std::cout << "start == "<<start<<" end == "<<end<<", Press \'Return\' to end." << std::endl;
			// 		cout<<"consCluster_serialNum: "<<consCluster_serialNum<<endl;
			// 		std::cin.get();					
			// }

			//stop at a certain cluster id that you want to find
			// int	toFindCluster = 179;
			// if(consCluster_serialNum == toFindCluster)
			// // if((*it).first ==9016 and (*it).second ==3876)
			// {
			// 	cout<<"loop in cluster "<<toFindCluster<<" is "<<(*it).first<<" "<<(*it).second<<endl;
			// 	// cout<<"element number in cluster 10 is "<<_clustersFound[10].positionserial.size()<<
			// 	// 	" the chi is "<< ele_chiInfo.first<<endl;
			// 	// std::cin.get();
			// }

		}

		// std::array<double,6> ty={1,1,1,1,1,1};

		// std::vector<std::pair<int, int> > clusterSizeVector;	
		std::pair<int,int > pair;	
	
		// using function as comp
		
		for(int clsuterNum = _clustersFound.size()-1; clsuterNum >= 0 ; clsuterNum--)
		{
			pair.first = clsuterNum;
			pair.second = getClusterByID(clsuterNum).size();
			clusterSizeVector.push_back(pair);
		}
		std::sort (clusterSizeVector.begin(), clusterSizeVector.end(), cmp); // 12 32 45 71(26 33 53 80)
		// for(int num = 0; num < clusterSizeVector.size();num++)
		// {
		// 	cout<<"cluster "<<clusterSizeVector[num].first<<" has "<<clusterSizeVector[num].second<<" loops"<<endl;
		// }
		// cout<<" "<<endl;
		// cout<<"loops.size is "<<loops.size()<<endl;
		// cin.get();

		// std::set<std::pair<int,int > >::iterator sta = conflictPairSet.begin(), enD = conflictPairSet.end();
		// for(int num = 0; num < clusterSizeVector.size();num++)
		// {
		// 	sta = conflictPairSet.begin(), enD = conflictPairSet.end();
		// 	for(; sta != enD; sta++)
		// 	{
		// 		if((*sta).first == clusterSizeVector[num].first)
		// 			cout<<"cluster "<<(*sta).first<<" has conflict "<<(*sta).second<<endl;
		// 	}

		// }
		// cin.get();
		// save the cluster information
		std::array<double,6> ty={1,1,1,1,1,1};

		// cout<<"get the loop's id when add to cluster "<<loopToClusterIDMap[std::pair<int,int> (9016, 3876)]<<endl;
		// std::cin.get();

		// fileStream.open("clusterFile.g2o",ios::trunc);
		
		std::pair<g2o::SE2, Matrix3d> tSave;
		double cluster_id;

		for(size_t i=0; i< _clustersFound.size(); i++)
		{
			cluster_id = clusterSizeVector[i].first;
			// fileStreamr<<"clusterID "<<cluster_id<<" "<<clusterSizeVector[i].second<<" "<<i+1<<"\n";
			bool mixed = 0, initialized = 0, privious_state;
			for(std::vector<std::array<double,6>>::const_iterator itVertex = _clustersFound[cluster_id].positionserial.begin(), 
				lendVertex = _clustersFound[cluster_id].positionserial.end();itVertex!=lendVertex;itVertex++)
			{
				ty = *itVertex;

				tSave = LP_Trans_Covar_Map[std::pair<int,int> (ty[0], ty[3])];
				bool bad  = 0;

				for(int i = 0; i <bad_loops_vector.size(); i++)
				{
					if(bad_loops_vector[i].first == ty[0] and bad_loops_vector[i].second == ty[3])
					{
						bad = 1;
						break;
					}
				}
				if(initialized == 0)
				{
					privious_state = bad;
					initialized = 1;
				}
				if(initialized == 1 and privious_state!= bad)
				{
					mixed = 1;
					mixed_clusters.insert(cluster_id);
					break;
				}
			}	
		}


		fileStreamr.open("clusters.txt",ios::trunc);
	
		// // cout<<"clusters:"<<endl;
		// fileStreamr.open("cluster-all.txt",ios::trunc);

		for(size_t i=0; i< _clustersFound.size(); i++)
		{
			cluster_id = clusterSizeVector[i].first;
			if(mixed_clusters.find(cluster_id) != mixed_clusters.end())
				fileStreamr<<"clusterID "<<cluster_id<<" "<<clusterSizeVector[i].second<<" "<<i+1<<" mixed"<<"\n";
			else
				fileStreamr<<"clusterID "<<cluster_id<<" "<<clusterSizeVector[i].second<<" "<<i+1<<"\n";
			for(std::vector<std::array<double,6>>::const_iterator itVertex = _clustersFound[cluster_id].positionserial.begin(), 
				lendVertex = _clustersFound[cluster_id].positionserial.end();itVertex!=lendVertex;itVertex++)
			{
				ty = *itVertex;

				tSave = LP_Trans_Covar_Map[std::pair<int,int> (ty[0], ty[3])];

				bool bad  = 0;

				for(int i = 0; i <bad_loops_vector.size(); i++)
				{
					if(bad_loops_vector[i].first == ty[0] and bad_loops_vector[i].second == ty[3])
					{
						bad = 1;
						break;
					}
				}

				fileStreamr<<"EDGE_SE2 "<<ty[0]<<" "<<ty[3]<<" "<<tSave.first[0]<<" "<<tSave.first[1]<<" "<<tSave.first[2]
					<<" "<<1.0/tSave.second(0,0)<<" "<<0<<" "<<0<<" "<<1.0/tSave.second(1,1)<<" "<<0<<" "<<1.0/tSave.second(2,2)<<" "<<bad<<"\n";
				// fileStream<<trystdarray[0]<<"\n";
			}	
		}
		fileStreamr.close();		

		fileStreamr.open("clusters_to_use.txt",ios::trunc);
	
		// // cout<<"clusters:"<<endl;
		// fileStreamr.open("cluster-all.txt",ios::trunc);

		for(size_t i=0; i< _clustersFound.size(); i++)
		{
			cluster_id = clusterSizeVector[i].first;
			if(mixed_clusters.find(cluster_id) != mixed_clusters.end())
				fileStreamr<<"clusterID "<<cluster_id<<" "<<clusterSizeVector[i].second<<" "<<i+1<<" mixed"<<"\n";
			else
				fileStreamr<<"clusterID "<<cluster_id<<" "<<clusterSizeVector[i].second<<" "<<i+1<<"\n";
			for(std::vector<std::array<double,6>>::const_iterator itVertex = _clustersFound[cluster_id].positionserial.begin(), 
				lendVertex = _clustersFound[cluster_id].positionserial.end();itVertex!=lendVertex;itVertex++)
			{
				ty = *itVertex;

				tSave = LP_Trans_Covar_Map[std::pair<int,int> (ty[0], ty[3])];

				fileStreamr<<"EDGE_SE2 "<<ty[0]<<" "<<ty[3]<<" "<<tSave.first[0]<<" "<<tSave.first[1]<<" "<<tSave.first[2]
					<<" "<<1.0/tSave.second(0,0)<<" "<<0<<" "<<0<<" "<<1.0/tSave.second(1,1)<<" "<<0<<" "<<1.0/tSave.second(2,2)<<"\n";
				// fileStream<<trystdarray[0]<<"\n";
			}	
		}
		fileStreamr.close();	

		cout<<"finish cluster"<<endl;
		cin.get()	;
		// exit(0);
#if 0
		if(0)
		{
			std::cout<<" \% Clusters formed "<<_clustersFound.size()<<std::endl;
			std::cout<<"limits = [ "<<std::endl;
			for(size_t i=0 ; i< _clustersFound.size() ; i++)
			{
				std::cout<<i<<" -> sz "<<_clustersFound[i].size<<" :: ";
				std::cout<<" "<<_clustersFound[i].startLow<<" "<<_clustersFound[i].startHigh<<" ";
				std::cout<<" "<<_clustersFound[i].endLow<<" "<<_clustersFound[i].endHigh<<std::endl;;

			}
			std::cout<<std::endl;
			std::cout<<"]; "<<std::endl;


			std::cout<<"membership =[ ";
			for(size_t i=0; i<membership.size();i++)
			{
				std::cout<<membership[i]<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"]; "<<std::endl;
		}
#endif

	}
	

	bool ifExistInMultiCluster(std::pair<int, int> &loop)
	{
		if(LoopMultiClusterIDMap.find(loop) == LoopMultiClusterIDMap.end())
			return 0;
		else
			return 1;
	}
	bool ifExistSingleCluster(std::pair<int, int> &loop)
	{
		if(loopToClusterIDMap.find(loop) == loopToClusterIDMap.end())
			return 0;
		else
			return 1;
	}
	void updateLoopID(const std::pair<int, int> &loop, int & id)
	{
		if(LoopMultiClusterIDMap.find(loop) == LoopMultiClusterIDMap.end())//not exist in multi 
		{
			if(loopToClusterIDMap.find(loop) == loopToClusterIDMap.end()) //not exist in single
				loopToClusterIDMap[loop] = id;
			else
			{
				std::vector<int> mutiC;
				mutiC.clear();
				mutiC.push_back(loopToClusterIDMap[loop]);
				if(mutiC[0] == id)
				{
					cout<<"exist id = new id"<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
				mutiC.push_back(id);
				LoopMultiClusterIDMap[loop] = mutiC;
			}
		}
		else
		{
			// std::vector<int> temporaryMlti = LoopMultiClusterIDMap[loop];
			if(find_ele.find(LoopMultiClusterIDMap[loop], id) != -1)
			// if((temporaryMlti).find(id) != temporaryMlti.end())
			{
				cout<<"exist id = new id"<<endl;
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			LoopMultiClusterIDMap[loop].push_back(id);
		}
	}
	
	int collect_vertexAndLP(const char* filename, std::map<std::pair<int, int>, std::array<double,9> > &LC_Inf,
		std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	& LP_Trans_Covar_Map)
	{

		ifstream fileStream;  

	    string tmp,temp1;  
	    std::array<double,4> verT={0, 0, 0, 0};
	    std::array<double,11> odoedge_element;
	    std::array<double,11> savemid;
	    std::array<double,9> lcedge_element;
	    std::pair<int, int>  lc_vertex_pair;
	    std::pair<g2o::SE2, Matrix3d> ele_lp_trans_covar;
	    // char* seg;
	    int count = 0,dddlndgsfdgj=0;// 行数计数器  
	    int position = 0, position2 = 0;  
	    double nul;

	    // fileStream.open("B25b_0.500.g2o",ios::in);//ios::in 表示以只读的方式读取文件
	    fileStream.open(filename,ios::in);  
	    if(fileStream.fail())//文件打开失败:返回0  
	    { 
	    	cout<<"open file failed"<<endl; 
	    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
	    	exit(0);
	        return 0;  
	    }  
	    else//文件存在  
	    {  
		    while(getline(fileStream,tmp,'\n'))//读取一行  
		    {  
		    	position = tmp.find("VERTEX_SE2");
		    	position2 = tmp.find("EDGE_SE2");
		    	// cout<<tmp<<endl; 	
		    	// exit(0);
		
		    	istringstream stringin(tmp);
		    	if (tmp.size() > 0 )
		    	{
					if(position!=string::npos)  
			    	{			    		
			    		for (dddlndgsfdgj=0;dddlndgsfdgj<5;dddlndgsfdgj++)
			    		{
			    			switch(dddlndgsfdgj)
			    			{
		                    	case 0:
			    					stringin >> temp1;
			    					// cout<<temp1<<endl;
			    					break;
			    				case 1:
			    				case 2:
			    				case 3:
			    				case 4:
			    					stringin >> verT[dddlndgsfdgj-1];
			    					break;
			    				default:
			    					break;
			    			}
			    		}
			    	
			    		VertexInf.push_back(verT);
			    		// cout<<verT[0]<<" "<<verT[1]<<" "<<verT[2]<<" "<<verT[3]<<endl;
			    		// cout<<"vertex size: " <<VertexInf.size()<<endl;
			    		// cout<<VertexInf[0][0]<<" "<<VertexInf[0][1]<<" "<<VertexInf[0][2]<<" "<<VertexInf[0][3]<<endl;
			    		// exit(0);
			    	}
			    	else if(position2 !=string::npos)
			    	{
			    		bool odobit = 0;

						for (dddlndgsfdgj = 0; dddlndgsfdgj < 12; dddlndgsfdgj++)
			    		{	
			    			switch(dddlndgsfdgj)
			    			{
		                    	case 0:
			    					stringin >> temp1;
			    					// cout<<temp1<<endl;
			    					break;
			    				case 1://start vertex numver
			    					stringin >> lc_vertex_pair.first;
			    					break;
			    				case 2://end vertex number
			    					stringin >> lc_vertex_pair.second;
			    					if(std::abs(lc_vertex_pair.first - lc_vertex_pair.second) == 1)
			    						odobit = 1;
			    					break;
			    				case 3://x of the edge
			    					if(odobit)//if its a odometry edge put in vector[8]
			    					{
			    						odoedge_element[0] = lc_vertex_pair.first;
			    						odoedge_element[1] = lc_vertex_pair.second;
			    						stringin >> odoedge_element[2];
			    					}
			    					else
			    						stringin >> lcedge_element[0];
			    					break;
			    				case 4: //y of the edge
			    				case 5: //andlge of the adge
			    				case 6: //information matrix (0,0)	
			    				case 7: //information matrix (0,1)	
			    				case 8: //information matrix (0,2)	
			    				case 10://information matrix (1,1)	
			    				case 9: //information matrix (1,2)
			    				case 11://information matrix (2,2)
			    					if(odobit)//if its a odometry edge put in vector[8]
			    						stringin >> odoedge_element[dddlndgsfdgj-1];
			    					else
			    						stringin >> lcedge_element[dddlndgsfdgj-3];
			    					break;	
			    				default:
			    				{
			    					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
									exit(0);
			    					break;
			    				}
			    			}
			    		}
			    		
			    		if(odobit)
			    			OdoInf.push_back(odoedge_element);
			    		else
			    		{
			    			count++; 
			    			LC_Inf[lc_vertex_pair] = lcedge_element;

							Matrix3d m1 = Matrix3d::Identity();
							g2o::Vector3 mid_vector3; 
							g2o::SE2 edge1;

							m1(0,0) = lcedge_element[3]; 
							m1(0,1) = lcedge_element[4]; 
							m1(1,0) = lcedge_element[4]; 

							m1(0,2) = lcedge_element[5];
							m1(2,0) = lcedge_element[5];

							m1(1,1 )= lcedge_element[6]; 

							m1(1,2) = lcedge_element[7]; 
							m1(2,1) = lcedge_element[7]; 

							m1(2,2) = lcedge_element[8];

							Matrix3d mx;
							mx = m1.inverse();

							mid_vector3[0] = lcedge_element[0];
							mid_vector3[1] = lcedge_element[1];
							mid_vector3[2] = lcedge_element[2];
							edge1.fromVector(mid_vector3);

							ele_lp_trans_covar.first  = edge1;
							ele_lp_trans_covar.second = mx;

							LP_Trans_Covar_Map[lc_vertex_pair] = ele_lp_trans_covar;

			    			if(lc_vertex_pair.first <lc_vertex_pair.second){
			    				// printf("loop closure from node id < to node id");
			    				// printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
								// exit(0);
			    			}
			    		}

			    	}
		    	} 
		    } 

		    // check if odo and vertex vector are sequential
		    bool esx = 0;
		    for(int i_ = 0; i_ < VertexInf.size()-1; i_++)
		    {
		    	if(VertexInf[i_][0] != i_)
		    	{
		    		esx = 1;
		    		cout<<"vertex info not match "<<"i = "<<i_<<", but vertex info[0] is "<<VertexInf[i_][0]<<endl;
		    	}
		    	if(OdoInf[i_][0] != OdoInf[i_][1]-1 or OdoInf[i_][0] != i_)
		    	{
		    		esx = 1;
		    		cout<<"odo info[0] [3] i are "<<OdoInf[i_][0]<<" "<<OdoInf[i_][1]<<" "<<OdoInf[i_][0]<<" "<<i_<<endl;
		    	}
		    }
		    if(esx == 1)
		    	exit(0);
		    if(OdoInf.size() != VertexInf.size()-1)
		    {
		    	cout<<"odo size: "<<OdoInf.size()<<" vertex size: "<<VertexInf.size()<<endl;
		    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }
		    for(int i=0;i<OdoInf.size();i++)
		    {
		    	if(std::abs(OdoInf[i][0] - OdoInf[i][1]) != 1 )
		    	{
		    		cout<<"OdoInf[i][0]: "<<OdoInf[i][0]<<" OdoInf[i][1]: "<<OdoInf[i][1]<<endl;
		    		cout<<"i: "<<i<<endl;
			    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
		    	}
		    	else if(OdoInf[i][0] != i)
		    	{
			    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);	
		    	}
		    	// cout<<OdoInf[i][0]<<" "<<OdoInf[i][1]<<" "<<OdoInf[i][2]<<" "<<OdoInf[i][3]<<" "<<OdoInf[i][4]<<" "
		    	// 	<<OdoInf[i][5]<<" "<<OdoInf[i][6]<<" "<<OdoInf[i][7]<<endl;
		    }

		    printf("loop count:%d \n",count);
		    printf("sum:%d \n",int(VertexInf.size()+LC_Inf.size()+OdoInf.size()));//+LC_Inf.size() +OdoInf.size()
		    cout<<"VertexInf.size():"<<VertexInf.size()<<endl;
		    fileStream.close();  

		    cout<<filename<<endl;
		    string calLineOfTLC, TLCFileName;  
		    calLineOfTLC = filename;
		    cout<<calLineOfTLC<<endl;

		    // int ret3 = calLineOfTLC.find("_alpha");
		    // cout<<calLineOfTLC<<endl;
		    // cout<<calLineOfTLC.find("_alpha")<<endl;

	    	TLCFileName = calLineOfTLC;//delete the last four chars and then add the part of GT LC.txt
	    	int pos = calLineOfTLC.find(".g2o");
	    	TLCFileName.replace(pos,4,"_GT_LC.txt");
				fileStream.open(TLCFileName,ios::in);  
		    if(fileStream.fail())//文件打开失败:返回0  
		    {  
		    	cout<<"GT Loop txt file exists but can not be open !"<<endl;
		    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
		    	// exit(0);
		    }  
		    else
		    {
		   //  	NumOfRealLoops;//number of the true loop closures , namely the line numbers in the txt file
		   //  	int pos = calLineOfTLC.find(".g2o");
		   //  	if(pos == -1)
		   //  	{
		   //  		cout<<"can not find \' .g2o\' in input file name ,so cant't calculate the number of true loops"<<endl;
			  //   	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					// exit(0);
		   //  	}
		   //  	TLCFileName = calLineOfTLC;//delete the last four chars and then add the part of GT LC.txt
		   //  	TLCFileName.replace(pos,4,"_GT_LC.txt");
 				
 				// fileStream.open(TLCFileName,ios::in);  
			  //   if(fileStream.fail())//文件打开失败:返回0  
			  //   {  
			  //   	cout<<"GT Loop txt file exists but can not be open !"<<endl;
			  //   	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			  //   	exit(0);
			  //   }  
			  //   else//文件存在  
			    {  
			    	NumOfRealLoops = 0;
			    	std::pair<int, int> ele_GTL;
				    while(getline(fileStream,tmp,'\n'))//读取一行  
				    { 
				    	istringstream stringin(tmp);
				    	if (tmp.size() > 0 )
				    	{
					    		for (int i=0; i < 2;i++)
					    		{
					    			// cout<<"tmp["<<i<<"] = "<<tmp[i]<<endl;
					    			switch(i)
					    			{
				                    	case 0:
					    					stringin >> ele_GTL.first;
					    					// cout<<temp1<<endl;
					    					break;
					    				case 1:
					    					stringin >> ele_GTL.second;
					    					break;
					    			}
					    		}
					    }
					    // exit(0);
				    	set4GTL.insert(ele_GTL);
						NumOfRealLoops = NumOfRealLoops+1;
						// cout<<ele_GTL.first<<" "<<ele_GTL.second<<" "<< NumOfRealLoops<<" "<< set4GTL.size()<<endl;
				    } 		    	
				}
				fileStream.close(); 
			} 
			if(NumOfRealLoops != set4GTL.size())
			{
				if(NumOfRealLoops == -1)
					cout<<"no ground truth LC file available"<<endl;
				else
				{
			    	cout<<"the number of true loops does not consistent!"<<endl;
			    	cout<<"NumOfRealLoops: "<< NumOfRealLoops<<" set4GTL.size(): "<<set4GTL.size() <<endl;
			    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			    	exit(0);	
				}
			}
		    // calLineOfTLC
		    // cout<<ret3<<endl; 
		    cout<<"final name: "<<TLCFileName<<" has "<< NumOfRealLoops<<" true loops."<<endl;
		    // exit(0);
		    return count ;  
		} 
	}

	void prepare_accelerateBySynthetic(std::vector<std::array<double,11>> & OdoInf, 
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo10_trans_varMatrix,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo100_trans_varMatrix,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo1000_trans_varMatrix,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo10000_trans_varMatrix,
		std::vector<double > & syntheticOdo10_dis, std::vector<double > & syntheticOdo100_dis,
		std::vector<double > & syntheticOdo1000_dis, std::vector<double > & syntheticOdo10000_dis)
	{
		int start, end;
		std::pair<g2o::SE2, Matrix3d> result_synthetic_odo;
		std::array<double, 5>  length;
		bool  fultileBit;
		for(int i = 0; i*10+10 <= OdoInf.size(); i++ )
		{
			start = i*10;
			end   = i*10+10;
			synthesize_odo_edges( start, end, OdoInf, result_synthetic_odo.first, 
				result_synthetic_odo.second, length, fultileBit);
			syntheticOdo10_trans_varMatrix.push_back(result_synthetic_odo);
			syntheticOdo10_dis.push_back(length[0]);

			if(fultileBit == 1)
			{
				futilePairVector.push_back(std::pair<int,int>(start,end));
				// cout<<"length:"<<length[0]<<" "<<length[1]<<" "<<length[2]<<" "<<length[3]<<" "<<length[4]<<" "<<length[5]<<" "<<length[6]<<" "<<endl;
				// cout<<"futile from "<<start<<" to "<<end<<endl;
			}			
		}
		for(int i = 0; i*100+100 <= OdoInf.size(); i++ )
		{
			start = i*100;
			end   = i*100+100;
			synthesize_odo_edges( start, end, OdoInf, result_synthetic_odo.first, 
				result_synthetic_odo.second, length, fultileBit);
			syntheticOdo100_trans_varMatrix.push_back(result_synthetic_odo);
			syntheticOdo100_dis.push_back(length[0]);

			if(fultileBit == 1)
			{
				futilePairVector.push_back(std::pair<int,int>(start,end));
				// cout<<"length:"<<length[0]<<" "<<length[1]<<" "<<length[2]<<" "<<length[3]<<" "<<length[4]<<" "<<length[5]<<" "<<length[6]<<" "<<endl;
				cout<<"futile from "<<start<<" to "<<end<<endl;
			}			
		}
		for(int i = 0; i*1000+1000 <= OdoInf.size(); i++ )
		{
			start = i*1000;
			end   = i*1000+1000;
			synthesize_odo_edges( start, end, OdoInf, result_synthetic_odo.first,
				result_synthetic_odo.second, length, fultileBit);
			syntheticOdo1000_trans_varMatrix.push_back(result_synthetic_odo);
			syntheticOdo1000_dis.push_back(length[0]);

			if(fultileBit == 1)
			{
				futilePairVector.push_back(std::pair<int,int>(start,end));
				// cout<<"length:"<<length[0]<<" "<<length[1]<<" "<<length[2]<<" "<<length[3]<<" "<<length[4]<<" "<<length[5]<<" "<<length[6]<<" "<<endl;
				cout<<"futile from "<<start<<" to "<<end<<endl;
			}			
		}
		for(int i = 0; i*10000+10000 <= OdoInf.size(); i++ )
		{
			start = i*10000;
			end   = i*10000+10000;
			synthesize_odo_edges( start, end, OdoInf, result_synthetic_odo.first,
				result_synthetic_odo.second, length, fultileBit);
			syntheticOdo10000_trans_varMatrix.push_back(result_synthetic_odo);
			syntheticOdo10000_dis.push_back(length[0]);

			if(fultileBit == 1)
			{
				futilePairVector.push_back(std::pair<int,int>(start,end));
				// cout<<"length:"<<length[0]<<" "<<length[1]<<" "<<length[2]<<" "<<length[3]<<" "<<length[4]<<" "<<length[5]<<" "<<length[6]<<" "<<endl;
				cout<<"futile from "<<start<<" to "<<end<<endl;
			}			
		}
	}

	void prepare_accelerateBySynthetic_adjoint(std::vector<std::array<double,11>> & OdoInf, 
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo10_trans_varMatrix_adjoint,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo100_trans_varMatrix_adjoint,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo1000_trans_varMatrix_adjoint,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo10000_trans_varMatrix_adjoint,
		std::vector<double > & syntheticOdo10_dis_adjoint, std::vector<double > & syntheticOdo100_dis_adjoint,
		std::vector<double > & syntheticOdo1000_dis_adjoint, std::vector<double > & syntheticOdo10000_dis_adjoint)
	{
		int start, end;
		std::pair<g2o::SE2, Matrix3d> result_synthetic_odo;
		std::array<double, 5>  length;
		bool  fultileBit;
		for(int i = 0; i*10+10 <= OdoInf.size(); i++ )
		{
			start = i*10;
			end   = i*10+10;
			synthesize_odo_edges_adjoint( start, end, OdoInf, result_synthetic_odo.first, 
				result_synthetic_odo.second, length, fultileBit);
			syntheticOdo10_trans_varMatrix_adjoint.push_back(result_synthetic_odo);
			syntheticOdo10_dis_adjoint.push_back(length[0]);

			if(fultileBit == 1)
			{
				futilePairVector.push_back(std::pair<int,int>(start,end));
				// cout<<"length:"<<length[0]<<" "<<length[1]<<" "<<length[2]<<" "<<length[3]<<" "<<length[4]<<" "<<length[5]<<" "<<length[6]<<" "<<endl;
				// cout<<"futile from "<<start<<" to "<<end<<endl;
			}			
		}
		for(int i = 0; i*100+100 <= OdoInf.size(); i++ )
		{
			start = i*100;
			end   = i*100+100;
			synthesize_odo_edges_adjoint( start, end, OdoInf, result_synthetic_odo.first, 
				result_synthetic_odo.second, length, fultileBit);
			syntheticOdo100_trans_varMatrix_adjoint.push_back(result_synthetic_odo);
			syntheticOdo100_dis_adjoint.push_back(length[0]);

			if(fultileBit == 1)
			{
				futilePairVector.push_back(std::pair<int,int>(start,end));
				// cout<<"length:"<<length[0]<<" "<<length[1]<<" "<<length[2]<<" "<<length[3]<<" "<<length[4]<<" "<<length[5]<<" "<<length[6]<<" "<<endl;
				cout<<"futile from "<<start<<" to "<<end<<endl;
			}			
		}
		for(int i = 0; i*1000+1000 <= OdoInf.size(); i++ )
		{
			start = i*1000;
			end   = i*1000+1000;
			synthesize_odo_edges_adjoint( start, end, OdoInf, result_synthetic_odo.first,
				result_synthetic_odo.second, length, fultileBit);
			syntheticOdo1000_trans_varMatrix_adjoint.push_back(result_synthetic_odo);
			syntheticOdo1000_dis_adjoint.push_back(length[0]);

			if(fultileBit == 1)
			{
				futilePairVector.push_back(std::pair<int,int>(start,end));
				// cout<<"length:"<<length[0]<<" "<<length[1]<<" "<<length[2]<<" "<<length[3]<<" "<<length[4]<<" "<<length[5]<<" "<<length[6]<<" "<<endl;
				cout<<"futile from "<<start<<" to "<<end<<endl;
			}			
		}
		for(int i = 0; i*10000+10000 <= OdoInf.size(); i++ )
		{
			start = i*10000;
			end   = i*10000+10000;
			synthesize_odo_edges_adjoint( start, end, OdoInf, result_synthetic_odo.first,
				result_synthetic_odo.second, length, fultileBit);
			syntheticOdo10000_trans_varMatrix_adjoint.push_back(result_synthetic_odo);
			syntheticOdo10000_dis_adjoint.push_back(length[0]);

			if(fultileBit == 1)
			{
				futilePairVector.push_back(std::pair<int,int>(start,end));
				// cout<<"length:"<<length[0]<<" "<<length[1]<<" "<<length[2]<<" "<<length[3]<<" "<<length[4]<<" "<<length[5]<<" "<<length[6]<<" "<<endl;
				cout<<"futile from "<<start<<" to "<<end<<endl;
			}			
		}
	}

	void second_prepare_accelerated_synthetic_odo(std::vector<std::array<double,11>> & OdoInf, 
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdoAccumulate_trans_varMatrix,
		std::vector<double> & syntheticOdoAccumulate_dis)
	{
		int i = 0;
		std::array<double, 4> length_ele;
		std::pair<g2o::SE2, Matrix3d>  result_synthetic_odo, result_synthetic_odo2;
		std::array<double, 5>  length;
		bool fultileBit;
		int odoSize = OdoInf.size();

		while(i <= OdoInf.size())
		{
			
			accelerated_synthetic_odo(0,  i, OdoInf, 
			result_synthetic_odo,  length, fultileBit);

			syntheticOdoAccumulate_dis.push_back(length[0]);
			syntheticOdoAccumulate_trans_varMatrix.push_back(result_synthetic_odo);

			synthesize_odo_edges( 0, i, OdoInf, result_synthetic_odo2.first, 
				result_synthetic_odo2.second, length, fultileBit);

			judgeTwoSE2Equal(result_synthetic_odo.first, result_synthetic_odo2.first);

			if(judgeTwoMatrix3dEqual(result_synthetic_odo.second,result_synthetic_odo2.second))
			{
				cout<<"start: 0 "<<" end: "<<i<<endl;
				cout<<"cov2:"<<endl;
				cout<<result_synthetic_odo.second<<endl;
				cout<<"cov1_:"<<endl;
				cout<<result_synthetic_odo2.second<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);			
			}
			i++;
		}
		// ;
		// result_synhetic_odo2.second;

	}



	// void second_accelerated_synthetic_odo(int start, int end,
	// 	std::pair<g2o::SE2, Matrix3d> & result_synthetic_odo, std::array<double, 5>  & length,
	// 	std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdoAccumulate_trans_varMatrix,
	// 	std::vector<double > & syntheticOdoAccumulate_dis,
	// 	bool & fultileBit)
	// {	
	// 	g2o::SE2 midTrans, nextTrans, SE; 
	// 	Matrix3d nextVarMatrix, midVarMatrix, j1, j2, cov;		

	// 	g2o::SE2 toFindOutWhyAcceleratedNotEqualSlow_SE2;
	// 	Matrix3d toFindOutWhyAcceleratedNotEqualSlow_VAR;
	// 	bool toFindOutWhyAcceleratedNotEqualSlow_fultileBit;
	// 	std::array<double,5> toFindOutWhyAcceleratedNotEqualSlow_length;

	// 	if(start == 0)
	// 	{
	// 		result_synthetic_odo = syntheticOdoAccumulate_trans_varMatrix[end];
	// 		SE = syntheticOdoAccumulate_trans_varMatrix[end].first;
	// 	}
	// 	else
	// 	{
	// 		midTrans     = syntheticOdoAccumulate_trans_varMatrix[start].first.inverse();
	// 		SE = midTrans;
	// 		midVarMatrix = syntheticOdoAccumulate_trans_varMatrix[start].second;
		    
	// 	    // cerr<<"vitualStart: "<<vitualStart <<endl;
	// 		Jacobian_4_edge_propagate(midTrans, syntheticOdoAccumulate_trans_varMatrix[end].first, j1, j2);
	// 		covariance_propagate(midVarMatrix, syntheticOdoAccumulate_trans_varMatrix[end].second, j1, j2, nextVarMatrix);
	// 		result_synthetic_odo.first = midTrans*(syntheticOdoAccumulate_trans_varMatrix[end].first);
	// 		result_synthetic_odo.second = nextVarMatrix;
	// 		length[0] = syntheticOdoAccumulate_dis[end] - syntheticOdoAccumulate_dis[start];
	// 		length[2] = syntheticOdoAccumulate_dis[end] - syntheticOdoAccumulate_dis[start];
	// 		length[3] = syntheticOdoAccumulate_dis[end] - syntheticOdoAccumulate_dis[start];
	// 		length[4] = syntheticOdoAccumulate_dis[end] - syntheticOdoAccumulate_dis[start];			
	// 	}
	// 	// fultileBit = 0;
	// 	fultileBit = judgeVarianceAndDistance_Small(length, midVarMatrix);

	// 	cout<<"midVarMatrix: "<<endl;
	// 	cout<<midVarMatrix;
	// 	cout<<"nextVarMatrix: "<<endl;	
	// 	cout<<nextVarMatrix<<endl;
		
	// }

	void accelerated_synthetic_odo_adjoint(int start, int end, std::vector<std::array<double,11>> & OdoInf, 
		std::pair<g2o::SE2, Matrix3d> & result_synthetic_odo, std::array<double, 5>  & length,
		bool & fultileBit)
	{	
		g2o::SE2 midTrans; 
		Matrix3d j1, j2,j3, cov;		
		std::pair<g2o::SE2, Matrix3d> next_pair;
		double distance_segment; 
		int a ,b ,c;

		if (start > 99999 or end > 99999)
		{
			cout<<"start "<<start<<" end "<<end<<endl;
			cout<<"the number of vectex is too big"<<endl;
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		if(start > OdoInf.size() or end > OdoInf.size())
		{
			cout<<"start: "<<start<<" end:"<<end <<" OdoInf.size:"<<OdoInf.size() <<endl;
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}

		int startWan = 0, startQian= 0, startBai = 0, startShi = 0, startGe  = 0, 
			endWan   = 0,  endQian = 0,  endBai  = 0,  endShi  = 0,  endGe   = 0;

		get_digits(start, end, startWan, startQian, startBai, startShi, startGe, endWan, endQian, endBai, endShi, endGe);		

		int vitualStart = start , virtualEnd =  end;	
		g2o::SE2 toFindOutWhyAcceleratedNotEqualSlow_SE2;
		Matrix3d toFindOutWhyAcceleratedNotEqualSlow_VAR;
		bool toFindOutWhyAcceleratedNotEqualSlow_fultileBit;
		std::array<double,5> toFindOutWhyAcceleratedNotEqualSlow_length;

		// cout<<"start: "<<start<<endl;
		//ge
		if ((start - startGe +20) > end) 
		    synthesize_odo_edges_adjoint( start, end, OdoInf, result_synthetic_odo.first, result_synthetic_odo.second, length, fultileBit);
		//shi
		else if ((start - startGe- startShi*10 +200) > end) 
		{
			int increaseToTen ;
			if(startGe != 0)
				increaseToTen = start - startGe +10;
			else
				increaseToTen = start;
			vitualStart = increaseToTen;
			// cerr<<"vitualStart: "<<vitualStart <<endl;
			synthesize_odo_edges_adjoint( start, increaseToTen, OdoInf, result_synthetic_odo.first, 
				result_synthetic_odo.second, length, fultileBit);	
			int decreaseToTen;
			decreaseToTen = (end - endGe)/10 -1;
			//should be an integer multiple of ten
			if((end - endGe -(start - startGe +10)) % 10 != 0)
			{
				cout<<"the number of vectex is too big"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			//from start to the nearest  integel multiple of ten
		    for(int tenSerial = (increaseToTen) / 10; tenSerial <= decreaseToTen; tenSerial++)
		    { 
		    	vitualStart = vitualStart + 10;  	
		    	// cerr<<"vitualStart: "<<vitualStart <<endl;
				adjoint_update(result_synthetic_odo, syntheticOdo10_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo10_dis_adjoint[tenSerial]);
		    }
		    //from the last  integel multiple of ten to end
		    g2o::Vector3 mid_vector3;
		    std::array<double,11> startNodeInfo;
		    for(int tenSerial = end - endGe; tenSerial < end; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1;   
		    	// cout<<"vitualStart: "<<vitualStart <<endl;
		    	startNodeInfo = OdoInf[tenSerial];

				cov(0,0) = startNodeInfo[5]; 
				cov(0,1) = startNodeInfo[6];
				cov(1,0) = startNodeInfo[6];  
				cov(0,2) = startNodeInfo[7];
				cov(2,0) = startNodeInfo[7];
				cov(1,1) = startNodeInfo[8]; 
				cov(1,2) = startNodeInfo[9]; 
				cov(2,1) = startNodeInfo[9]; 
				cov(2,2) = startNodeInfo[10];

				mid_vector3[0] = startNodeInfo[2];
				mid_vector3[1] = startNodeInfo[3];
				mid_vector3[2] = startNodeInfo[4];

				next_pair.first.fromVector(mid_vector3);	
				next_pair.second = cov.inverse();

				distance_segment = sqrt(startNodeInfo[2]*startNodeInfo[2] + startNodeInfo[3]*startNodeInfo[3] );
				adjoint_update(result_synthetic_odo, next_pair, j1, j2, j3, midTrans, length, distance_segment);
		    }
		    if(vitualStart != end)
		    {
				cout<<"the final node does not equal to end vertex"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }
		    fultileBit = judgeVarianceAndDistance_Small(length, result_synthetic_odo.second);
		    // fultileBit = judgeVarianceAndDistance_wholeMatrix(length, cov);
		}
		//bai wei 
		else if ((start - startGe- startShi*10 - startBai*100+2000) > end)
		{
			int increaseToTen ;
			if(startGe != 0)
				increaseToTen = start - startGe +10;
			else
				increaseToTen = start;
			//bai ge
			synthesize_odo_edges_adjoint( start, increaseToTen, OdoInf, result_synthetic_odo.first, 
				result_synthetic_odo.second, length, fultileBit);

			vitualStart = increaseToTen; 
			// cout<<"vitualStart: "<<vitualStart <<endl;
			// cout<<"cov:"<<endl<<result_synthetic_odo.second<<endl;

			//should be an integer multiple of ten
			if((end - endGe -(start - startGe +10)) % 10 != 0)
			{
				cout<<"the number of vectex is too big"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			//bai shi from start to the nearest  integel multiple of ten
			a = increaseToTen / 10;
			b = (start -startGe - startShi*10 +100) /10 -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10;  	
		
				adjoint_update(result_synthetic_odo, syntheticOdo10_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo10_dis_adjoint[tenSerial]);

				// cout<<"vitualStart: "<<vitualStart <<endl;
				// cout<<"cov:"<<endl<<result_synthetic_odo.second<<endl;
				// cout<<"cov_incre_sements:"<<endl<<syntheticOdo10_trans_varMatrix_adjoint[tenSerial].second<<endl;
		    }
		    //bai to bai
		    a = (start -startGe - startShi*10 +100) /100;
		    b = ((end -endGe - endShi*10) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 100;  	
				adjoint_update(result_synthetic_odo, syntheticOdo100_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo100_dis_adjoint[tenSerial]);
				// cout<<"vitualStart: "<<vitualStart <<endl;
				// cout<<"cov:"<<endl<<result_synthetic_odo.second<<endl;
		    }
		    //bai to shi
		    a = (end - endGe - endShi*10) /10;
		    b = ((end - endGe) /10) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10;  	
				adjoint_update(result_synthetic_odo, syntheticOdo10_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo10_dis_adjoint[tenSerial]);
				// cout<<"vitualStart: "<<vitualStart <<endl;
				// cout<<"cov:"<<endl<<result_synthetic_odo.second<<endl;				
		    }
		    //from the last  integel multiple of ten to end
		    g2o::Vector3 mid_vector3;
		    std::array<double,11> startNodeInfo;
		    for(int tenSerial = end - endGe; tenSerial < end; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1; 
		    	startNodeInfo = OdoInf[tenSerial];

				cov(0,0) = startNodeInfo[5]; 
				cov(0,1) = startNodeInfo[6];
				cov(1,0) = startNodeInfo[6];  
				cov(0,2) = startNodeInfo[7];
				cov(2,0) = startNodeInfo[7];
				cov(1,1) = startNodeInfo[8]; 
				cov(1,2) = startNodeInfo[9]; 
				cov(2,1) = startNodeInfo[9]; 
				cov(2,2) = startNodeInfo[10];

				mid_vector3[0] = startNodeInfo[2];
				mid_vector3[1] = startNodeInfo[3];
				mid_vector3[2] = startNodeInfo[4];

				next_pair.first.fromVector(mid_vector3);	
				next_pair.second = cov.inverse();

				distance_segment = sqrt(startNodeInfo[2]*startNodeInfo[2] + startNodeInfo[3]*startNodeInfo[3] );
				adjoint_update(result_synthetic_odo, next_pair, j1, j2, j3, midTrans, length, distance_segment);
				// cout<<"vitualStart: "<<vitualStart <<endl;
				// cout<<"cov:"<<endl<<result_synthetic_odo.second<<endl;					
		    }
		    if(vitualStart != end)
		    {
				cout<<"the final node does not equal to end vertex"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }
		    fultileBit = judgeVarianceAndDistance_Small(length, result_synthetic_odo.second);
		}
		//thousand 
		else if ((start - startGe- startShi*10 - startBai*100 - startQian*1000+20000) > end)
		{
			int increaseToTen ;
			if(startGe != 0)
				increaseToTen = start - startGe +10;
			else
				increaseToTen = start;
			synthesize_odo_edges_adjoint( start, increaseToTen, OdoInf, result_synthetic_odo.first, result_synthetic_odo.second, length, fultileBit);
			vitualStart = increaseToTen;
			// cout<<"start: "<<start<<endl;
			// cout<<"ge shi bai qian: "<<startGe<<" "<<startShi<<" "<<startBai<<" "<<startQian<<endl;
			// cout<<"vitualStart: "<<vitualStart <<endl;

			//should be an integer multiple of ten
			if((end - endGe -(start - startGe +10)) % 10 != 0)
			{
				cout<<"the number of vectex is too big"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			// cout<<"result_synthetic_odo.second: "<< result_synthetic_odo.second<<endl;
			// synthesize_odo_edges_adjoint( start, vitualStart, OdoInf, toFindOutWhyAcceleratedNotEqualSlow_SE2, 
			// 	toFindOutWhyAcceleratedNotEqualSlow_VAR, toFindOutWhyAcceleratedNotEqualSlow_length, 
			// 	toFindOutWhyAcceleratedNotEqualSlow_fultileBit);

			//shi to bai
			a = increaseToTen / 10;
			b = (start -startGe - startShi*10 +100) /10 -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {   
		    	vitualStart = vitualStart + 10;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;							
				adjoint_update(result_synthetic_odo, syntheticOdo10_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo10_dis_adjoint[tenSerial]);
		    }
		    //bai to qian
		    a = (start -startGe - startShi*10 +100) /100;
		    b = ((start -startGe - startShi*10 - startBai*100 +1000 ) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {   
		    	vitualStart = vitualStart + 100;	
				adjoint_update(result_synthetic_odo, syntheticOdo100_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo100_dis_adjoint[tenSerial]);
		    }
		    //qian to qian
		    a = (start -startGe - startShi*10 - startBai*100 +1000) /1000;
		    b = ((end -endGe - endShi*10 - endBai*100) /1000) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1000;	
		    	// cout<<"a: "<<a<<endl; 
				adjoint_update(result_synthetic_odo, syntheticOdo1000_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo1000_dis_adjoint[tenSerial]);
		    }
		    //qian to bai
		    a = ((end -endGe - endShi*10 - endBai*100) /100) ;
		    b = ((end -endGe - endShi*10 ) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 100;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl; 								
				adjoint_update(result_synthetic_odo, syntheticOdo100_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo100_dis_adjoint[tenSerial]);
		    }
		    //bai to shi
		    a = ((end -endGe - endShi*10) /10) ;
		    b = ((end -endGe ) /10) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10;	 
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				adjoint_update(result_synthetic_odo, syntheticOdo10_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo10_dis_adjoint[tenSerial]);
		    }
		    //from the last  integel multiple of ten to end
		    g2o::Vector3 mid_vector3;
		    std::array<double,11> startNodeInfo;
		    for(int tenSerial = end - endGe; tenSerial < end; tenSerial++)
		    { 
		    	vitualStart = vitualStart + 1;  
		    	// cout<<"vitualStart: "<<vitualStart <<endl;	
		    	startNodeInfo = OdoInf[tenSerial];

				cov(0,0) = startNodeInfo[5]; 
				cov(0,1) = startNodeInfo[6];
				cov(1,0) = startNodeInfo[6];  
				cov(0,2) = startNodeInfo[7];
				cov(2,0) = startNodeInfo[7];
				cov(1,1) = startNodeInfo[8]; 
				cov(1,2) = startNodeInfo[9]; 
				cov(2,1) = startNodeInfo[9]; 
				cov(2,2) = startNodeInfo[10];

				mid_vector3[0] = startNodeInfo[2];
				mid_vector3[1] = startNodeInfo[3];
				mid_vector3[2] = startNodeInfo[4];

				next_pair.first.fromVector(mid_vector3);	
				next_pair.second = cov.inverse();

				distance_segment = sqrt(startNodeInfo[2]*startNodeInfo[2] + startNodeInfo[3]*startNodeInfo[3] );
				adjoint_update(result_synthetic_odo, next_pair, j1, j2, j3, midTrans, length, distance_segment);
		    }
		    if(vitualStart != end)
		    {
		    	cout<<"vitualStart: "<<vitualStart <<endl;
		    	cout<<"end: "<< end<<endl;
				cout<<"the final node does not equal to end vertex"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }

		    fultileBit = judgeVarianceAndDistance_Small(length, result_synthetic_odo.second);
		}
		//wan
		else
		{
			int increaseToTen ;
			if(startGe != 0)
				increaseToTen = start - startGe +10;
			else
				increaseToTen = start;
			synthesize_odo_edges_adjoint( start, increaseToTen, OdoInf, result_synthetic_odo.first, result_synthetic_odo.second, length, fultileBit);
			// int a = (end - endGe -(start - startGe +10)) /10;
			int end_ten = (end - endGe) /10 -1;

			//should be an integer multiple of ten
			if((end - endGe -(start - startGe +10)) % 10 != 0)
			{
				cout<<"the number of vectex is too big"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			vitualStart = increaseToTen;
			// cout<<"vitualStart: "<<vitualStart <<endl;
			//shi to bai
			a = increaseToTen / 10;
			b = (start -startGe - startShi*10 +100) /10 -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {   
		    	vitualStart = vitualStart + 10;		
		    	// cout<<"vitualStart: "<<vitualStart <<endl;						
				adjoint_update(result_synthetic_odo, syntheticOdo10_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo10_dis_adjoint[tenSerial]);
		    }
		    //bai to qian
		    a = (start -startGe - startShi*10 +100) /100;
		    b = ((start -startGe - startShi*10 - startBai*100 +1000 ) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {   
		    	vitualStart = vitualStart + 100;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				adjoint_update(result_synthetic_odo, syntheticOdo100_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo100_dis_adjoint[tenSerial]);
		    }
		    //qian to wan
		    a = (start - startGe - startShi*10 - startBai*100 +1000) /1000;
		    b = ((start - startGe - startShi*10 - startBai*100 - startQian*1000 + 10000) /1000) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1000;	 	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;							
				adjoint_update(result_synthetic_odo, syntheticOdo1000_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo1000_dis_adjoint[tenSerial]);
		    }
		    //wan to wan
		    a = (start - startGe - startShi*10 - startBai*100 - startQian*1000 +10000) /10000;
		    b = ((end - endGe - endShi*10 - endBai*100 - endQian*1000) /10000) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10000;	 
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				adjoint_update(result_synthetic_odo, syntheticOdo10000_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo10000_dis_adjoint[tenSerial]);
		    }
		    //wan to qian
		    a = (end - endGe - endShi*10 - endBai*100 - endQian*1000) /1000;
		    b = ((end -endGe - endShi*10 - endBai*100) /1000) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1000;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl; 								
				adjoint_update(result_synthetic_odo, syntheticOdo1000_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo1000_dis_adjoint[tenSerial]);
		    }
		    //qian to bai
		    a = ((end -endGe - endShi*10 - endBai*100) /100) ;
		    b = ((end -endGe - endShi*10 ) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 100;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl; 								
				adjoint_update(result_synthetic_odo, syntheticOdo100_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo100_dis_adjoint[tenSerial]);
		    }
		    //bai to shii
		    a = ((end -endGe - endShi*10) /10) ;
		    b = ((end -endGe ) /10) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl; 								
				adjoint_update(result_synthetic_odo, syntheticOdo10_trans_varMatrix_adjoint[tenSerial], j1, j2, j3, midTrans, length, 
					syntheticOdo10_dis_adjoint[tenSerial]);
		    }
		    //from the last  integel multiple of ten to end
		    g2o::Vector3 mid_vector3;
		    std::array<double,11> startNodeInfo;
		    for(int tenSerial = end - endGe; tenSerial < end; tenSerial++)
		    { 
		    	vitualStart = vitualStart + 1;
		    	// cout<<"vitualStart: "<<vitualStart <<endl;  	
		    	startNodeInfo = OdoInf[tenSerial];

				cov(0,0) = startNodeInfo[5]; 
				cov(0,1) = startNodeInfo[6];
				cov(1,0) = startNodeInfo[6];  
				cov(0,2) = startNodeInfo[7];
				cov(2,0) = startNodeInfo[7];
				cov(1,1) = startNodeInfo[8]; 
				cov(1,2) = startNodeInfo[9]; 
				cov(2,1) = startNodeInfo[9]; 
				cov(2,2) = startNodeInfo[10];

				mid_vector3[0] = startNodeInfo[2];
				mid_vector3[1] = startNodeInfo[3];
				mid_vector3[2] = startNodeInfo[4];

				next_pair.first.fromVector(mid_vector3);	
				next_pair.second = cov.inverse();

				distance_segment = sqrt(startNodeInfo[2]*startNodeInfo[2] + startNodeInfo[3]*startNodeInfo[3] );
				adjoint_update(result_synthetic_odo, next_pair, j1, j2, j3, midTrans, length, distance_segment);
		    }
		    if(vitualStart != end)
		    {
				cout<<"the final node does not equal to end vertex"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }
		    
		    fultileBit = judgeVarianceAndDistance_Small(length, result_synthetic_odo.second);
		}
		// cout<<"length[0]"<<length[0]<<endl;
		// cout<<"end: "<<end<<endl;
	}

	void accelerated_synthetic_odo(int start, int end, std::vector<std::array<double,11>> & OdoInf, 
		std::pair<g2o::SE2, Matrix3d> & result_synthetic_odo, std::array<double, 5>  & length,
		bool & fultileBit)
	{	
		g2o::SE2 midTrans, nextTrans; 
		Matrix3d nextVarMatrix, midVarMatrix, j1, j2,j3, cov;		
		int a ,b ,c;

		if (start > 99999 or end > 99999)
		{
			cout<<"start "<<start<<" end "<<end<<endl;
			cout<<"the number of vectex is too big"<<endl;
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		if(start > OdoInf.size() or end > OdoInf.size())
		{
			cout<<"start: "<<start<<" end:"<<end <<" OdoInf.size:"<<OdoInf.size() <<endl;
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}

		int startWan = 0, startQian= 0, startBai = 0, startShi = 0, startGe  = 0;
		int endWan   = 0,  endQian = 0,  endBai  = 0,  endShi  = 0,  endGe   = 0;

		get_digits(start, end, startWan, startQian, startBai, startShi, startGe, endWan, endQian, endBai, endShi, endGe);

		int vitualStart = start , virtualEnd =  end;	
		g2o::SE2 toFindOutWhyAcceleratedNotEqualSlow_SE2;
		Matrix3d toFindOutWhyAcceleratedNotEqualSlow_VAR;
		bool toFindOutWhyAcceleratedNotEqualSlow_fultileBit;
		std::array<double,5> toFindOutWhyAcceleratedNotEqualSlow_length;

		// cout<<"start: "<<start<<endl;
		//ge
		if ((start - startGe +20) > end) 
		    synthesize_odo_edges( start, end, OdoInf, result_synthetic_odo.first, result_synthetic_odo.second, length, fultileBit);
		//shi
		else if ((start - startGe- startShi*10 +200) > end) 
		{
			int increaseToTen ;
			if(startGe != 0)
				increaseToTen = start - startGe +10;
			else
				increaseToTen = start;
			vitualStart = increaseToTen;
			// cerr<<"vitualStart: "<<vitualStart <<endl;
			synthesize_odo_edges( start, increaseToTen, OdoInf, result_synthetic_odo.first, 
				result_synthetic_odo.second, length, fultileBit);	
			int decreaseToTen;
			decreaseToTen = (end - endGe)/10 -1;
			//should be an integer multiple of ten
			if((end - endGe -(start - startGe +10)) % 10 != 0)
			{
				cout<<"the number of vectex is too big"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}

			//from start to the nearest  integel multiple of ten
		    for(int tenSerial = (increaseToTen) / 10; tenSerial <= decreaseToTen; tenSerial++)
		    { 
		    	vitualStart = vitualStart + 10;  	
		    	// cerr<<"vitualStart: "<<vitualStart <<endl;
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //from the last  integel multiple of ten to end
		    g2o::Vector3 mid_vector3;
		    std::array<double,11> startNodeInfo;
		    for(int tenSerial = end - endGe; tenSerial < end; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1;   
		    	// cout<<"vitualStart: "<<vitualStart <<endl;
		    	startNodeInfo = OdoInf[tenSerial];

				cov(0,0) = startNodeInfo[5]; 
				cov(0,1) = startNodeInfo[6];
				cov(1,0) = startNodeInfo[6];  
				cov(0,2) = startNodeInfo[7];
				cov(2,0) = startNodeInfo[7];
				cov(1,1) = startNodeInfo[8]; 
				cov(1,2) = startNodeInfo[9]; 
				cov(2,1) = startNodeInfo[9]; 
				cov(2,2) = startNodeInfo[10];

				mid_vector3[0] = startNodeInfo[2];
				mid_vector3[1] = startNodeInfo[3];
				mid_vector3[2] = startNodeInfo[4];

				nextTrans.fromVector(mid_vector3);	
				nextVarMatrix = cov.inverse();
					
				Jacobian_4_edge_propagate(result_synthetic_odo.first, nextTrans, j1, j2);
				covariance_propagate(result_synthetic_odo.second, nextVarMatrix, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * nextTrans;
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + sqrt(startNodeInfo[2]*startNodeInfo[2] + startNodeInfo[3]*startNodeInfo[3] );
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    if(vitualStart != end)
		    {
				cout<<"the final node does not equal to end vertex"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }
		    fultileBit = judgeVarianceAndDistance_Small(length, midVarMatrix);
		    // fultileBit = judgeVarianceAndDistance_wholeMatrix(length, cov);
		}
		//bai wei 
		else if ((start - startGe- startShi*10 - startBai*100+2000) > end)
		{
			int increaseToTen ;
			if(startGe != 0)
				increaseToTen = start - startGe +10;
			else
				increaseToTen = start;
			//bai ge
			synthesize_odo_edges( start, increaseToTen, OdoInf, result_synthetic_odo.first, 
				result_synthetic_odo.second, length, fultileBit);

			vitualStart = increaseToTen; 
			// cout<<"vitualStart: "<<vitualStart <<endl;	

			//should be an integer multiple of ten
			if((end - endGe -(start - startGe +10)) % 10 != 0)
			{
				cout<<"the number of vectex is too big"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}

			//bai shi from start to the nearest  integel multiple of ten
			a = increaseToTen / 10;
			b = (start -startGe - startShi*10 +100) /10 -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10;  	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);

		    }
		    //bai to bai
		    a = (start -startGe - startShi*10 +100) /100;
		    b = ((end -endGe - endShi*10) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 100;  	

				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo100_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo100_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo100_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo100_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);

		    }
		    //bai to shi
		    a = (end - endGe - endShi*10) /10;
		    b = ((end - endGe) /10) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10;  	
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }

		    //from the last  integel multiple of ten to end
		    g2o::Vector3 mid_vector3;
		    std::array<double,11> startNodeInfo;
		    for(int tenSerial = end - endGe; tenSerial < end; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1; 
		    	startNodeInfo = OdoInf[tenSerial];

				cov(0,0) = startNodeInfo[5]; 
				cov(0,1) = startNodeInfo[6];
				cov(1,0) = startNodeInfo[6];  
				cov(0,2) = startNodeInfo[7];
				cov(2,0) = startNodeInfo[7];
				cov(1,1) = startNodeInfo[8]; 
				cov(1,2) = startNodeInfo[9]; 
				cov(2,1) = startNodeInfo[9]; 
				cov(2,2) = startNodeInfo[10];

				mid_vector3[0] = startNodeInfo[2];
				mid_vector3[1] = startNodeInfo[3];
				mid_vector3[2] = startNodeInfo[4];

				nextTrans.fromVector(mid_vector3);	
				nextVarMatrix = cov.inverse();
					
				Jacobian_4_edge_propagate(result_synthetic_odo.first, nextTrans, j1, j2);
				covariance_propagate(result_synthetic_odo.second, nextVarMatrix, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * nextTrans;
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + sqrt(startNodeInfo[2]*startNodeInfo[2] + startNodeInfo[3]*startNodeInfo[3] );
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);

		    }
		    if(vitualStart != end)
		    {
				cout<<"the final node does not equal to end vertex"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }
		    fultileBit = judgeVarianceAndDistance_Small(length, midVarMatrix);
		}
		//thousand 
		else if ((start - startGe- startShi*10 - startBai*100 - startQian*1000+20000) > end)
		{
			int increaseToTen ;
			if(startGe != 0)
				increaseToTen = start - startGe +10;
			else
				increaseToTen = start;
			synthesize_odo_edges( start, increaseToTen, OdoInf, result_synthetic_odo.first, result_synthetic_odo.second, length, fultileBit);
			vitualStart = increaseToTen;
			// cout<<"start: "<<start<<endl;
			// cout<<"ge shi bai qian: "<<startGe<<" "<<startShi<<" "<<startBai<<" "<<startQian<<endl;
			// cout<<"vitualStart: "<<vitualStart <<endl;

			//should be an integer multiple of ten
			if((end - endGe -(start - startGe +10)) % 10 != 0)
			{
				cout<<"the number of vectex is too big"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			// cout<<"result_synthetic_odo.second: "<< result_synthetic_odo.second<<endl;
			// synthesize_odo_edges( start, vitualStart, OdoInf, toFindOutWhyAcceleratedNotEqualSlow_SE2, 
			// 	toFindOutWhyAcceleratedNotEqualSlow_VAR, toFindOutWhyAcceleratedNotEqualSlow_length, 
			// 	toFindOutWhyAcceleratedNotEqualSlow_fultileBit);

			//shi to bai
			a = increaseToTen / 10;
			b = (start -startGe - startShi*10 +100) /10 -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {   
		    	vitualStart = vitualStart + 10;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;							
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //bai to qian
		    a = (start -startGe - startShi*10 +100) /100;
		    b = ((start -startGe - startShi*10 - startBai*100 +1000 ) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {   
		    	vitualStart = vitualStart + 100;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo100_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo100_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo100_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo100_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //qian to qian
		    a = (start -startGe - startShi*10 - startBai*100 +1000) /1000;
		    b = ((end -endGe - endShi*10 - endBai*100) /1000) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1000;	
		    	// cout<<"a: "<<a<<endl; 
		    	// cout<<"start: "<<start<<"startGe: "<<startGe<<"startShi: "<<startShi<<"start: "<<start<<"start: "<<start<<endl; 
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo1000_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo1000_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo1000_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo1000_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //qian to bai
		    a = ((end -endGe - endShi*10 - endBai*100) /100) ;
		    b = ((end -endGe - endShi*10 ) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 100;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl; 								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo100_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo100_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo100_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo100_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);

		    }
		    //bai to shi
		    a = ((end -endGe - endShi*10) /10) ;
		    b = ((end -endGe ) /10) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10;	 
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //from the last  integel multiple of ten to end
		    g2o::Vector3 mid_vector3;
		    std::array<double,11> startNodeInfo;
		    for(int tenSerial = end - endGe; tenSerial < end; tenSerial++)
		    { 
		    	vitualStart = vitualStart + 1;  
		    	// cout<<"vitualStart: "<<vitualStart <<endl;	
		    	startNodeInfo = OdoInf[tenSerial];

				cov(0,0) = startNodeInfo[5]; 
				cov(0,1) = startNodeInfo[6];
				cov(1,0) = startNodeInfo[6];  
				cov(0,2) = startNodeInfo[7];
				cov(2,0) = startNodeInfo[7];
				cov(1,1) = startNodeInfo[8]; 
				cov(1,2) = startNodeInfo[9]; 
				cov(2,1) = startNodeInfo[9]; 
				cov(2,2) = startNodeInfo[10];

				mid_vector3[0] = startNodeInfo[2];
				mid_vector3[1] = startNodeInfo[3];
				mid_vector3[2] = startNodeInfo[4];

				nextTrans.fromVector(mid_vector3);	
				nextVarMatrix = cov.inverse();
					
				Jacobian_4_edge_propagate(result_synthetic_odo.first, nextTrans, j1, j2);
				covariance_propagate(result_synthetic_odo.second, nextVarMatrix, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * nextTrans;
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + sqrt(startNodeInfo[2]*startNodeInfo[2] + startNodeInfo[3]*startNodeInfo[3] );
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    if(vitualStart != end)
		    {
		    	cout<<"vitualStart: "<<vitualStart <<endl;
		    	cout<<"end: "<< end<<endl;
				cout<<"the final node does not equal to end vertex"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }

		    fultileBit = judgeVarianceAndDistance_Small(length, midVarMatrix);
		}
		//wan
		else
		{
			int increaseToTen ;
			if(startGe != 0)
				increaseToTen = start - startGe +10;
			else
				increaseToTen = start;
			synthesize_odo_edges( start, increaseToTen, OdoInf, result_synthetic_odo.first, result_synthetic_odo.second, length, fultileBit);
			// int a = (end - endGe -(start - startGe +10)) /10;
			int end_ten = (end - endGe) /10 -1;

			//should be an integer multiple of ten
			if((end - endGe -(start - startGe +10)) % 10 != 0)
			{
				cout<<"the number of vectex is too big"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			vitualStart = increaseToTen;
			// cout<<"vitualStart: "<<vitualStart <<endl;
			//shi to bai
			a = increaseToTen / 10;
			b = (start -startGe - startShi*10 +100) /10 -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {   
		    	vitualStart = vitualStart + 10;		
		    	// cout<<"vitualStart: "<<vitualStart <<endl;						
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //bai to qian
		    a = (start -startGe - startShi*10 +100) /100;
		    b = ((start -startGe - startShi*10 - startBai*100 +1000 ) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {   
		    	vitualStart = vitualStart + 100;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo100_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo100_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo100_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo100_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		   
		    //qian to wan
		    a = (start - startGe - startShi*10 - startBai*100 +1000) /1000;
		    b = ((start - startGe - startShi*10 - startBai*100 - startQian*1000 + 10000) /1000) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1000;	 	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;							
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo1000_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo1000_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo1000_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo1000_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //wan to wan
		    a = (start - startGe - startShi*10 - startBai*100 - startQian*1000 +10000) /10000;
		    b = ((end - endGe - endShi*10 - endBai*100 - endQian*1000) /10000) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10000;	 
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10000_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10000_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10000_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10000_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }

		    //wan to qian
		    a = (end - endGe - endShi*10 - endBai*100 - endQian*1000) /1000;
		    b = ((end -endGe - endShi*10 - endBai*100) /1000) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1000;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl; 								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo1000_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo1000_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo1000_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo1000_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //qian to bai
		    a = ((end -endGe - endShi*10 - endBai*100) /100) ;
		    b = ((end -endGe - endShi*10 ) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 100;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl; 								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo100_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo100_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo100_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo100_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //bai to shii
		    a = ((end -endGe - endShi*10) /10) ;
		    b = ((end -endGe ) /10) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl; 								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //from the last  integel multiple of ten to end
		    g2o::Vector3 mid_vector3;
		    std::array<double,11> startNodeInfo;
		    for(int tenSerial = end - endGe; tenSerial < end; tenSerial++)
		    { 
		    	vitualStart = vitualStart + 1;
		    	// cout<<"vitualStart: "<<vitualStart <<endl;  	
		    	startNodeInfo = OdoInf[tenSerial];

				cov(0,0) = startNodeInfo[5]; 
				cov(0,1) = startNodeInfo[6];
				cov(1,0) = startNodeInfo[6];  
				cov(0,2) = startNodeInfo[7];
				cov(2,0) = startNodeInfo[7];
				cov(1,1) = startNodeInfo[8]; 
				cov(1,2) = startNodeInfo[9]; 
				cov(2,1) = startNodeInfo[9]; 
				cov(2,2) = startNodeInfo[10];

				mid_vector3[0] = startNodeInfo[2];
				mid_vector3[1] = startNodeInfo[3];
				mid_vector3[2] = startNodeInfo[4];

				nextTrans.fromVector(mid_vector3);	
				nextVarMatrix = cov.inverse();
					
				Jacobian_4_edge_propagate(result_synthetic_odo.first, nextTrans, j1, j2);
				covariance_propagate(result_synthetic_odo.second, nextVarMatrix, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * nextTrans;
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + sqrt(startNodeInfo[2]*startNodeInfo[2] + startNodeInfo[3]*startNodeInfo[3] );
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    if(vitualStart != end)
		    {
				cout<<"the final node does not equal to end vertex"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }
		    
		    fultileBit = judgeVarianceAndDistance_Small(length, midVarMatrix);
		}
		// cout<<"length[0]"<<length[0]<<endl;
		// cout<<"end: "<<end<<endl;
	}

	//input: the start and end node serial number correspond to the segment you want to synthesize
	//output: the tranform and covariance info of the synthesized segment, in Trans and cov
	void synthesize_odo_edges(int start, int end, std::vector<std::array<double,11>> & OdoInf, g2o::SE2 & Trans, 
		Matrix3d & cov, std::array<double, 5> & length, bool & fultileBit)
	{
		fultileBit = 0;
		std::array<double,11> startNodeInfo = OdoInf[start];
		g2o::Vector3  mid_vector3;
		Matrix3d m_m, m2, J1, J2, J_container, J_container_1, J_container_2;
		g2o::SE2 edge2;
		int lengthNode = end -start;
		// cout<<"start node is "<<start<<endl;
		// cout<<"end   node is "<<end<<endl;
		if(startNodeInfo[0] != start){
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		if(lengthNode < 0){
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		cov = Matrix3d::Identity(); 
		if(lengthNode == 0){
			cov(0,0)= 0; cov(1,1)= 0; cov(2,2) = 0;
			mid_vector3[0] = 0;mid_vector3[1] = 0;mid_vector3[2] = 0;
			Trans.fromVector(mid_vector3);
			length[0] = 0;
			length[1] = 0;
			length[2] = 0;

			length[3] = 0;
			length[4] = 0;			
		}
		else
		{
			if(lengthNode == 1 and (startNodeInfo[1] != end))
			{
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}			
			cov(0,0) = startNodeInfo[5]; 
			cov(0,1) = startNodeInfo[6];
			cov(1,0) = startNodeInfo[6];  
			cov(0,2) = startNodeInfo[7];
			cov(2,0) = startNodeInfo[7];
			cov(1,1) = startNodeInfo[8]; 
			cov(1,2) = startNodeInfo[9]; 
			cov(2,1) = startNodeInfo[9]; 
			cov(2,2) = startNodeInfo[10];

			Matrix3d cov_mid;
			cov_mid = cov.inverse();
			cov = cov_mid;

			mid_vector3[0] = startNodeInfo[2];
			mid_vector3[1] = startNodeInfo[3];
			mid_vector3[2] = startNodeInfo[4];
			Trans.fromVector(mid_vector3);
			// cout<<"node "<<start<<" "<<Trans[0]<<" "<<Trans[1]<<" "<<Trans[2]<<endl;
		

			length[0] = sqrt(Trans[0]*Trans[0]+Trans[1]*Trans[1]);
			length[1] = abs(Trans[0]);
			length[2] = abs(Trans[1]);

			length[3] = cov(0,0);
			length[4] = cov(1,1);

			int j=start, minus = start;
			for(j=start+1; j<end; j++)
			{
				m2 = Matrix3d::Identity();
				m2(0,0 ) = OdoInf[j][5];
				m2(0,1 ) = OdoInf[j][6];
				m2(0,2 ) = OdoInf[j][7];
				m2(1,0 ) = m2(0,1 );
				m2(2,0 ) = m2(0,2 );

				m2(1,1) = OdoInf[j][8]; 
				m2(1,2) = OdoInf[j][9]; 
				m2(2,1 ) = m2(1,2 );
				m2(2,2) = OdoInf[j][10];

				cov_mid = m2.inverse();
				m2 = cov_mid;

				mid_vector3[0] = OdoInf[j][2]; mid_vector3[1] = OdoInf[j][3]; mid_vector3[2] = OdoInf[j][4];
				edge2.fromVector(mid_vector3);

				// cout<<"node "<<j<<" "<<edge2[0]<<" "<<edge2[1]<<" "<<edge2[2]<<endl;

				length[0] = length[0] + sqrt(edge2[0]*edge2[0]+edge2[1]*edge2[1]);
				length[1] = length[1]+abs(edge2[0]);
				length[2] = length[2]+abs(edge2[1]);
				Jacobian_4_edge_propagate(Trans, edge2, J1, J2);
				covariance_propagate(cov, m2, J1, J2, m_m);

				// cout<<"**&********************   progress to "<<j<<"*********************"<<endl;
				// cout<<"J1"<<endl;
				// cout<<J1<<endl;

				// cout<<"J2"<<endl;
				// cout<<J2<<endl;

				// cout<<"cov"<<endl;
				// cout<<cov<<endl;

				// cout<<"m2"<<endl;
				// cout<<m2<<endl;

				// cout<<"j1 cov j1"<<endl;
				// cout<<J1*cov*(J1.transpose())<<endl;

				J_container = J1.transpose();
				J_container_1 = J1*cov*(J1.transpose());
				J_container_2 = J1*cov*(J_container);
				if(judgeTwoMatrix3dEqual(J_container_1, J_container_2) == 1)
				{
					cout<<"two matrix 3d not equal"<<endl;
					cout<<endl<<J_container_1<<endl;
					cout<<" "<<endl;
					cout<<J_container_2<<endl;
					cin.get();
				}

				// cout<<J1*cov*(J1.transpose())<<endl;

				// cout<<"j2 m2 j2"<<endl;
				// cout<<J2*m2*(J2.transpose())<<endl;

				// cout<<"m_m"<<endl;
				// cout<<m_m<<endl;	
				double norml2_cov = sqrt(m_m(0, 0)*m_m(0, 0)+m_m(1, 1)*m_m(1, 1));
				// cout<<"covx "<<m_m(0, 0) <<" covy "<<m_m(1, 1)<<" norm_cov "<<norml2_cov<<" length "<<  length[0]<<" angle_cov "<<m_m(2, 2) <<endl;	


				cov = m_m;

				length[3] = cov(0,0);
				length[4] = cov(1,1);

				double xa = Trans[0], ya = Trans[1], thetaa = Trans[2];
				double xb = edge2[0], yb = edge2[1], thetab = edge2[2];

				for(int i = 0; i < 3; i++)
				{
					if(abs(edge2[i] - mid_vector3[i]) > 0.00000001)
					{
						cout<<"wrong to construct edge form vector"<<endl;
						cout<<"vector:"<<endl<<mid_vector3[i]<<" edge[i]:"<<edge2[i]<<endl;
						printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);
					}
				}	

				Trans *= edge2;//update transform
				// cout<<"trans synthesis: "<<j<<" "<<Trans[0]<<" "<<Trans[1]<<" "<<Trans[2]<<endl;

				double xc = cos(thetaa)*xb - sin(thetaa)*yb + xa;
				double yc = sin(thetaa)*xb + cos(thetaa)*yb + ya;
				double thetac = thetaa + thetab;


			  	double const_pi =  3.14159265358979323846;
			  	if (thetac >= -const_pi && thetac < const_pi)
			  	{

			  	}
			  	else
			  	{
					 number_t multiplier = std::floor(thetac / (2*const_pi));
				  	thetac = thetac - multiplier*2*const_pi;
				  	if (thetac >= const_pi)
				    	thetac -= 2*const_pi;
				  	if (thetac < -const_pi)
				    	thetac += 2*const_pi;
			  	}


				if(abs(xc - Trans[0]) > 0.00000000001)
				{
					cout<<"xc: "<<xc<<" trans[0]:"<<Trans[0]<<endl;
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}

				if(abs(yc - Trans[1]) > 0.00000000001)
				{
					cout<<"yc: "<<yc<<" trans[0]:"<<Trans[1]<<endl;
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}

				if(abs(thetac - Trans[2]) > 0.00000000001)
				{
					cout<<"thetaa: "<<thetaa<<" thetab:"<<thetab<<endl;
					cout<<"thetac: "<<thetac<<" trans[0]:"<<Trans[2]<<endl;
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
				if(j - minus == 30)
				{
					fultileBit = judgeVarianceAndDistance_Small(length, cov);
					if(fultileBit == 1)
						return;
				}
			}

			fultileBit = judgeVarianceAndDistance_Small(length, cov);

			if(OdoInf[j-1][1] != end)
			{
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
		}
	}

	void synthesize_odo_edges_adjoint(int start, int end, std::vector<std::array<double,11>> & OdoInf, g2o::SE2 & Trans, 
		Matrix3d & cov, std::array<double, 5> & length, bool & fultileBit)
	{
		fultileBit = 0;
		std::array<double,11> startNodeInfo = OdoInf[start];
		g2o::Vector3  mid_vector3;
		Matrix3d m_m, m2, J1, J2, J_container, J_container_1, J_container_2, adjoint;
		g2o::SE2 edge2;
		int lengthNode = end -start;
		// cout<<"start node is "<<start<<endl;
		// cout<<"end   node is "<<end<<endl;
		if(startNodeInfo[0] != start){
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		if(lengthNode < 0){
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		cov = Matrix3d::Identity(); 
		if(lengthNode == 0){
			cov(0,0)= 0; cov(1,1)= 0; cov(2,2) = 0;
			mid_vector3[0] = 0;mid_vector3[1] = 0;mid_vector3[2] = 0;
			Trans.fromVector(mid_vector3);
			length[0] = 0;
			length[1] = 0;
			length[2] = 0;

			length[3] = 0;
			length[4] = 0;			
		}
		else
		{
			if(lengthNode == 1 and (startNodeInfo[1] != end))
			{
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}			
			cov(0,0) = startNodeInfo[5]; 
			cov(0,1) = startNodeInfo[6];
			cov(1,0) = startNodeInfo[6];  
			cov(0,2) = startNodeInfo[7];
			cov(2,0) = startNodeInfo[7];
			cov(1,1) = startNodeInfo[8]; 
			cov(1,2) = startNodeInfo[9]; 
			cov(2,1) = startNodeInfo[9]; 
			cov(2,2) = startNodeInfo[10];

			Matrix3d cov_mid;
			cov_mid = cov.inverse();
			cov = cov_mid;

			mid_vector3[0] = startNodeInfo[2];
			mid_vector3[1] = startNodeInfo[3];
			mid_vector3[2] = startNodeInfo[4];
			Trans.fromVector(mid_vector3);

			length[0] = sqrt(Trans[0]*Trans[0]+Trans[1]*Trans[1]);
			length[1] = abs(Trans[0]);
			length[2] = abs(Trans[1]);

			length[3] = cov(0,0);
			length[4] = cov(1,1);

			int j=start, minus = start;
			for(j=start+1; j<end; j++)
			{
				m2 = Matrix3d::Identity();
				m2(0,0 ) = OdoInf[j][5];
				m2(0,1 ) = OdoInf[j][6];
				m2(0,2 ) = OdoInf[j][7];
				m2(1,0 ) = m2(0,1 );
				m2(2,0 ) = m2(0,2 );

				m2(1,1) = OdoInf[j][8]; 
				m2(1,2) = OdoInf[j][9]; 
				m2(2,1 ) = m2(1,2 );
				m2(2,2) = OdoInf[j][10];

				cov_mid = m2.inverse();
				m2 = cov_mid;

				mid_vector3[0] = OdoInf[j][2]; mid_vector3[1] = OdoInf[j][3]; mid_vector3[2] = OdoInf[j][4];
				edge2.fromVector(mid_vector3);
				length[0] = length[0] + sqrt(edge2[0]*edge2[0]+edge2[1]*edge2[1]);
				length[1] = length[1]+abs(edge2[0]);
				length[2] = length[2]+abs(edge2[1]);

				get_adjoint(adjoint, Trans);
				J2 = adjoint.transpose();
				m_m = cov + adjoint*m2*(J2) ;
				cov = m_m;
				// Jacobian_4_edge_propagate(Trans, edge2, J1, J2);
				// covariance_propagate(cov, m2, J1, J2, m_m);

				// cout<<"**&********************   progress to "<<j<<"*********************"<<endl;
				// cout<<"J1"<<endl;
				// cout<<J1<<endl;

				// cout<<"J2"<<endl;
				// cout<<J2<<endl;

				// cout<<"cov"<<endl;
				// cout<<cov<<endl;
				// cout<<"m2"<<endl;
				// cout<<m2<<endl;

				// cout<<"j1 cov j1"<<endl;
				// cout<<J1*cov*(J1.transpose())<<endl;

				// J_container = J1.transpose();
				// J_container_1 = J1*cov*(J1.transpose());
				// J_container_2 = J1*cov*(J_container);
				// if(judgeTwoMatrix3dEqual(J_container_1, J_container_2) == 1)
				// {
				// 	cout<<"two matrix 3d not equal"<<endl;
				// 	cout<<endl<<J_container_1<<endl;
				// 	cout<<" "<<endl;
				// 	cout<<J_container_2<<endl;
				// 	cin.get();
				// }

				// cout<<J1*cov*(J1.transpose())<<endl;

				// cout<<"j2 m2 j2"<<endl;
				// cout<<J2*m2*(J2.transpose())<<endl;

				// cout<<"m_m"<<endl;
				// cout<<m_m<<endl;	

				// m_m = cov.transpose();
				// cout<<"m_m: "<<endl<<m_m<<endl;	
				// (cov*m_m).trace
				// double norml2_cov = sqrt(m_m(0, 0)*m_m(0, 0)+m_m(1, 1)*m_m(1, 1));
				// cout<<"covx "<<m_m(0, 0) <<" covy "<<m_m(1, 1)<<" norm_cov "<<norml2_cov<<" length "<<  length[0]<<" angle_cov "<<m_m(2, 2) <<endl;	

				length[3] = cov(0,0);
				length[4] = cov(1,1);

	

				Trans *= edge2;//update transform

			}
			fultileBit = judgeVarianceAndDistance_Small(length, cov);

			if(OdoInf[j-1][1] != end)
			{
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
		}
	}

	void get_adjoint(Matrix3d & adj, const g2o::SE2 & trans)
	{
		adj.block(0,0,2,2) = trans.rotation().toRotationMatrix();
		// cout<<"trans: "<<trans[0]<<" "<<trans[1]<<" "<<trans[2]<<endl;
		// cout<<"adj.block(0,0,2,2)"<<endl<<adj.block(0,0,2,2)<<endl;
		adj(2,0) = 0;
		adj(2,1) = 0;
		adj(0,2) = trans[1];
		adj(1,2) = -trans[0];
		adj(2,2) = 1;
		// cout<<"adj"<<endl<<adj<<endl;
		// cin.get();

	}
	void Jacobian_4_edge_propagate(g2o::SE2 & TransA, g2o::SE2 & TransB, Matrix3d & J1, Matrix3d & J2)
	{
		J1 = Matrix3d::Identity();
		J2 = Matrix3d::Identity();
		J1(0,2) = -TransB[0]*sin(TransA[2])-TransB[1]*cos(TransA[2]);
		J1(1,2) = TransB[0]*cos(TransA[2])-TransB[1]*sin(TransA[2]);

		
		J2(0,0) = cos(TransA[2]);
		J2(0,1) = -sin(TransA[2]);
		J2(1,0) = sin(TransA[2]);
		J2(1,1) = cos(TransA[2]);
	}

	void covariance_propagate(Matrix3d & cov1, Matrix3d & cov2, Matrix3d & J1, Matrix3d & J2, Matrix3d & result)
	{
		result = J1*cov1*(J1.transpose())+J2*cov2*(J2.transpose());
	}

	std::pair<bool, double> check_single_loop(int loop_to_check, cluster & _clustersFoundi, int loop_checked,
		std::vector<std::vector<std::pair<g2o::SE2, Matrix3d> > > &transSequence_cluster, int clus_num, int group_num)//, double& statis
	{
		g2o::SE2 loop1, edgeGo, loop2, edgeBack, Edge_midd, Edge_go_midd, Edge_back_midd, transform_interator;
		Matrix3d Cov_loop1, Cov_edgeGo, Cov_loop2, Cov_edgeBack, Cov_midd, Cov_go_midd, Cov_back_midd, J1, J2, Cov_interator;
		MatrixXd T(1,3),T_inverse(3,1); 
		std::pair<bool, double> returnV;

		// std::array<int,8> nodes;
		int go_odo_num = 0, back_odo_num = 0;
		// //nodes loop check
		// nodes[0] = _clustersFoundi.positionserial[loop_checked][0];
		// nodes[1] = _clustersFoundi.positionserial[loop_checked][3];
		if(loop_to_check <= loop_checked)
		{
			cout<<"loop_to_check: "<<loop_to_check<<" loop_checked: "<<loop_checked<<endl;
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}

		loop1        =  transSequence_cluster[2][loop_checked].first;
		Cov_loop1    =  transSequence_cluster[2][loop_checked].second;
		loop2        =  transSequence_cluster[2][loop_to_check].first;
		Cov_loop2    =  transSequence_cluster[2][loop_to_check].second;
		if(loop_to_check == loop_checked+1)
		{ 
			edgeGo       =  transSequence_cluster[1][loop_checked].first;
			Cov_edgeGo  =  transSequence_cluster[1][loop_checked].second;
			edgeBack     =  transSequence_cluster[0][loop_checked].first;
			Cov_edgeBack =  transSequence_cluster[0][loop_checked].second;	
		}
		else
		{
			edgeGo       =  transSequence_cluster[1][loop_checked].first;
			Cov_edgeGo  =  transSequence_cluster[1][loop_checked].second;
			edgeBack     =  transSequence_cluster[0][loop_checked].first;
			Cov_edgeBack =  transSequence_cluster[0][loop_checked].second;	

			for(int sumodo = 1; sumodo < loop_to_check; sumodo++)
			{

				Matrix3d m2 = Matrix3d::Identity(), m_m, J1_2_odo, J2_2_odo;

				Edge_go_midd  =  transSequence_cluster[1][loop_checked+sumodo].first;
				Cov_go_midd   =  transSequence_cluster[1][loop_checked+sumodo].second;	
				Edge_back_midd  =  transSequence_cluster[0][loop_checked+sumodo].first;
				Cov_back_midd   =  transSequence_cluster[0][loop_checked+sumodo].second;	

				Jacobian_4_edge_propagate(edgeGo, Edge_go_midd, J1, J2);
				covariance_propagate(Cov_edgeGo, Cov_go_midd, J1, J2, m_m);
				Cov_edgeGo = m_m;
				edgeGo *= Edge_go_midd;//update transform

				Jacobian_4_edge_propagate(edgeBack, Edge_back_midd, J1, J2);
				covariance_propagate(Cov_edgeBack, Cov_back_midd, J1, J2, m_m);
				Cov_edgeBack = m_m;
				edgeBack *= Edge_go_midd;//update transform

				if(loop_checked+sumodo+1 == loop_to_check)
					break;
				else if(sumodo == loop_to_check-1)
				{
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
			}
		}
		Jacobian_4_edge_propagate(loop1, edgeGo, J1, J2);//generate jacobian 
		covariance_propagate(Cov_loop1, Cov_edgeGo, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		Cov_interator = Cov_midd;
		transform_interator = loop1 * edgeGo;//update transform

		Edge_midd = loop2.inverse();
		Jacobian_4_edge_propagate(transform_interator, Edge_midd, J1, J2);//generate jacobian 
		covariance_propagate(Cov_interator, Cov_loop2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		Cov_interator = Cov_midd;
		transform_interator *=  Edge_midd;//update transform

		Edge_midd = edgeBack.inverse();
		Jacobian_4_edge_propagate(transform_interator, Edge_midd, J1, J2);//generate jacobian 
		covariance_propagate(Cov_interator, Cov_edgeBack, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		Cov_interator = Cov_midd;
		transform_interator *=  Edge_midd;//update transform

		// T = transform_interator.toVector();
		Matrix3d mmmm =  Cov_interator.inverse();
		T(0)= transform_interator[0];
		T(1) = transform_interator[1];
		T(2) = transform_interator[2];	
		T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		T_inverse(1) = T(1) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		T_inverse(2) = T(2) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		double transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);

		// double transformDistance = (T *) * T_inverse;
		// double transformDistance = T(0)*T(0)/Cov_interator(0,0) +  T(1)*T(1)/Cov_interator(1,1) +  T(2)*T(2)/Cov_interator(2,2);
		cout<<"cluster num: "<<clus_num<<" group_num: "<<group_num<<endl<<" loop1: "<<loop_checked<<" loop2: "<<loop_to_check<<endl;
		cout<<"cov: "<<endl<<Cov_interator<<endl;
		cout<<"transformDistance: "<<transformDistance<<endl;
		cout<<"transform_interator: "<<transform_interator[0]<<" "<<transform_interator[1]<<" "<<transform_interator[2]<<endl;
		cout<<" "<<endl;
		returnV.second = transformDistance;
		if (transformDistance < 7.81)
			returnV.first = 1;
		else
			returnV.first = 0;
		return returnV;
	}	

//input is four pair<transform, covariance> of two loop and two odo edge segments
//return is one pair<pass_check_or_not, transfrom_distance>
	std::pair<bool, double> check_single_loop_inter(std::array<std::pair<g2o::SE2, Matrix3d>, 4 > &transSequence_cluster_inter, 
	double& covX, double & covY, Matrix3d & displayCov,double & transX_residual, double & transY_residual, double & transA_residual,
	double & length, bool & futileb)//, double& statis
	{
		g2o::SE2 loop1, loop2, Edge_midd;
		Matrix3d Cov1, Cov2, Cov_midd, J1, J2;
		MatrixXd T(1,3),T_inverse(3,1); 
		std::pair<bool, double> returnV;

		loop1        =  transSequence_cluster_inter[0].first;
		Cov1    =  transSequence_cluster_inter[0].second;

		for(int i = 1; i<4; i++)
		{
			loop2   =  transSequence_cluster_inter[i].first;
			Cov2    =  transSequence_cluster_inter[i].second;

			Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
			covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
			Cov1 = Cov_midd;
			// //dispaly the covariacne matrix
			// cout<<Cov2.inverse()<<endl;
			// cout<<" "<<endl;

			loop1 = loop1 * loop2;//update transform
		}
		// exit(0);

		// double SNR = (length*length)/sqrt(Cov1(0, 0)*Cov1(0, 0) + Cov1(1, 1)*Cov1(1, 1));
		double SNR = (length)/sqrt(Cov1(0, 0)*Cov1(0, 0) + Cov1(1, 1)*Cov1(1, 1));

		// cout<<"length[0]: "<<length[0]<<" length[3]: "<<length[3]<<" length[4]: "<<length[4]<<endl;
		// cout<<"sqrt(length[3]*length[3] + length[4]*length[4]): "<<sqrt(length[3]*length[3] + length[4]*length[4])<<endl;
		if(SNR < snrThres)//3.16 corresponding to 5 dB
		// if(SNR < 1)
		// if(SNR < 10)
		{
			futileb = 1;
		}
		else
			futileb = 0;


		displayCov = Cov1;
		covX = Cov1(0, 0);
		covY = Cov1(1, 1);
		// T = transform_interator.toVector();
		Matrix3d mmmm =  Cov1.inverse();
		T(0)= loop1[0];
		T(1) = loop1[1];
		T(2) = loop1[2];	
		T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		T_inverse(1) = T(0) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		T_inverse(2) = T(0) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		double transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);
		transX_residual = T(0);
		transY_residual = T(1);
		transA_residual = T(2);
		// cout<<"transX_residual: "<<transX_residual<<"transY_residual: "<<transY_residual<<"transA_residual: "<<transA_residual<<endl;
		// cout<<"transformDistance: "<<transformDistance<<endl;
		// exit(0);

		// double transformDistance = (T *) * T_inverse;
		// double transformDistance = T(0)*T(0)/Cov_interator(0,0) +  T(1)*T(1)/Cov_interator(1,1) +  T(2)*T(2)/Cov_interator(2,2);

		// cout<<"cluster num: "<<clus_num<<" group_num: "<<group_num<<endl<<" loop1: "<<loop_checked<<" loop2: "<<loop_to_check<<endl;
		// cout<<"cov: "<<endl<<Cov_interator<<endl;
		// cout<<"transformDistance: "<<transformDistance<<endl;
		// cout<<"transform_interator: "<<transform_interator[0]<<" "<<transform_interator[1]<<" "<<transform_interator[2]<<endl;
		// cout<<" "<<endl;
		returnV.second = transformDistance;
		if (abs(transformDistance) < utils::chi2(3, 0.7))//11.34 correspond to P=0.01, which 7.81 to 0.05
			returnV.first = 1;
		else
			returnV.first = 0;
		return returnV;

	}
	std::pair<bool, double> check_single_loop_inter_varying_belief(std::array<std::pair<g2o::SE2, Matrix3d>, 4 > &transSequence_cluster_inter, 
	double& covX, double & covY, Matrix3d & displayCov,double & transX_residual, double & transY_residual, double & transA_residual,
	double & length, bool & futileb, double & belief)//, double& statis
	{
		g2o::SE2 loop1, loop2, Edge_midd, loop_two_odo;
		Matrix3d Cov1, Cov2, Cov_midd, J1, J2;
		MatrixXd T(1,3),T_inverse(3,1); 
		std::pair<bool, double> returnV;


		loop1   =  transSequence_cluster_inter[0].first;
		Cov1    =  transSequence_cluster_inter[0].second;


		// // check two odo segments futile state
		// loop2   =  (transSequence_cluster_inter[2].first).inverse();
		// Cov2    =  transSequence_cluster_inter[2].second;

		// Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
		// covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		// Cov_two_odo = Cov_midd;

		// loop_two_odo = loop1 * loop2;//update transform

		// cout<<"0: "<<loop1[0]<<" "<<loop1[1]<<" "<<loop1[2]<<endl;
		for(int i = 1; i<4; i++)
		{
			loop2   =  transSequence_cluster_inter[i].first;
			Cov2    =  transSequence_cluster_inter[i].second;

			Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
			covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
			Cov1 = Cov_midd;
			// //dispaly the covariacne matrix
			// cout<<Cov2.inverse()<<endl;
			// cout<<" "<<endl;

			loop1 = loop1 * loop2;//update transform
			// cout<<i<<": "<<loop1[0]<<" "<<loop1[1]<<" "<<loop1[2]<<endl;
			// cout<< "    "<<loop2[0]<<" "<<loop2[1]<<" "<<loop2[2]<<endl;
		}
		// exit(0);

		// // double SNR = (length*length)/sqrt(Cov1(0, 0)*Cov1(0, 0) + Cov1(1, 1)*Cov1(1, 1));
		// double SNR = sqrt(Cov1(0, 0)*Cov1(0, 0) + Cov1(1, 1)*Cov1(1, 1))/(length);

		// // cout<<"length[0]: "<<length[0]<<" length[3]: "<<length[3]<<" length[4]: "<<length[4]<<endl;
		// // cout<<"sqrt(length[3]*length[3] + length[4]*length[4]): "<<sqrt(length[3]*length[3] + length[4]*length[4])<<endl;
		// if(SNR > snrThres)//3.16 corresponding to 5 dB
		// // if(SNR < 1)
		// // if(SNR < 10)
		// {
		// 	futileb = 1;
		// }
		// else
		// 	futileb = 0;


		displayCov = Cov1;
		covX = Cov1(0, 0);
		covY = Cov1(1, 1);
		// T = transform_interator.toVector();
		LOG(INFO)<<"COV"<<Cov1;
		Matrix3d mmmm =  Cov1.inverse();
		LOG(INFO)<<"COV"<<Cov1;
		LOG(INFO)<<"COV"<<mmmm;
		T(0)= loop1[0];
		T(1) = loop1[1];
		T(2) = loop1[2];	

		LOG(INFO)<<"T "<<T;
		LOG(INFO)<<"COV of T: "<<Cov1;
		// LOG(INFO)<<"T transpose "<<T.transpose();
		// LOG(INFO)<<"T "<<T;
		// exit(0);
		// length = T[0]*T[0] + T[1]*T[1];
		futileb = judgeVarianceAndDistance_Small(length, Cov1);

		T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		T_inverse(1) = T(0) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		T_inverse(2) = T(0) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		double transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);
		LOG(INFO)<<"transformDistance: "<<transformDistance;
		LOG(INFO)<<"TRANS dis with transpose :"<<T*mmmm*(T.transpose());

		// loop1   =  transSequence_cluster_inter[1].first;
		// Cov1    =  transSequence_cluster_inter[1].second;
		// loop2   =  transSequence_cluster_inter[2].first;
		// Cov2    =  transSequence_cluster_inter[2].second;
		// Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
		// covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		// Cov1 = Cov_midd;
		// loop1 = loop1 * loop2;//update transform

		// loop2   =  transSequence_cluster_inter[3].first;
		// Cov2    =  transSequence_cluster_inter[3].second;
		// Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
		// covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		// Cov1 = Cov_midd;
		// loop1 = loop1 * loop2;//update transform

		// loop2   =  transSequence_cluster_inter[0].first;
		// Cov2    =  transSequence_cluster_inter[0].second;
		// Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
		// covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		// Cov1 = Cov_midd;
		// loop1 = loop1 * loop2;//update transform

		// mmmm =  Cov1.inverse();
		// T(0)= loop1[0];
		// T(1) = loop1[1];
		// T(2) = loop1[2];
		// T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		// T_inverse(1) = T(0) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		// T_inverse(2) = T(0) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		//  transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);
		//  LOG(INFO)<<"T "<<T;
		// LOG(INFO)<<"transformDistance: "<<transformDistance;
		// LOG(INFO)<<"TRANS dis with transpose :"<<T*mmmm*(T.transpose());


		// loop1   =  transSequence_cluster_inter[2].first;
		// Cov1    =  transSequence_cluster_inter[2].second;
		// loop2   =  transSequence_cluster_inter[3].first;
		// Cov2    =  transSequence_cluster_inter[3].second;
		// Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
		// covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		// Cov1 = Cov_midd;
		// loop1 = loop1 * loop2;//update transform

		// loop2   =  transSequence_cluster_inter[0].first;
		// Cov2    =  transSequence_cluster_inter[0].second;
		// Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
		// covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		// Cov1 = Cov_midd;
		// loop1 = loop1 * loop2;//update transform

		// loop2   =  transSequence_cluster_inter[1].first;
		// Cov2    =  transSequence_cluster_inter[1].second;
		// Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
		// covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		// Cov1 = Cov_midd;
		// loop1 = loop1 * loop2;//update transform

		// mmmm =  Cov1.inverse();
		// T(0)= loop1[0];
		// T(1) = loop1[1];
		// T(2) = loop1[2];
		// T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		// T_inverse(1) = T(0) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		// T_inverse(2) = T(0) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		//  transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);
		//  LOG(INFO)<<"T "<<T;
		// LOG(INFO)<<"transformDistance: "<<transformDistance;
		// LOG(INFO)<<"TRANS dis with transpose :"<<T*mmmm*(T.transpose());



		// loop1   =  transSequence_cluster_inter[3].first;
		// Cov1    =  transSequence_cluster_inter[3].second;
		// loop2   =  transSequence_cluster_inter[0].first;
		// Cov2    =  transSequence_cluster_inter[0].second;
		// Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
		// covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		// Cov1 = Cov_midd;
		// loop1 = loop1 * loop2;//update transform

		// loop2   =  transSequence_cluster_inter[1].first;
		// Cov2    =  transSequence_cluster_inter[1].second;
		// Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
		// covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		// Cov1 = Cov_midd;
		// loop1 = loop1 * loop2;//update transform

		// loop2   =  transSequence_cluster_inter[2].first;
		// Cov2    =  transSequence_cluster_inter[2].second;
		// Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
		// covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		// Cov1 = Cov_midd;
		// loop1 = loop1 * loop2;//update transform

		// mmmm =  Cov1.inverse();
		// T(0)= loop1[0];
		// T(1) = loop1[1];
		// T(2) = loop1[2];
		// T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		// T_inverse(1) = T(0) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		// T_inverse(2) = T(0) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		// transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);
		// LOG(INFO)<<"T "<<T;
		// LOG(INFO)<<"transformDistance: "<<transformDistance;
		// LOG(INFO)<<"TRANS dis with transpose :"<<T*mmmm*(T.transpose());


		// exit(0);
		transX_residual = T(0);
		transY_residual = T(1);
		transA_residual = T(2);
		double transAngle = T_inverse(2)*T(2);
		double transTrans = T_inverse(0)*T(0) + T_inverse(1)*T(1);
		// cout<<"transX_residual: "<<T_inverse(0)*T(0)<<"transY_residual: "<<T_inverse(1)*T(1)<<"transA_residual: "<<T_inverse(2)*T(2)<<endl;

		// cout<<"loop1[0]: "<<loop1[0]<<" loop1[1]:"<<loop1[1]<<" loop1[2]:"<<loop1[2]<<endl;
		// cout<<Cov1<<endl;
		// exit(0);

		// double transformDistance = (T *) * T_inverse;
		// double transformDistance = T(0)*T(0)/Cov_interator(0,0) +  T(1)*T(1)/Cov_interator(1,1) +  T(2)*T(2)/Cov_interator(2,2);

		// cout<<"cluster num: "<<clus_num<<" group_num: "<<group_num<<endl<<" loop1: "<<loop_checked<<" loop2: "<<loop_to_check<<endl;
		// cout<<"cov: "<<endl<<Cov_interator<<endl;
		// cout<<"transformDistance: "<<transformDistance<<endl;
		// cout<<"transform_interator: "<<transform_interator[0]<<" "<<transform_interator[1]<<" "<<transform_interator[2]<<endl;
		// cout<<" "<<endl;
		returnV.second = transformDistance;
		if (abs(transformDistance) < utils::chi2_continuous(3, belief))
		// if (abs(transAngle) < utils::chi2_continuous(1, belief))//11.34 correspond to P=0.01, which 7.81 to 0.05
		// if((transformDistance < utils::chi2_continuous(3, belief)) and (transAngle < utils::chi2_continuous(1, belief)))
		// if((transTrans < utils::chi2_continuous(2, belief)) and (transAngle < utils::chi2_continuous(1, belief)))
			returnV.first = 1;
		else
			returnV.first = 0;
		return returnV;

	}

	std::pair<bool, double> check_single_loop_odo(const std::pair<int, int> & loop_node_pair, 
		 std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	& LP_Trans_Covar_Map, 
		bool & fultileBit1, const double & belief)//, double& statis
	{
		//input variable: std::pair<int, int>  and full info of four segments.
		int start, end;
		double length;
		g2o::SE2  Trans1, loop_two_odo;
		std::pair<g2o::SE2, Matrix3d > result_synthetic_odo;
		Matrix3d  cov1, J1, J2, midCov;
		std::array<double ,5> length_;

		start = loop_node_pair.second;
		end   = (loop_node_pair).first;
		if(start > end)
		{
			// synthesize_odo_edges( end,  start,  OdoInf, Trans1,  cov1, length[0], fultileBit);///
			accelerated_synthetic_odo( end, start, OdoInf, 
				result_synthetic_odo,  length_, 
				fultileBit1);
			Trans1 = result_synthetic_odo.first.inverse();
			cov1   = result_synthetic_odo.second;
		}
		else
		{
			// synthesize_odo_edges( start,  end,  OdoInf, Trans1,  cov1, length[0], fultileBit);
			accelerated_synthetic_odo(start, end, OdoInf, 
				result_synthetic_odo,  length_, 
				fultileBit1);
			Trans1 = result_synthetic_odo.first;
			cov1   = result_synthetic_odo.second;
		}
		length = (LP_Trans_Covar_Map[loop_node_pair].first[0])*(LP_Trans_Covar_Map[loop_node_pair].first[0]) + 
				(LP_Trans_Covar_Map[loop_node_pair].first[1])*(LP_Trans_Covar_Map[loop_node_pair].first[1]) +
				length_[0];

		MatrixXd T(1,3),T_inverse(3,1); 
		std::pair<bool, double> returnV;


		Jacobian_4_edge_propagate(Trans1, LP_Trans_Covar_Map[(loop_node_pair)].first, J1, J2);//generate jacobian 
		covariance_propagate(cov1, LP_Trans_Covar_Map[(loop_node_pair)].second, J1, J2, midCov);// update covariance from two covs and two Jacobians

		LOG(INFO)<<"Synthesis trans "<<Trans1[0]<<" "<<Trans1[1]<<" "<<Trans1[2];
		LOG(INFO)<<"single loop check, loop "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" is "<<LP_Trans_Covar_Map[(loop_node_pair)].first[0]<<" "<<
			LP_Trans_Covar_Map[(loop_node_pair)].first[1]<<" "<<LP_Trans_Covar_Map[(loop_node_pair)].first[2];
		LOG(INFO)<<"VAR MATRIX:"<<endl<<cov1;
		loop_two_odo = Trans1 * LP_Trans_Covar_Map[(loop_node_pair)].first;//update transform

		fultileBit1 = judgeVarianceAndDistance_Small(length, midCov);

		displayCov = midCov;

		// T = transform_interator.toVector();
		Matrix3d mmmm =  midCov.inverse();
		T(0)= loop_two_odo[0];
		T(1) = loop_two_odo[1];
		T(2) = loop_two_odo[2];	
		T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		T_inverse(1) = T(0) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		T_inverse(2) = T(0) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		double transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);
		double transAngle = T_inverse(2)*T(2);
		double transTrans = T_inverse(0)*T(0) + T_inverse(1)*T(1);
		// cout<<"transX_residual: "<<transX_residual<<"transY_residual: "<<transY_residual<<"transA_residual: "<<transA_residual<<endl;
		// cout<<"transformDistance: "<<transformDistance<<" length:"<<length<<endl;
		// cout<<"cov: "<<endl;
		// cout<<Cov1<<endl;
		// exit(0);

		// double transformDistance = (T *) * T_inverse;
		// double transformDistance = T(0)*T(0)/Cov_interator(0,0) +  T(1)*T(1)/Cov_interator(1,1) +  T(2)*T(2)/Cov_interator(2,2);

		// cout<<"cluster num: "<<clus_num<<" group_num: "<<group_num<<endl<<" loop1: "<<loop_checked<<" loop2: "<<loop_to_check<<endl;
		// cout<<"cov: "<<endl<<Cov_interator<<endl;
		// cout<<"transformDistance: "<<transformDistance<<endl;

		// cout<<" "<<endl;
		// cout<<"loop trans: "<<LP_Trans_Covar_Map[(loop_node_pair)].first[0]<<" "<<LP_Trans_Covar_Map[(loop_node_pair)].first[1]<<
		// 	" "<<LP_Trans_Covar_Map[(loop_node_pair)].first[2]<<endl;
		// cout<<"Trans1: "<<Trans1[0]<<" "<<Trans1[1]<<" "<<Trans1[2]<<endl;
		// cout<<"Trans1.inverse: "<<Trans1.inverse()[0]<<" "<<Trans1.inverse()[1]<<" "<<Trans1.inverse()[2]<<endl;
		// cout<<"loop transform:"<<loop_two_odo[0]<<" "<<loop_two_odo[1]<<" "<<loop_two_odo[2]<<endl;
		// cout<<"covariance: "<<endl<<mmmm<<endl;
		// cout<<" "<<endl;
		returnV.second = transformDistance;
		if (abs(transformDistance) < utils::chi2_continuous(3, belief))
		// if (abs(transAngle) < utils::chi2_continuous(1, belief))//11.34 correspond to P=0.01, which 7.81 to 0.05
		// if((transformDistance < utils::chi2_continuous(3, belief)) and (transAngle < utils::chi2_continuous(1, belief)))
		// if((transTrans < utils::chi2_continuous(2, belief)) and (transAngle < utils::chi2_continuous(1, belief)))
			returnV.first = 1;
		else
			returnV.first = 0;
		return returnV;
	}

	std::pair<bool, double> check_single_loop_inter_varying_belief_adjoint(std::array<std::pair<g2o::SE2, Matrix3d>, 4 > &transSequence_cluster_inter, 
	double& covX, double & covY, Matrix3d & displayCov,double & transX_residual, double & transY_residual, double & transA_residual,
	double & length, bool & futileb, double & belief)//, double& statis
	{
		g2o::SE2 midTrans;
		Matrix3d j1, j2, j3;
		MatrixXd T(1,3),T_inverse(3,1); 
		std::pair<bool, double> returnV;
		std::pair<g2o::SE2, Matrix3d> result_synthetic_odo;
		double distance_segment;
		std::array<double, 5> length_vector;

		result_synthetic_odo.first     =  transSequence_cluster_inter[0].first;
		result_synthetic_odo.second    =  transSequence_cluster_inter[0].second;

		// check two odo segments futile state
		// cout<<"first cov"<<endl<<result_synthetic_odo.second<<endl;
		for(int i = 1; i<4; i++)
		{
			adjoint_update(result_synthetic_odo, transSequence_cluster_inter[i], j1, j2, j3, midTrans, length_vector, distance_segment);
			// cout<<i<<"ed cov:"<<transSequence_cluster_inter[i].second<<endl;
		}

		futileb = judgeVarianceAndDistance_Small(length, result_synthetic_odo.second);

		displayCov = result_synthetic_odo.second;
		covX = result_synthetic_odo.second(0, 0);
		covY = result_synthetic_odo.second(1, 1);
		// T = transform_interator.toVector();
		Matrix3d mmmm =  result_synthetic_odo.second.inverse();
		T(0)= result_synthetic_odo.first[0];
		T(1) = result_synthetic_odo.first[1];
		T(2) = result_synthetic_odo.first[2];	
		T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		T_inverse(1) = T(0) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		T_inverse(2) = T(0) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		double transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);

		if(std::isnan(transformDistance))
		{
			cout<<"the trans dis of adjoint is nan, the information matrix:"<<endl<<mmmm<<endl<<"T 0 1 2 : "<<
				result_synthetic_odo.first[0]<<" "<<result_synthetic_odo.first[1]<<" "<<result_synthetic_odo.first[0]<<endl;
			cin.get();
		}
		transX_residual = T(0);
		transY_residual = T(1);
		transA_residual = T(2);
		double transAngle = T_inverse(2)*T(2);
		double transTrans = T_inverse(0)*T(0) + T_inverse(1)*T(1);

		returnV.second = transformDistance;
		if (abs(transformDistance) < utils::chi2_continuous(3, belief))
		// if (abs(transAngle) < utils::chi2_continuous(1, belief))//11.34 correspond to P=0.01, which 7.81 to 0.05
		// if((transformDistance < utils::chi2_continuous(3, belief)) and (transAngle < utils::chi2_continuous(1, belief)))
		// if((transTrans < utils::chi2_continuous(2, belief)) and (transAngle < utils::chi2_continuous(1, belief)))
			returnV.first = 1;
		else
			returnV.first = 0;
		return returnV;

	}

	IntPairSet& getClusterByID(int id){
		// cout<<"go into the get cluster by id program"<<endl;
		if(clusterIDtoLoopsMap.find(id) == clusterIDtoLoopsMap.end())
		{
			cout<<"no such a cluster id exist"<<endl;
		}
		// else
		// 	cout<<"the cluster corresponding to this id exists"<<endl;
		return clusterIDtoLoopsMap[id];
	}

	int getClusterSize(int id){
		// cout<<"go into the get cluster by id program"<<endl;
		if(clusterIDtoLoopsMap.find(id) == clusterIDtoLoopsMap.end())
		{
			cout<<"no such a cluster id exist"<<endl;
			return -1;
		}

		int size = clusterIDtoLoopsMap[id].size();
		return size;
	}

	size_t clusterCount()
	{
		return _clustersFound.size();
	}

	bool deleteCluster(int clusterID)
	{
		clusterIDtoLoopsMap.erase(clusterID);

		for(IntPairIDMap::iterator it= loopToClusterIDMap.begin();
				it!=loopToClusterIDMap.end(); it++)
		{
			if(it->second == clusterID)
				loopToClusterIDMap.erase(it->first);
		}

		return true;
	}
	void merge_cluster(std::vector<std::vector<int> > & consistent_pair_clusterr_real)
	{
	}
	inline void get_digits(int & start, int & end, int & startWan, int & startQian, int & startBai, int & startShi, int & startGe,
						   int & endWan,  int & endQian,  int & endBai,  int & endShi,  int & endGe)
	{
		if(start > 9999)
		{
			startWan =  start / 10000;          startQian= (start % 10000) /1000;
			startBai = (start % 1000) /100;     startShi = (start % 100) /10;       startGe  =  start % 10;
		}
		else if(start > 999)
		{
			startQian= (start % 10000) /1000;   startBai = (start % 1000) /100;
			startShi = (start % 100) /10;       startGe  =  start % 10;
		}
		else if(start > 99)
		{
			startBai = (start % 1000) /100;     startShi = (start % 100) /10;      startGe  =  start % 10;
		}
		else if(start > 9)
		{
			startShi =  (start % 100) /10;      startGe  =  start % 10;
		}
		else
			startGe  =  start;			


		if(end > 9999)
		{
			endWan =  end / 10000;         endQian= (end % 10000) /1000;
			endBai = (end % 1000) /100;    endShi = (end % 100) /10;       endGe  =  end % 10;
		}
		else if(end > 999)
		{
			endQian= (end % 10000) /1000;  endBai = (end % 1000) /100;
			endShi = (end % 100) /10;      endGe  =  end % 10;
		}
		else if(end > 99)
		{
			endBai = (end % 1000) /100;    endShi = (end % 100) /10;       endGe  =  end % 10;
		}
		else if(end > 9)
		{
			endShi =  (end % 100) /10;     endGe  =  end % 10;
		}
		else
			endGe  =  end;		
	}
    inline void adjoint_update(std::pair<g2o::SE2, Matrix3d> & result_synthetic_odo, std::pair<g2o::SE2, Matrix3d> & incremental,
    	Matrix3d & adjoint, Matrix3d & adjoint_transpose, Matrix3d & middle_cov, g2o::SE2 & middle_se2,
    	std::array<double, 5>  & length, double & dis)
    {
       	get_adjoint(adjoint, result_synthetic_odo.first);
       	// cout<<"adjoint:"<<endl<<adjoint<<endl;
       	adjoint_transpose = adjoint.transpose();
		middle_cov = result_synthetic_odo.second + adjoint*(incremental.second)*adjoint_transpose;
		result_synthetic_odo.second = middle_cov;
		middle_se2 = result_synthetic_odo.first*(incremental.first);
		result_synthetic_odo.first = middle_se2;
		length[0] = length[0] + dis;
		length[3] = result_synthetic_odo.second(0,0);
		length[4] = result_synthetic_odo.second(1,1);
    };	
		
};

#endif /* CLUSTER_HPP_ */
