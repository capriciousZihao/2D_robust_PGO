// RRR - Robust Loop Closing over Time
// Copyright (C) 2014 Y.Latif, C.Cadena, J.Neira


#ifndef CLUSTER_HPP_
#define CLUSTER_HPP_

#include <iostream>
#include <fstream> 
#include <string.h>
#include <vector>
#include <cstdlib>
#include <Eigen/Dense>
#include <math.h>
#include "utils.hpp"
using namespace Eigen; 
using namespace std;
using namespace utils;

struct cluster
{
	int startLow, startHigh;
	int endLow, endHigh;
	int size;
	std::vector<std::array<double,6>> positionserial;
	std::vector<double> dis_cluster_backup;
	std::vector<double> dis_cluster_start_end;

	cluster(): startLow(-1), startHigh(-1), endLow(-1), endHigh(-1), size(0) {}
	cluster(int start, int end) : startLow(start), startHigh(start), endLow(end), endHigh(end), size(1){}

	bool contains(int start, int end, int threshold)
	{
		return
				(
					std::abs(start-startHigh) < threshold or	std::abs(start-startLow)  < threshold
				)
				and
				(
					std::abs(end-endHigh) < threshold or std::abs(end-endLow) < threshold
				);

	}
	// std::vector<std::array<double,4>> VertexInf
	std::array<double,2> cal_distance(std::array<double,6> LoopPosition)
	{
		
		double startDIS,endDIS;
		int IDofNearestLC=0, IDDynamic=0;
		std::array<double,6> lastLI;
		std::array<double,2> ret;
		if(positionserial.size()==0)
			cout<<"can't calculate distance because positionserial has no information"<<endl;

		for(std::vector<std::array<double,6>>::const_iterator it = positionserial.begin(),
			lendVertex = positionserial.end();it!=lendVertex;it++)
		{
			lastLI = *it;
			startDIS = abs(lastLI[1]-LoopPosition[1])+abs(lastLI[2]-LoopPosition[2]);
			endDIS = abs(lastLI[4]-LoopPosition[4])+abs(lastLI[5]-LoopPosition[5]);
			if(it == positionserial.begin())
			{
				// ret[0] = startDIS;
				ret[0] = IDofNearestLC;
 				ret[1] = endDIS+startDIS;
			}
			else if (endDIS+startDIS < ret[1])
			{
				// ret[0] = startDIS;
				ret[0] = IDofNearestLC;
				ret[1] = endDIS+startDIS;
			}
			IDofNearestLC++;
		}
		return ret;
	}

	void update_distance()
	{
		std::array<double,6> LoopPosition = *(positionserial.end()-1);
		std::array<double,6> lastLI;
		double startDIS,endDIS;
		std::array<double,2> ret;
		if(positionserial.size()==0)
			cout<<"can't calculate distance because positionserial has no information"<<endl;
		//clear the backup dis vector
		dis_cluster_backup.clear();


		if(dis_cluster_start_end.size()!=positionserial.size()-1)
		{
			cout<<"update distance abnormal: dis size does not match positionserial size"<<endl;
			exit(0);
		}

		std::vector<double>::const_iterator pointer_dis_vector = dis_cluster_start_end.begin();

		for(std::vector<std::array<double,6>>::const_iterator it = positionserial.begin()+1,
			lendVertex = positionserial.end()-1;it!=lendVertex;it++)
		{
			lastLI = *it;
			startDIS = abs(lastLI[1]-LoopPosition[1])+abs(lastLI[2]-LoopPosition[2]);
			endDIS = abs(lastLI[4]-LoopPosition[4])+abs(lastLI[5]-LoopPosition[5]);

			if (endDIS+startDIS < *pointer_dis_vector)
			{
				dis_cluster_backup.push_back(endDIS+startDIS);
			}
			else
			{
				dis_cluster_backup.push_back(*pointer_dis_vector);
			}
			pointer_dis_vector++;
		}

		dis_cluster_backup.push_back(*pointer_dis_vector);

		if(dis_cluster_start_end.size()!=dis_cluster_backup.size())
		{
			cout<<"update distance abnormal: dis_cluster_start_end size != dis_cluster_backup size"<<endl;
			exit(0);
		}
	
	}
};

class Clusterizer
{
	typedef IntPairIDMap 			LoopToClusterIDMap;//loop is a pair of Vertex ID, this map is a int pair,LC, to a int, cluster ID 
	typedef IDintPairSetMap 		ClusterIDtoLoopsMap;
	

	
	std::vector<std::array<double,4>> VertexInf;
	std::vector<double> distance;
	std::vector<cluster> _clustersFound;
	ClusterIDtoLoopsMap	clusterIDtoLoopsMap;
	LoopToClusterIDMap  loopToClusterIDMap;

public:
	string     nameofclusterfile;
	// Assumes that in the vector the values are passed as (start_1,end_1), (start_2,end_2), ...

	int getClusterID(const IntPair& loop)//input a pair, lC, output the ID of the cluster to which the pair belongs
	{
		return loopToClusterIDMap[loop];
	}

	void setClusterID(const IntPair& loopClosure, int ID)
	{
		if(loopToClusterIDMap.find(loopClosure)!=loopToClusterIDMap.end())
		{
			int oldID = loopToClusterIDMap[loopClosure];

//			for(IntPairSet::iterator
//					it = clusterIDtoLoopsMap[oldID].begin(),
//					end = clusterIDtoLoopsMap[oldID].end();
//					it!=end;
//					it++)
//
//			{
				//if (*it == loopClosure){
					clusterIDtoLoopsMap[oldID].erase(loopClosure);
				//	loopToClusterIDMap[*it] = ID;
				//	break;
				//}
//			}
			clusterIDtoLoopsMap[ID].insert(loopClosure);
			loopToClusterIDMap[loopClosure] = ID;

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

	std::array<double,2> find_nearest_cluster(std::array<double,6> fullLoopInfo, std::vector<cluster> _clustersFound)
	{
		double dis = 0, ID = 0;
		std::array<double,2> returndis;

		for(size_t i=0; i<_clustersFound.size(); i++)
		{
			returndis=_clustersFound[i].cal_distance(fullLoopInfo);
			if(i == 0)
			{
				dis = returndis[1];
				ID = 0;
			}
			
			else if(returndis[1] < dis)
			{
				dis = returndis[1];
				ID = i;
			}
		}
		returndis[0] = ID;
		returndis[1] = dis;
		return returndis;
	}
		
	// std::array<double,2> chi2_test(std::vector<double> dis_clu)
	std::array<double,2> chi2_test(cluster _clustersFound)
	{
		// std::array<double,2> p_value;

		// double sum = std::accumulate(std::begin(dis_clu), std::end(dis_clu), 0.0);  
		// double mean =  sum / dis_clu.size(); //均值  
						  
		// double accum  = 0.0;  
		// std::for_each (std::begin(dis_clu), std::end(dis_clu), [&](const double d) {  
		// 	accum  += (d-mean)*(d-mean);  
		// 	});  
						  
		// double chi_statis = (accum/mean); //卡方统计量
		// p_value[0] = 1-chi2_p(dis_clu.size()-1, chi_statis);	
		std::vector<double> dis_clu = _clustersFound.dis_cluster_start_end;
		std::vector<double> dis_clu_back = _clustersFound.dis_cluster_backup;
		std::array<double,2> p_value;

		double sum = std::accumulate(std::begin(dis_clu_back), std::end(dis_clu_back), 0.0);  
		double mean =  sum / dis_clu_back.size(); //均值  
						  
		double accum  = 0.0;  
		std::for_each (std::begin(dis_clu_back), std::end(dis_clu_back), [&](const double d) {  
			accum  += (d-mean)*(d-mean);  
			});  
						  
		double chi_statis = (accum/mean); //卡方统计量
		p_value[0] = 1-chi2_p(dis_clu.size()-1, chi_statis);	
		if(dis_clu.size()>2)
		{
			dis_clu.pop_back();
			sum = std::accumulate(std::begin(dis_clu), std::end(dis_clu), 0.0);  
			mean =  sum / dis_clu.size(); //均值  
			accum  = 0.0;  
			std::for_each (std::begin(dis_clu), std::end(dis_clu), [&](const double d) {  
			accum  += (d-mean)*(d-mean);  
			});  
						  
			 chi_statis = (accum/mean); //卡方统计量
			p_value[1] = 1-chi2_p(dis_clu.size()-1, chi_statis);	
		}
		else
		{
			p_value[1] = p_value[0];
		}
		cout<<"p_value:"<<p_value[0]<<"  pre_p_value:"<<p_value[1]<<
			"  dof:"<<(dis_clu.size()-1)<<"  chi_statis:"<<
			chi_statis<<endl;
		return p_value;
	}		

	void clusterize_zihao( const IntPairSet& loops, const char* filename, double thres)//, std::vector<int>& membership, int& clusterCount)
	{
		double disSTART, disEND,  pre_p_value = 0.5;
		int num_loop=0;
		bool findStartVertexInMatrix=false,findEndVertexInMatrix=false, ready2chiTest=false;
		std::array<double,3> startPosition, endPosition;
		std::array<double,2> tem_dis, nearest_cluster, p_value;
		std::array<double,6> fullLoopInfo;//six elements:start ID,X,Y,end ID,X,Y		

		collect_vertex(filename);

		if(loops.empty())
		{
			std::cerr<<"clusterize(): "<<__LINE__<<" no loops to make clusters"<<std::endl;
			return;
		}
		_clustersFound.clear();
		for(IntPairSet::const_iterator it = loops.begin(), lend = loops.end();it!=lend;it++)
		{
			int start 	= std::max(it->first,it->second);
			int end 	= std::min(it->first,it->second);

			//print loop number and cluster id of the loop
			num_loop++;
			cout<<"loop "<<num_loop<<" "<<start<<" "<<end<<endl;

			fullLoopInfo = get_LC_Pos( start,  end);//get loop closure vextexes position


			if(_clustersFound.empty())
			{
				cluster s(start,end);
				s.positionserial.push_back(fullLoopInfo);

				_clustersFound.push_back(s);

			// cout<<"fullLoopInfo: "<<fullLoopInfo[0]<<" "<<fullLoopInfo[1]<<" "<<fullLoopInfo[2]<<" "<<fullLoopInfo[3]<<
			// " "<<fullLoopInfo[4]<<" "<<fullLoopInfo[5]<<endl;

						cout<<"fullLoopInfo: "<<s.positionserial[0][0]<<" "<<s.positionserial[0][1]<<" "
						<<s.positionserial[0][2]<<" "<<s.positionserial[0][3]<<
			" "<<s.positionserial[0][4]<<" "<<s.positionserial[0][5]<<endl;


				clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
				loopToClusterIDMap[*it] = _clustersFound.size()-1;
			}
			else
			{
				cluster* currentCluster = NULL;

				//search for the nearest cluster to the loop
				nearest_cluster = find_nearest_cluster(fullLoopInfo, _clustersFound);
				
				//use the backup dis vector update_distance dis_cluster_backup
				//_clustersFound[nearest_cluster[0]].dis_cluster_backup.

				//add dis information, ad fullinfo of loop closure, update map of ID and LC
				_clustersFound[nearest_cluster[0]].dis_cluster_start_end.push_back(nearest_cluster[1]);
				_clustersFound[nearest_cluster[0]].positionserial.push_back(fullLoopInfo);
				

				std::array<double,6> tryeveryposserial;
				tryeveryposserial=_clustersFound[nearest_cluster[0]].positionserial.back();
				
				cout<<"tryeveryposserial: "<<tryeveryposserial[0]<<" "<<tryeveryposserial[1]<<" "
						<<tryeveryposserial[2]<<" "<<tryeveryposserial[3]<<
						" "<<tryeveryposserial[4]<<" "<<tryeveryposserial[5]<<endl;


				//if the size of dis_cluster is bigger than one, we do chi2 test
				if(_clustersFound[nearest_cluster[0]].dis_cluster_start_end.size()==1)
				{
					clusterIDtoLoopsMap[nearest_cluster[0]].insert(*it);
					loopToClusterIDMap[*it] = nearest_cluster[0];

					continue;
				}
				else
				{
					_clustersFound[nearest_cluster[0]].update_distance();
					 
					 // p_value = chi2_test(_clustersFound[nearest_cluster[0]].dis_cluster_start_end);
					 p_value = chi2_test(_clustersFound[nearest_cluster[0]]);
					 cout<<"thres: "<<thres<<endl;
					 if((p_value[0] > 0.05) and (abs(p_value[0] - p_value[1]) < thres))//
					 {

					 	clusterIDtoLoopsMap[nearest_cluster[0]].insert(*it);
						loopToClusterIDMap[*it] = nearest_cluster[0];
						// update dis: move elements in dis_backup to dis

					 	continue;
					 }
					 else
					 {
					 	int num = _clustersFound[nearest_cluster[0]].dis_cluster_start_end.size();
					 	double dis1 = *(_clustersFound[nearest_cluster[0]].dis_cluster_start_end.end()-1);
					 	double dis2 = *(_clustersFound[nearest_cluster[0]].dis_cluster_start_end.end()-2);
					 	if(num > 2 or ((num == 2) and (dis2 < dis1)))
					 	{
					 		_clustersFound[nearest_cluster[0]].dis_cluster_start_end.pop_back();
					 		_clustersFound[nearest_cluster[0]].positionserial.pop_back();

					 		cluster s(start,end);
							s.positionserial.push_back(fullLoopInfo);

							_clustersFound.push_back(s);
							clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
							loopToClusterIDMap[*it] = _clustersFound.size()-1;
					 	}
					 	else if((num == 2) and (dis2 > dis1))
					 	{
					 		//create new cluster
					 		cluster s(start,end);
					 		//add distance to new cluster
							s.dis_cluster_start_end.push_back(nearest_cluster[1]);
							//add posiinfo to new cluster
							s.positionserial.push_back(*(_clustersFound[nearest_cluster[0]].positionserial.end()-2));
							s.positionserial.push_back(*(_clustersFound[nearest_cluster[0]].positionserial.end()-1));
							//add new cluster to _clustersFound
							_clustersFound.push_back(s);
							//delete old position infor and dis info from old cluster
							cout<<"_clustersFound[nearest_cluster[0]].positionserial.size: "<<_clustersFound[nearest_cluster[0]].positionserial.size()<<endl;
							cout<<"_clustersFound[nearest_cluster[0]].dis_cluster_start_end.size: "<<
								_clustersFound[nearest_cluster[0]].dis_cluster_start_end.size()<<endl;
							//set LC and ID map
							clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
							loopToClusterIDMap[*it] = _clustersFound.size()-1;
							//handle the ID and LC map of another loop 
							std::pair<int,int> tochange;
							int startID = (*(_clustersFound[nearest_cluster[0]].positionserial.end()-2))[0];
							int endID = (*(_clustersFound[nearest_cluster[0]].positionserial.end()-2))[3];
							tochange.first = startID;
							tochange.second = endID;
							//delete the last two elements in positonserial and two all elements in dis_serial
							_clustersFound[nearest_cluster[0]].positionserial.pop_back();
							_clustersFound[nearest_cluster[0]].positionserial.pop_back();
							_clustersFound[nearest_cluster[0]].dis_cluster_start_end.pop_back();
							//debug
							std::array<double,6> mammamlai;
							mammamlai = (*(_clustersFound[nearest_cluster[0]].positionserial.end()-2));
							cout<<mammamlai[0]<<" "<<mammamlai[1]<<" "<<mammamlai[2]<<endl;
							cout<<"nearest_cluster:"<<nearest_cluster[0]<<" "<<nearest_cluster[1]<<endl;
							cout<<"higher id:"<<startID<<" "<<endID<<endl;

							if(clusterIDtoLoopsMap[nearest_cluster[0]].count(tochange))
							{
								clusterIDtoLoopsMap[nearest_cluster[0]].erase(tochange);
								clusterIDtoLoopsMap[_clustersFound.size()-1].insert(tochange);
								//need to chage: delete map from loop to cluster ID
								loopToClusterIDMap[tochange] = _clustersFound.size()-1;
								cout<<"loop "<< tochange.first<<" and "<<tochange.second<<"to cluster id map: "<<
									loopToClusterIDMap[tochange]<<endl;
							}
							else
							{
								tochange.first = endID;
								tochange.second = startID;


								
								cout<<"else  loop "<< tochange.first<<" and "<<tochange.second<<"to cluster id map: "<<
									loopToClusterIDMap[tochange]<<endl;
								// std::for_each (std::begin(clusterIDtoLoopsMap[0]), std::end(clusterIDtoLoopsMap[0]), 
								// 	[&](const std::pair<int int> d) {  
									 
								// 	});  

								// for(int i=0;i<clusterIDtoLoopsMap[0].size();i++)
								// {
								// 	std::pair<int int> d;
								// 	d= *clusterIDtoLoopsMap[0][i];
								// 	cout<<"lc: "<d.first<<" "<<d.second<<endl; 
								// }
 
							    for(std::set<std::pair<int,int >>::iterator it = clusterIDtoLoopsMap[0].begin();
							    	it!=clusterIDtoLoopsMap[0].end();it++)  
							    {  
			
							            cout << (*it).first <<" "<< (*it).second << "}\n";  
							    } 

							
								if(!clusterIDtoLoopsMap[nearest_cluster[0]].count(tochange))
								{
									cout<<"can't find the previous loop in mapset"<<endl;
									cout<<"id:"<<startID<<" "<<endID<<endl;
									cout<<clusterIDtoLoopsMap[nearest_cluster[0]].size()<<endl;
									cout<<"total cluster size:"<<_clustersFound.size()<<endl;
									// cout<<(*(_clustersFound[nearest_cluster[0]].positionserial.end()-2))[0]<<" "<<
									//    (*(_clustersFound[nearest_cluster[0]].positionserial.end()-2))[1]<<" "
									// (*(_clustersFound[nearest_cluster[0]].positionserial.end()-2))[2]<<" "<<
									// (*(_clustersFound[nearest_cluster[0]].positionserial.end()-2))[3]<<" "<<
									// (*(_clustersFound[nearest_cluster[0]].positionserial.end()-2))[4]<<" "<<
									// (*(_clustersFound[nearest_cluster[0]].positionserial.end()-2))[5]<<" "<<endl;
									exit(0);
								}
								clusterIDtoLoopsMap[nearest_cluster[0]].erase(tochange);
								clusterIDtoLoopsMap[_clustersFound.size()-1].insert(tochange);
								loopToClusterIDMap[tochange] = _clustersFound.size()-1;

							}
						}
							// _clustersFound[nearest_cluster[0]].dis_cluster_start_end.size();
					 }
					 
				}
			}



		}

		ofstream fileStream;  
		// fileStream.open("clusterFile.g2o",ios::trunc);
		fileStream.open(nameofclusterfile,ios::trunc);

		std::array<double,6> ty={1,1,1,1,1,1};


		for(size_t i=0; i< _clustersFound.size(); i++)
		{
			for(std::vector<std::array<double,6>>::const_iterator itVertex = _clustersFound[i].positionserial.begin(), 
				lendVertex = _clustersFound[i].positionserial.end();itVertex!=lendVertex;itVertex++)
			{
				ty = *itVertex;

				fileStream<<i<<"  "<<ty[0]<<"  "<<ty[1]<<"  "<<ty[2]<<"  "<<ty[3]<<"  "<<ty[4]<<"  "<<ty[5]<<"\n";
				// fileStream<<trystdarray[0]<<"\n";
			}	
		}
		fileStream.close();
		cout<<"origianl cluster file has been saved"<<endl;

		
	

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

// 	/* collect vertex ID and position and angle in vertexInfo
	
// 	 */

	int collect_vertex(const char* filename)
	{
		ifstream fileStream;  
	    string tmp,temp1;  
	    std::array<double,4> verT={0, 0, 0, 0};
	    char* seg;
	    int count = 0,dddlndgsfdgj=0;// 行数计数器  
	    int position = 0;  
        if (position != string::npos)  


	    // fileStream.open("B25b_0.500.g2o",ios::in);//ios::in 表示以只读的方式读取文件
	    fileStream.open(filename,ios::in);  
	    if(fileStream.fail())//文件打开失败:返回0  
	    {  
	        return 0;  
	    }  
	    else//文件存在  
	    {  
		    while(getline(fileStream,tmp,'\n'))//读取一行  
		    {  
		    	position = tmp.find("VERTEX_SE2");
		    	// cout<<tmp<<endl; 	
		    	// exit(0);

		    	if (tmp.size() > 0 and position!=string::npos)  
		    	{
		    		 // printf(typeid(position));
		    		// std::cout<<typeid(position)<<std::endl;
                    // cout<<tmp<<endl;
		    		// printf("%s",tmp);
		    		stringstream stringin(tmp);
		    		for (dddlndgsfdgj=0;dddlndgsfdgj<5;dddlndgsfdgj++)
		    		{
		    			switch(dddlndgsfdgj)
		    			{
	                    	case 0:
		    					stringin >> temp1;
		    					// cout<<temp1<<endl;
		    					break;
		    				case 1:
		    					stringin >> verT[0];
		    					break;
		    				case 2:
		    					stringin >> verT[1];
		    					break;
		    				case 3:
		    					stringin >> verT[2];
		    					break;
		    				case 4:
		    					stringin >> verT[3];
		    					break;
		    				default:
		    					break;
		    			}
		    		}
		    		count++; 
		    		VertexInf.push_back(verT);
		    		// cout<<verT[0]<<" "<<verT[1]<<" "<<verT[2]<<" "<<verT[3]<<endl;
		    		// cout<<VertexInf.size()<<endl;
		    		// cout<<VertexInf[0][0]<<" "<<VertexInf[0][1]<<" "<<VertexInf[0][2]<<" "<<VertexInf[0][3]<<endl;
		    		// exit(0);
		    	}
		    }  
		    printf("number of vertexes:%d \n",count);
		    cout<<"VertexInf.size():"<<VertexInf.size()<<endl;
		    fileStream.close();  
		    // exit(0);
		    return count ;  
		}  
	}


// std::ifstream infile_feat(loadFeatList); //加载数据文件
// 	std::string feature; //存储读取的每行数据
// 	float feat_onePoint;  //存储每行按空格分开的每一个float数据
// 	std::vector<float> lines; //存储每行数据
// 	std::vector<vector<float>> lines_feat; //存储所有数据
// 	lines.clear();
// 	lines_feat.clear();

// 	while(!infile_feat.eof()) 
// 	{	
// 		getline(infile_feat, feature); //一次读取一行数据
// 		stringstream stringin(feature); //使用串流实现对string的输入输出操作
// 		while (stringin >> feat_onePoint) {  //按空格一次读取一个数据存入feat_onePoint
//       			lines.push_back(feat_onePoint);  //存储每行按空格分开的数据
//    		}
// 		lines_feat.push_back(lines);  //存储所有数据
// 	}

// 	infile_feat.close();



	IntPairSet& getClusterByID(int id){
		return clusterIDtoLoopsMap[id];
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

};

#endif /* CLUSTER_HPP_ */
