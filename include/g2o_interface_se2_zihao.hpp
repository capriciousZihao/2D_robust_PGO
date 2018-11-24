#ifndef G2O_INTERFACE_SE2_ZIHAO_HPP_
#define G2O_INTERFACE_SE2_ZIHAO_HPP_

#include"iostream"
#include"ctime"
#include <algorithm> 
#include <Eigen/Dense>
#include "g2o/types/slam2d/se2.h"

using namespace Eigen;
using namespace std; 



static bool cmp(  const std::pair<double, int> p,  const std::pair<double, int> q)
{
	return p.first < q.first;
}

static bool cmp_1(  const  std::pair<double, std::pair<std::pair<int, int>, int> > p,  const std::pair<double, std::pair<std::pair<int, int>, int> > q)
{
	return p.first < q.first;
}	
static bool cmp_vertor_int_pair(  const  std::pair<int, std::vector<int> > p,  const std::pair<int, std::vector<int> > q)
{
	return p.first > q.first;
}

void rand_of_n(std::vector<int> & a,int n);  //产生 1-n 的随机排列并存到 a[] 中 


void rand_of_n(std::vector<int> & a,int n)
{
	int i;
	std::pair<double, int> dis_ele;
	std::vector<std::pair<double, int> > sequence_dis_interClusters;
	
	srand((int)time(0));  // 初始化随机数种子 
	for(i=0;i<n;i++){
		dis_ele.first=rand();  // 随机生成一个数 
		dis_ele.second=i;
		sequence_dis_interClusters.push_back(dis_ele);
	}

	std::sort(sequence_dis_interClusters.begin(), sequence_dis_interClusters.end(), cmp);  //排序 
	for(i=0;i<n;i++){
		a[i]= sequence_dis_interClusters[i].second;
	}
}

bool judgeTwoMatrix3dEqual(Matrix3d & j1, Matrix3d & j2)
{
	for(int i =0; i<3; i++)
	{
		for(int p =0; p<3; p++)
		{
			if ( abs( ( j1(i,p) - j2(i,p) )/j1(i,p) )>0.01)
			// if ( abs( ( j1(i,p) - j2(i,p) ) )>0.00000000001)
			{ 
				cout<<"abs(j1(i,p) - j2(i,p)): "<<abs(j1(i,p) - j2(i,p))<<endl;
				cout<<"j1(i,p): "<<j1(i,p)<<endl;
				cout<<"j2(i,p): "<<j2(i,p)<<endl;
				return 1;
			}
		}
	}
	return 0;
}
bool judgeTwoSE2Equal(g2o::SE2 & SE1, g2o::SE2 & SE2)
{
	for(int i =0; i<3; i++)
	{
		// if ( abs( ( j1(i,p) - j2(i,p) )/j1(i,p) )>0.00000000001)
		if ( abs( ( SE1[i] - SE2[i] ) )>0.00000000001)
		{ 
			cout<<SE1[0] <<" "<<SE1[1] <<" "<<SE1[2]<<endl;
			cout<<SE2[0] <<" "<<SE2[1] <<" "<<SE2[2]<<endl;
			return 1;
		}
	}
	return 0;
}


class find_element{
   	public:
    int find(std::vector<int > & vector, int obj)
	{
		for(int i = 0; i < vector.size(); i++)
		{
			if(vector[i] == obj)
				return i;
		}
		return -1;
	}

    int find(std::set<int > & vector, int obj)
	{
		auto start = vector.begin();
		int j = 0;
		for(; start != vector.end(); start++)
		{
			if(*start == obj)
			{
				return j;
			}
			else
				j = j + 1;
		}
		return -1;
	}

    int find_SecondBiggest(std::vector<double > & vector, double & obj)
	{
		double biggest = *(std::max_element(vector.begin(), vector.end()));
		int serial;
		obj = 0.0;
		for(int start = 0; start < vector.size(); start++)
		{
			if(vector[start] < 0)
			{
				printf("the value is negative. This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			if(vector[start] < biggest and vector[start] > obj)
			{
				serial = start;
				obj = vector[start];
			}
		}
		return serial;
	}


	int find(std::vector<double > & vector, double obj)
	{
		for(int i = 0; i < vector.size(); i++)
		{
			if(vector[i] == obj)
				return i;
		}
		return -1;
	}

	int find(std::vector<std::vector<int > > & vector, std::pair<int,int > & p) 
	{
		std::pair<int, int> p1, p2 ;
		for(int i = 0; i < vector.size(); i++)
		{
			if(vector[i].size() != 2)
			{
				cerr<<"this function compare a pair to a vector<vector<int> >, howerer, the size of vector<int> is not two"<<endl;
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			p1.first  = vector[i][0];
			p1.second = vector[i][1];
			p2.first  = vector[i][1];
			p2.second = vector[i][0];		

			if(p1 == p or p2== p)
				return i;
		}
		return -1;
	}

	int find(std::vector<std::vector<int > > & vector, std::vector<int > & p) 
	{
		for(int i = 0; i < vector.size(); i++)
		{
			if(vector[i].size() != p.size())
			{
				cerr<<"size shoudl equal but not"<<endl;
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			if(vector[i] == p)
				return i;
		}
		return -1;
	}

	int find(std::vector<std::pair<std::pair<int,int>, int>> & suvivedClusterRelationship, std::pair<int, int> & p) 
	{
		std::pair<int, int> reversePair ;
		reversePair.first = p.second;
		reversePair.second = p.first;
		for(int i = 0; i < suvivedClusterRelationship.size(); i++)
		{
			// cout<<"i: "<<i<<"   vector[i].first: "<<vector[i].first<<"  the value to find:"<<p<<endl;
			if(suvivedClusterRelationship[i].first == p or suvivedClusterRelationship[i].first == reversePair)
				return i;
		}
		return -1;
	}	
	int find(std::vector<std::pair<int, std::pair<int,int>>> & suvivedClusterRelationship, std::pair<int, int> & p) 
	{
		std::pair<int, int> reversePair ;
		reversePair.first = p.second;
		reversePair.second = p.first;
		for(int i = 0; i < suvivedClusterRelationship.size(); i++)
		{
			// cout<<"i: "<<i<<"   vector[i].first: "<<vector[i].first<<"  the value to find:"<<p<<endl;
			if(suvivedClusterRelationship[i].second == p or suvivedClusterRelationship[i].second == reversePair)
				return i;
		}
		return -1;
	}

	int find(std::vector<std::pair<int, int> > & suvivedClusterRelationship, std::pair<int, int> & p) 
	{
		std::pair<int, int> reversePair ;
		reversePair.first = p.second;
		reversePair.second = p.first;
		for(int i = 0; i < suvivedClusterRelationship.size(); i++)
		{
			// cout<<"i: "<<i<<"   vector[i].first: "<<vector[i].first<<"  the value to find:"<<p<<endl;
			if(suvivedClusterRelationship[i] == p or suvivedClusterRelationship[i] == reversePair)
				return i;
		}
		return -1;
	}

	int find(std::vector<std::pair<int, std::pair<int,int>>> & suvivedClusterRelationship, int & p) 
	{

		for(int i = 0; i < suvivedClusterRelationship.size(); i++)
		{
			// cout<<"i: "<<i<<"   vector[i].first: "<<vector[i].first<<"  the value to find:"<<p<<endl;
			if(suvivedClusterRelationship[i].first == p)
				return i;
		}
		return -1;
	}
	int find(std::vector<std::pair<int,std::vector<int> > > & vector, int p) 
	{
		for(int i = 0; i < vector.size(); i++)
		{
			// cout<<"i: "<<i<<"   vector[i].first: "<<vector[i].first<<"  the value to find:"<<p<<endl;
			if(vector[i].first == p)
				return i;
		}
		return -1;
	}

	void updateRelationship(int consistentClusterID, std::vector<int> & conflictVector, 
		std::vector<int> & uncertainVector, std::vector<int> & consistentVector,
		std::vector<std::pair<int,std::vector<int> > > & allConflictVector, 
		std::vector<std::pair<int,std::vector<int> > > & allUncertainVector, bool & exitBit,
		std::vector<int> & existInTwoSide)
	{
		std::pair<int,std::vector<int> > toAddRelationship;
		existInTwoSide.clear();

		int conflictExist = find(allConflictVector, consistentClusterID),
			uncertainExist = find(allUncertainVector, consistentClusterID);
		
		int Lctt = conflictExist;

		//make sure find result is correct
		if(conflictExist != -1)
		{
			if(allConflictVector[conflictExist].first != consistentClusterID)
			{
				cout<<"consistentClusterID: "<<consistentClusterID<<" allConflictVector[conflictExist].first"<<allConflictVector[conflictExist].first<<endl;
				printf("the two value should be equal.");
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
		}
		if(uncertainExist != -1)
		{
			if(allUncertainVector[uncertainExist].first != consistentClusterID)
			{
				cout<<"consistentClusterID: "<<consistentClusterID<<" allUncertainVector[uncertainExist].first"<<allUncertainVector[uncertainExist].first<<endl;
				printf("the two value should be equal.");
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
		}

		//handle conflict info
		if(conflictExist == -1)
		{
			cerr<<"the cluster which is the owner of the conflict set does not exist in allConflictSet"<<endl;
			if(conflictVector.size() >0)
			{
				cerr<<"add a new conflict element to allConflictSet with consistentClusterID as "<<consistentClusterID<<endl;
				toAddRelationship.first  = consistentClusterID;
				toAddRelationship.second = conflictVector;	
				allConflictVector.push_back(toAddRelationship);
				Lctt = find(allConflictVector, consistentClusterID);
				if(Lctt != (allConflictVector.size() -1))	
				{
					cerr<<"add a new ele in allConflictVector so the serial of conflictExist should be the size - 1"<<endl;
					cerr<<"but size -1 is "<<allConflictVector.size() -1<<" and serial of consistentClusterID is "<<Lctt<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);	
				}				
			}
			else
			{
				cerr<<"the conflictVector is empty but it has some info in allConflictVector"<<endl;
			}
			// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			// cout<<"debug seg fault"<<endl;

		}
		else
		{

			cerr<<"the cluster which is the owner of the conflict set exist in allConflictSet"<<endl;
			if(conflictVector.size() >0)
			{
				Lctt = conflictExist;
				cerr<<"add to old conflict info in allConflictSet"<<endl;
				add_vector_to_vector(conflictVector, allConflictVector[conflictExist].second);	
			}
			else
			{
				Lctt = conflictExist;
				cerr<<"exist old coflict info but current conflict vector is empty"<<endl;

			}
		}

		// //reverse update  conflict info
		// if(conflictVector.size() >0)
		// {
		// 	int reverseConflictExist;
		// 	std::vector<int> toSECOND;
		// 	for(int reverseUpdateConflict = 0; reverseUpdateConflict < conflictVector.size(); reverseUpdateConflict++)
		// 	{	
		// 		toSECOND.clear();	
		// 		reverseConflictExist = find(allConflictVector, conflictVector[reverseUpdateConflict]);
		// 		if(reverseConflictExist == -1)//
		// 		{
		// 			toAddRelationship.first  = conflictVector[reverseUpdateConflict];
		// 			toSECOND.push_back(consistentClusterID);	
		// 			toAddRelationship.second = toSECOND;	
		// 			allConflictVector.push_back(toAddRelationship);
		// 		}
		// 		else
		// 		{
		// 			if(find(allConflictVector[reverseConflictExist].second, consistentClusterID) != -1)//reverse info already exist
		// 			{

		// 			} 
		// 			else//reverse info doesnt exist
		// 			{
		// 				allConflictVector[reverseConflictExist].second.push_back(consistentClusterID);
		// 			}
		// 		}
		// 	} 	
		// }

		// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
		// cout<<"debug seg fault"<<endl;

		//update existing uncertain relationship
		// if((conflictVector.size() >0) and (uncertainExist != -1))
		if( (Lctt != -1) and (uncertainExist != -1) )
		{
			int FindExistingUncertain;
			for(int Li = allConflictVector[Lctt].second.size() -1; Li >= 0 ; Li--)
			{
				// cout<<"Lctt: "<<Lctt<<endl;
				// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				// cout<<"debug seg fault"<<endl;

				FindExistingUncertain = find(allUncertainVector[uncertainExist].second, allConflictVector[Lctt].second[Li]);
				if(FindExistingUncertain != -1)
				{
					cerr<<"delete an existing uncertain cluster "<<allConflictVector[Lctt].second[Li] <<" after conflict relationship update"<<endl;
					allUncertainVector[uncertainExist].second.erase(allUncertainVector[uncertainExist].second.begin() + FindExistingUncertain);
				}
			}
		}

		// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
		// cout<<"debug seg fault"<<endl;

		//handle uncertain info
		if(uncertainVector.size() > 0)
		{
			int deleteBit;
			int secondLevelFind;
			if(Lctt != -1)//the conflict info already existed 
			{
				for(int LUC = uncertainVector.size() -1; LUC >= 0; LUC--)
				{
					// cerr<<"uncertain cluster "<<uncertainVector[LUC]<<" is under check"<<endl;
					deleteBit = find(allConflictVector[Lctt].second, uncertainVector[LUC]);
					if(deleteBit != -1)
					{
						cerr<<"find one uncertain element alreay exieted in conflict set, so delete this uncertain one"<<endl;
						secondLevelFind = find(consistentVector, uncertainVector[LUC]);
						if(secondLevelFind == -1)
							uncertainVector.erase(uncertainVector.begin() + LUC);
						else
						{
							cerr<<"the uncertain cluster alreay exists in conflict info comes from cosistentVector, which is abnormal, so we stop the program."<<endl;
							cerr<<"the existing conflixt clusters are: ";
							for(int AM1101 = 0; AM1101 < allConflictVector[conflictExist].second.size(); AM1101++)
							{
								cerr<<allConflictVector[conflictExist].second[AM1101]<<" ";
							}
							cerr<<endl;
							cout<<"at last, zihao think the uncertain should be delete and continue to run the program"<<endl;
							printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
							existInTwoSide.push_back(uncertainVector[LUC]);
							// exit(0);
							exitBit = 1;
							// sleep(3);
							uncertainVector.erase(uncertainVector.begin() + LUC);
						}
					}
					else
					{
						// cerr<<"pass"<<endl;
					}
				}
			}
			else
			{
				if(find(allConflictVector, consistentClusterID) != -1)
				{
					cout<<"find conflict info in an opposite condition"<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
			}

			if(uncertainExist == -1)
			{
				cerr<<"add a new conflict element to allConflictSet"<<endl;
				toAddRelationship.first  = consistentClusterID;
				toAddRelationship.second = uncertainVector;	
				allUncertainVector.push_back(toAddRelationship);
				int Lctta = find(allUncertainVector, consistentClusterID);
				if(Lctta != (allUncertainVector.size() -1))	
				{
					cerr<<"add a new ele in allUNcertainVector so the serial of consistentClusterID should be the size - 1"<<endl;
					cerr<<"but size -1 is "<<allUncertainVector.size() -1<<" and serial of consistentClusterID is "<<Lctta<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);	
				}

			}
			else
			{
				add_vector_to_vector(uncertainVector, allUncertainVector[uncertainExist].second);	
			}

		}
		cout<<"consistentClusterID: "<<consistentClusterID<<endl;
		
		// if(consistentClusterID == 95)
		// {
		// 	cout<<"current cluster is 95"<<endl;
		// 	if(Lctt != -1)
		// 	{
		// 		for(int findConflict = 0; findConflict < allConflictVector[Lctt].second.size(); findConflict++)
		// 		{
		// 			if(allConflictVector[Lctt].second[findConflict] == 97)
		// 			{
		// 				cout<<"find cluster 95 has conflict 97"<<endl;
		// 				int uncertain95 = find(allUncertainVector, 95);
		// 				if(uncertain95 != -1)
		// 				{
		// 					cout<<"find uncertain info of: "<<allUncertainVector[uncertain95].first<<endl;
		// 					int find97 =  find(allUncertainVector[uncertain95].second, 97);
		// 					if(find97 != -1)
		// 						exit(0);		
		// 				}

		// 			}
		// 		}
		// 	}
			
		// }

	}


	void add_vector_to_vector(std::vector<int> & fromVector, std::vector<int> & toVector)
	{
		int repeat;
		for(int i=0; i < fromVector.size(); i++)
		{
			repeat = find(toVector, fromVector[i]);
			if(repeat != -1)
				continue;
			else
				toVector.push_back(fromVector[i]);
		}
	}

	void add_int_to_vector(int & fromINt, std::vector<int> & toVector)
	{
		if(find(toVector, fromINt) != -1)
			return;
		else
			toVector.push_back(fromINt);
	}
};


#endif