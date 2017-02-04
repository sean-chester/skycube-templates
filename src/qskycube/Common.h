#ifndef _COMMON_H
#define _COMMON_H

#define TDS_THRESHOLD		10
#define RAND_PARAMETER		1
#define MEASURE

// Synthetic dataset parameters
#define NUM_OF_DATA		5
#define NUM_OF_DIM		4
#define NUM_OF_DIST		2
#define DEFAULT_CAR		200000
#define DEFAULT_DIM		8
#define INITIAL_CAR			100000
#define INITIAL_DIM			4
#define TOTAL_REPEAT	5

#define NUM_OF_REAL		2

#include <stdint.h>

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <algorithm>
//#include <cmath>
#include <vector>
#include <list>
#include <bitset>
#include <stack>
#include <queue>
#include <map>
#include <set>

using namespace std;

typedef float* Point;

extern int nGlobalAtt;
extern unsigned long long nGMeasure;
extern vector<int> GlobalAttList;
extern vector<Point> GlobalData;

bool ComparePointID(Point FirPoint, Point SecPoint);
bool CompareAtt(Point FirPoint, Point SecPoint);
bool CompareAtt_pointer(int FirPoint, int SecPoint);
bool CompareMultipleAtt(Point FirPoint, Point SecPoint);

struct compareAttStruct {
	compareAttStruct(int Att){ this->Att = Att; }
	bool operator () (Point FirPoint, Point SecPoint)
	{
		if (FirPoint[Att] - SecPoint[Att] < 0)
			return true;
		else if (FirPoint[Att] == SecPoint[Att] && FirPoint[0] < SecPoint[0])
			return true;
		else
			return false;
	}


	int Att;
};

struct compareAttStruct_pointer {
	compareAttStruct_pointer(int Att){ this->Att = Att; }
	bool operator () (int FirPoint, int SecPoint)
	{
		if (GlobalData[FirPoint][Att] - GlobalData[SecPoint][Att] < 0)
			return true;
		else if (GlobalData[FirPoint][Att] == GlobalData[SecPoint][Att] && FirPoint < SecPoint)
			return true;
		else
			return false;
	}


	int Att;
};

template<int BITSIZE> void SetSubspaceList(int nNumAtt, vector<int>* SubspaceList)
{
	bitset<BITSIZE> Subspace;
	int nNumCuboid = (1 << nNumAtt) - 1;

	for (int nCuboid = 0; nCuboid < nNumCuboid; nCuboid++)
	{
		Subspace = nCuboid;
//		cout << "cuboid: " << Subspace << endl;
		for (int nAttID=0; nAttID<nNumAtt; nAttID++)
		{
			if (!Subspace.test(nAttID))
				SubspaceList[nCuboid].push_back(nAttID+1);
		}
	}
}

void SortPointList(int nNumAtt, vector<Point>& PointList, vector<vector<Point> >& SPointList);
void SortPointList_pointer(int nNumAtt, vector<int>& PointList, vector<vector<int> >& SPointList);

#endif
