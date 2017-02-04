#include "Common.h"

int nGlobalAtt;
unsigned long long nGMeasure;
vector<int> GlobalAttList;
vector<Point> GlobalData;


bool ComparePointID(Point FirPoint, Point SecPoint)
{
	if (FirPoint[0] < SecPoint[0])
		return true;
	else
		return false;
}


bool CompareAtt(Point FirPoint, Point SecPoint)
{
	if (FirPoint[nGlobalAtt] - SecPoint[nGlobalAtt] < 0)
		return true;
	else if (FirPoint[nGlobalAtt] == SecPoint[nGlobalAtt] && FirPoint[0] < SecPoint[0])
		return true;
	else
		return false;
}

bool CompareAtt_pointer(int FirPoint, int SecPoint)
{
	if (GlobalData[FirPoint][nGlobalAtt] - GlobalData[SecPoint][nGlobalAtt] < 0)
		return true;
	else if (GlobalData[FirPoint][nGlobalAtt] == GlobalData[SecPoint][nGlobalAtt] && FirPoint < SecPoint)
		return true;
	else
		return false;
}

bool CompareMultipleAtt(Point FirPoint, Point SecPoint)
{
	for (vector<int>::iterator it = GlobalAttList.begin(); it != GlobalAttList.end(); it++)
	{
		if (FirPoint[*it] - SecPoint[*it] < 0)
			return true;
		else if (FirPoint[*it] - SecPoint[*it] > 0)
			return false;
	}

	if (FirPoint[0] < SecPoint[0])
		return true;
	else
		return false;
}

void SortPointList(int nNumAtt, vector<Point>& PointList, vector<vector<Point> >& SPointList)
{
	for (int nAttID = 1; nAttID <= nNumAtt; nAttID++)
	{
		nGlobalAtt = nAttID;
		sort(PointList.begin(), PointList.end(), CompareAtt);
		SPointList[nAttID] = PointList;
	}
}

void SortPointList_pointer(int nNumAtt, vector<int>& PointList, vector<vector<int> >& SPointList)
{
	for (int nAttID = 1; nAttID <= nNumAtt; nAttID++)
	{
		nGlobalAtt = nAttID;
		sort(PointList.begin(), PointList.end(), CompareAtt_pointer);
		SPointList[nAttID] = PointList;
	}
}

