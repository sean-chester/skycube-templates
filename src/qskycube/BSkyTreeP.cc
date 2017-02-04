#include "BSkyTreeP.h"


void ExecuteBSkyTree(vector<int>& AttList, vector<Point>& PointList, SNode& SkyTree)
{
	SkyTree.nLatticeID = 0;
	vector<Point> CPointList = PointList;

	ComputeSubBSkyTree(AttList, CPointList, SkyTree);
}

void ExecuteBSkyTree_pointer(vector<int>& AttList, vector<int>& PointList, SNode_pointer& SkyTree)
{
	SkyTree.nLatticeID = 0;
	vector<int> CPointList = PointList;

	ComputeSubBSkyTree_pointer(AttList, CPointList, SkyTree);
}


void ExecuteBSkyTree(vector<int>& AttList, vector<Point>& PointList, vector<Point>& Skyline)
{
	SNode SkyTree;
	SkyTree.nLatticeID = 0;

	vector<Point> CPointList = PointList;

	ComputeSubBSkyTree(AttList, CPointList, SkyTree);
	InsertSkyline(Skyline, SkyTree);
	ClearSkyTree(SkyTree);
}


void ComputeSubBSkyTree(vector<int>& AttList, vector<Point>& PointList, SNode& SkyTree)
{
	int nLatticeID, nNumChild = 0;
	vector<Point> CPointList;
	map<int, vector<Point> > PointMap;

	SelectPivotPoint(AttList, PointList);											// Pivot selection
	MapPointToRegion(AttList, PointList, PointMap, SkyTree);		// Map points to binary vectors representing subregions.

	if (!PointMap.empty())
		SkyTree.ChildNode = new SNode[PointMap.size()];

	for (map<int, vector<Point> >::iterator it = PointMap.begin(); it != PointMap.end(); it++)
	{
		nLatticeID = (*it).first;
		CPointList = (*it).second;

		if (nNumChild > 0)
		{
			SkyTree.nNumChild = nNumChild;
			PartialDominance(AttList, nLatticeID, CPointList, SkyTree);			// Partial dominance check
		}

		if (CPointList.size() == 1)
		{
			SkyTree.ChildNode[nNumChild].nLatticeID = nLatticeID;
			SkyTree.ChildNode[nNumChild].NodePointList = CPointList;
			SkyTree.ChildNode[nNumChild++].nNumChild = 0;
		}
		else if (CPointList.size() > 1)
		{
			SkyTree.ChildNode[nNumChild].nLatticeID = nLatticeID;
			ComputeSubBSkyTree(AttList, CPointList, SkyTree.ChildNode[nNumChild++]);	// Recursive call.
		}
	}

	SkyTree.nNumChild = nNumChild;
}

void ComputeSubBSkyTree_pointer(vector<int>& AttList, vector<int>& PointList, SNode_pointer& SkyTree)
{
	int nLatticeID, nNumChild = 0;
	vector<int> CPointList;
	map<int, vector<int> > PointMap;

	SelectPivotPoint_pointer(AttList, PointList);											// Pivot selection
	MapPointToRegion_pointer(AttList, PointList, PointMap, SkyTree);		// Map points to binary vectors representing subregions.

	if (!PointMap.empty())
		SkyTree.ChildNode = new SNode_pointer[PointMap.size()];

	for (map<int, vector<int> >::iterator it = PointMap.begin(); it != PointMap.end(); it++)
	{
		nLatticeID = (*it).first;
		CPointList = (*it).second;

		if (nNumChild > 0)
		{
			SkyTree.nNumChild = nNumChild;
			PartialDominance_pointer(AttList, nLatticeID, CPointList, SkyTree);			// Partial dominance check
		}

		if (CPointList.size() == 1)
		{
			SkyTree.ChildNode[nNumChild].nLatticeID = nLatticeID;
			SkyTree.ChildNode[nNumChild].NodePointList = CPointList;
			SkyTree.ChildNode[nNumChild++].nNumChild = 0;
		}
		else if (CPointList.size() > 1)
		{
			SkyTree.ChildNode[nNumChild].nLatticeID = nLatticeID;
			ComputeSubBSkyTree_pointer(AttList, CPointList, SkyTree.ChildNode[nNumChild++]);	// Recursive call.
		}
	}

	SkyTree.nNumChild = nNumChild;
}


void MapPointToRegion(vector<int>& AttList, vector<Point>& PointList, map<int, vector<Point> >& PointMap, SNode& SkyTree)
{
	int nCurAtt, nNumAtt = (int)AttList.size();
	int nNumPnt = (int)PointList.size();
	int nLattice, nEqlLattice, nPruned = (1 << nNumAtt)  - 1;

	Point CPoint, BasisPoint = PointList[0];
	SkyTree.NodePointList.push_back(PointList[0]);

	for (int nPnt = 1; nPnt < nNumPnt; nPnt++)
	{
#ifdef MEASURE
		nGMeasure++;
#endif

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		nLattice = 0, nEqlLattice = nPruned;
		CPoint = PointList[nPnt];
		for (int nAttID = 0; nAttID < nNumAtt; nAttID++)
		{
			nCurAtt = AttList[nAttID];
			if (BasisPoint[nCurAtt] < CPoint[nCurAtt])
				nLattice |= 1 << nAttID;
			else if (BasisPoint[nCurAtt] == CPoint[nCurAtt])
			{
				nLattice |= 1 << nAttID;
				nEqlLattice ^= 1 << nAttID;
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (nLattice < nPruned)
		{
			nLattice &= nEqlLattice;
			if (PointMap.find(nLattice) != PointMap.end())
				PointMap[nLattice].push_back(PointList[nPnt]);
			else
			{
				vector<Point> CPointList;
				CPointList.push_back(PointList[nPnt]);
				PointMap.insert( pair<int, vector<Point> >(nLattice, CPointList));
			}
		}
		else if (nEqlLattice == 0)
			SkyTree.NodePointList.push_back(PointList[nPnt]);
	}
}

void MapPointToRegion_pointer(vector<int>& AttList, vector<int>& PointList, map<int, vector<int> >& PointMap, SNode_pointer& SkyTree)
{
	int nCurAtt, nNumAtt = (int)AttList.size();
	int nNumPnt = (int)PointList.size();
	int nLattice, nEqlLattice, nPruned = (1 << nNumAtt)  - 1;

	int CPoint, BasisPoint = PointList[0];
	SkyTree.NodePointList.push_back(PointList[0]);

	for (int nPnt = 1; nPnt < nNumPnt; nPnt++)
	{
#ifdef MEASURE
		nGMeasure++;
#endif

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		nLattice = 0, nEqlLattice = nPruned;
		CPoint = PointList[nPnt];
		for (int nAttID = 0; nAttID < nNumAtt; nAttID++)
		{
			nCurAtt = AttList[nAttID];
			if (GlobalData[BasisPoint][nCurAtt] < GlobalData[CPoint][nCurAtt])
				nLattice |= 1 << nAttID;
			else if (GlobalData[BasisPoint][nCurAtt] == GlobalData[CPoint][nCurAtt])
			{
				nLattice |= 1 << nAttID;
				nEqlLattice ^= 1 << nAttID;
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (nLattice < nPruned)
		{
			nLattice &= nEqlLattice;
			if (PointMap.find(nLattice) != PointMap.end())
				PointMap[nLattice].push_back(PointList[nPnt]);
			else
			{
				vector<int> CPointList;
				CPointList.push_back(PointList[nPnt]);
				PointMap.insert( pair<int, vector<int> >(nLattice, CPointList));
			}
		}
		else if (nEqlLattice == 0)
			SkyTree.NodePointList.push_back(PointList[nPnt]);
	}
}

void PartialDominance(vector<int>& AttList, int nBase, vector<Point>& PointList, SNode& SkyTree)
{
	int nCLattice;
	int nNumPnt, nNumChild = SkyTree.nNumChild;

	for (int nChild = 0; nChild < nNumChild; nChild++)
	{
		nCLattice = SkyTree.ChildNode[nChild].nLatticeID;
		if (nCLattice <= nBase)
		{
			if ((nCLattice & nBase) == nCLattice)
			{
				nNumPnt = (int)PointList.size();
				for (int nPnt = 0; nPnt < nNumPnt; nPnt++)
				{
					if (!FilterPoint(PointList[nPnt], AttList, SkyTree.ChildNode[nChild]))
					{
						PointList[nPnt] = PointList[nNumPnt-1];
						PointList.pop_back();

						nPnt--, nNumPnt--;
					}
				}

				if (PointList.size() == 0)
					break;
			}
		}
		else
			break;
	}
}

void PartialDominance_pointer(vector<int>& AttList, int nBase, vector<int>& PointList, SNode_pointer& SkyTree)
{
	int nCLattice;
	int nNumPnt, nNumChild = SkyTree.nNumChild;

	for (int nChild = 0; nChild < nNumChild; nChild++)
	{
		nCLattice = SkyTree.ChildNode[nChild].nLatticeID;
		if (nCLattice <= nBase)
		{
			if ((nCLattice & nBase) == nCLattice)
			{
				nNumPnt = (int)PointList.size();
				for (int nPnt = 0; nPnt < nNumPnt; nPnt++)
				{
					if (!FilterPoint_pointer(PointList[nPnt], AttList, SkyTree.ChildNode[nChild]))
					{
						PointList[nPnt] = PointList[nNumPnt-1];
						PointList.pop_back();

						nPnt--, nNumPnt--;
					}
				}

				if (PointList.size() == 0)
					break;
			}
		}
		else
			break;
	}
}


bool FilterPoint(Point& CPoint, vector<int>& AttList, SNode& SkyNode)
{
	int nCurAtt, nNumAtt = AttList.size();
	Point SPoint = SkyNode.NodePointList[0];

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int nPruned = (1 << nNumAtt) - 1;
	int nLattice = 0, nEqlLattice = nPruned;
	for (int nAttID = 0; nAttID < nNumAtt; nAttID++)
	{
		nCurAtt = AttList[nAttID];
		if (SPoint[nCurAtt] < CPoint[nCurAtt])
			nLattice |= 1 << nAttID;
		else if (SPoint[nCurAtt] == CPoint[nCurAtt])
		{
			nLattice |= 1 << nAttID;
			nEqlLattice ^= 1 << nAttID;
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef MEASURE
	nGMeasure++;
#endif

	if (nLattice == nPruned)
		return false;
	else
	{
		nLattice &= nEqlLattice;
		if (SkyNode.nNumChild > 0)
		{
			int nCLattice, nNumChild = SkyNode.nNumChild;
			for (int nChild = 0; nChild < nNumChild; nChild++)
			{
				nCLattice = SkyNode.ChildNode[nChild].nLatticeID;
				if (nCLattice <= nLattice)
				{
					if ((nCLattice & nLattice) == nCLattice)
					{
						if (!FilterPoint(CPoint, AttList, SkyNode.ChildNode[nChild]))
							return false;
					}
				}
				else
					break;
			}
		}

		return true;
	}
}



bool FilterPoint_pointer(int& CPoint, vector<int>& AttList, SNode_pointer& SkyNode)
{
	int nCurAtt, nNumAtt = AttList.size();
	int SPoint = SkyNode.NodePointList[0];

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int nPruned = (1 << nNumAtt) - 1;
	int nLattice = 0, nEqlLattice = nPruned;
	for (int nAttID = 0; nAttID < nNumAtt; nAttID++)
	{
		nCurAtt = AttList[nAttID];
		if (GlobalData[SPoint][nCurAtt] < GlobalData[CPoint][nCurAtt])
			nLattice |= 1 << nAttID;
		else if (GlobalData[SPoint][nCurAtt] == GlobalData[CPoint][nCurAtt])
		{
			nLattice |= 1 << nAttID;
			nEqlLattice ^= 1 << nAttID;
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef MEASURE
	nGMeasure++;
#endif

	if (nLattice == nPruned)
		return false;
	else
	{
		nLattice &= nEqlLattice;
		if (SkyNode.nNumChild > 0)
		{
			int nCLattice, nNumChild = SkyNode.nNumChild;
			for (int nChild = 0; nChild < nNumChild; nChild++)
			{
				nCLattice = SkyNode.ChildNode[nChild].nLatticeID;
				if (nCLattice <= nLattice)
				{
					if ((nCLattice & nLattice) == nCLattice)
					{
						if (!FilterPoint_pointer(CPoint, AttList, SkyNode.ChildNode[nChild]))
							return false;
					}
				}
				else
					break;
			}
		}

		return true;
	}
}
