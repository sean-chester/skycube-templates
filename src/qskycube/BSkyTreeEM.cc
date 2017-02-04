#include "BSkyTreeEM.h"

void ExecuteBSkyTreeEM(vector<int>& AttList, vector<Point>& PointList, vector<Point>* EqlPointList, SNode& SkyTree)
{
	SkyTree.nLatticeID = 0;
	vector<Point> CPointList = PointList;

	ComputeSubBSkyTreeEM(AttList, CPointList, EqlPointList, SkyTree);
}

void ExecuteBSkyTreeEM(vector<int>& AttList, vector<Point>& PointList, vector<Point>* EqlPointList, vector<Point>& Skyline)
{
	SNode SkyTree;
	SkyTree.nLatticeID = 0;

	vector<Point> CPointList = PointList;

	ComputeSubBSkyTreeEM(AttList, CPointList, EqlPointList, SkyTree);
	InsertSkyline(Skyline, SkyTree);
	ClearSkyTree(SkyTree);
}


void ExecuteBSkyTreeEM_pointer(vector<int>& AttList, vector<int>& PointList, vector<int>* EqlPointList, SNode_pointer& SkyTree)
{
	SkyTree.nLatticeID = 0;
	vector<int> CPointList = PointList;

	ComputeSubBSkyTreeEM_Pointer(AttList, CPointList, EqlPointList, SkyTree);
}


void ComputeSubBSkyTreeEM(vector<int>& AttList, vector<Point>& PointList, vector<Point>* EqlPointList, SNode& SkyTree)
{
	int nLatticeID, nNumChild = 0;
	vector<Point> CPointList;
	map<int, vector<Point> > PointMap;

	SelectPivotPointEM(AttList, PointList);																// Pivot selection
	MapPointToRegionEM(AttList, PointList, PointMap, EqlPointList, SkyTree);		// Map points to binary vectors representing subregions.

	if (!PointMap.empty())
		SkyTree.ChildNode = new SNode[PointMap.size()];

	for (map<int, vector<Point> >::iterator it = PointMap.begin(); it != PointMap.end(); it++)
	{
		nLatticeID = (*it).first;
		CPointList = (*it).second;

		if (nNumChild > 0)
		{
			SkyTree.nNumChild = nNumChild;
			PartialDominanceEM(AttList, nLatticeID, CPointList, EqlPointList, SkyTree);			// Partial dominance check
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
			ComputeSubBSkyTreeEM(AttList, CPointList, EqlPointList, SkyTree.ChildNode[nNumChild++]);	// Recursive call.
		}
	}

	SkyTree.nNumChild = nNumChild;
}

void ComputeSubBSkyTreeEM_Pointer(vector<int>& AttList, vector<int>& PointList, vector<int>* EqlPointList, SNode_pointer& SkyTree)
{
	int nLatticeID, nNumChild = 0;
	vector<int> CPointList;
	map<int, vector<int> > PointMap;

	SelectPivotPointEM_pointer(AttList, PointList);																// Pivot selection
	MapPointToRegionEM_pointer(AttList, PointList, PointMap, EqlPointList, SkyTree);		// Map points to binary vectors representing subregions.

	if (!PointMap.empty())
		SkyTree.ChildNode = new SNode_pointer[PointMap.size()];

	for (map<int, vector<int> >::iterator it = PointMap.begin(); it != PointMap.end(); it++)
	{
		nLatticeID = (*it).first;
		CPointList = (*it).second;

		if (nNumChild > 0)
		{
			SkyTree.nNumChild = nNumChild;
			PartialDominanceEM_pointer(AttList, nLatticeID, CPointList, EqlPointList, SkyTree);			// Partial dominance check
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
			ComputeSubBSkyTreeEM_Pointer(AttList, CPointList, EqlPointList, SkyTree.ChildNode[nNumChild++]);	// Recursive call.
		}
	}

	SkyTree.nNumChild = nNumChild;
}


void SelectPivotPointEM(vector<int>& AttList, vector<Point>& PointList)
{
	int nHead = 0, nTail = (int)PointList.size() - 1, nCPos = 1;
	int nCurAtt, nNumAtt = (int)AttList.size();

	Point CPoint, HPoint = PointList[nHead], Temp;
	double dCurDist, dMinDist = ComputeDistance(AttList, HPoint);

	while (nCPos <= nTail)
	{
		CPoint = PointList[nCPos];

#ifdef MEASURE
		nGMeasure++;
#endif

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		bool bDominate = false, bDominated = false;
		for (int nAttID = 0; nAttID < nNumAtt; nAttID++)
		{
			nCurAtt = AttList[nAttID];
			if (HPoint[nCurAtt] < CPoint[nCurAtt])
				bDominate = true;
			else if (HPoint[nCurAtt] > CPoint[nCurAtt])
				bDominated = true;

			if (bDominate && bDominated)
				break;
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (bDominate && !bDominated)
		{
			Temp = PointList[nCPos];
			PointList[nCPos] = PointList[nTail];
			PointList[nTail--] = Temp;
		}
		else if (!bDominate && bDominated)
		{
			Temp = PointList[nHead];
			PointList[nHead] = PointList[nCPos];
			PointList[nCPos] = PointList[nTail];
			PointList[nTail--] = Temp;

			HPoint = PointList[nHead];
			dMinDist = ComputeDistance(AttList, HPoint);
		}
		else
		{
			dCurDist = ComputeDistance(AttList, CPoint);

			if (dCurDist < dMinDist)
			{
				if (EvaluatePoint(AttList, nCPos, PointList))
				{
					Temp = PointList[nHead];
					PointList[nHead] = PointList[nCPos];
					PointList[nCPos] = Temp;

					HPoint = PointList[nHead];
					dMinDist = dCurDist;
					nCPos++;
				}
				else
				{
					Temp = PointList[nCPos];
					PointList[nCPos] = PointList[nTail];
					PointList[nTail--] = Temp;
				}
			}
			else
				nCPos++;
		}
	}
}

void SelectPivotPointEM_pointer(vector<int>& AttList, vector<int>& PointList)
{
	int nHead = 0, nTail = (int)PointList.size() - 1, nCPos = 1;
	int nCurAtt, nNumAtt = (int)AttList.size();

	int CPoint, HPoint = PointList[nHead], Temp;
	double dCurDist, dMinDist = ComputeDistance_pointer(AttList, HPoint);

	while (nCPos <= nTail)
	{
		CPoint = PointList[nCPos];

#ifdef MEASURE
		nGMeasure++;
#endif

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		bool bDominate = false, bDominated = false;
		for (int nAttID = 0; nAttID < nNumAtt; nAttID++)
		{
			nCurAtt = AttList[nAttID];
			if (GlobalData[HPoint][nCurAtt] < GlobalData[CPoint][nCurAtt])
				bDominate = true;
			else if (GlobalData[HPoint][nCurAtt] > GlobalData[CPoint][nCurAtt])
				bDominated = true;

			if (bDominate && bDominated)
				break;
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (bDominate && !bDominated)
		{
			Temp = PointList[nCPos];
			PointList[nCPos] = PointList[nTail];
			PointList[nTail--] = Temp;
		}
		else if (!bDominate && bDominated)
		{
			Temp = PointList[nHead];
			PointList[nHead] = PointList[nCPos];
			PointList[nCPos] = PointList[nTail];
			PointList[nTail--] = Temp;

			HPoint = PointList[nHead];
			dMinDist = ComputeDistance_pointer(AttList, HPoint);
		}
		else
		{
			dCurDist = ComputeDistance_pointer(AttList, CPoint);

			if (dCurDist < dMinDist)
			{
				if (EvaluatePoint_pointer(AttList, nCPos, PointList))
				{
					Temp = PointList[nHead];
					PointList[nHead] = PointList[nCPos];
					PointList[nCPos] = Temp;

					HPoint = PointList[nHead];
					dMinDist = dCurDist;
					nCPos++;
				}
				else
				{
					Temp = PointList[nCPos];
					PointList[nCPos] = PointList[nTail];
					PointList[nTail--] = Temp;
				}
			}
			else
				nCPos++;
		}
	}
}


void MapPointToRegionEM(vector<int>& AttList, vector<Point>& PointList, map<int, vector<Point> >& PointMap, vector<Point>* EqlPointList, SNode& SkyTree)
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
		else if (nEqlLattice < nPruned)
		{
			// Map non-skyline points with equal values.
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (nEqlLattice == 0)
				SkyTree.NodePointList.push_back(PointList[nPnt]);
			else
				EqlPointList[nEqlLattice].push_back(PointList[nPnt]);
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}
}

void MapPointToRegionEM_pointer(vector<int>& AttList, vector<int>& PointList, map<int, vector<int> >& PointMap, vector<int>* EqlPointList, SNode_pointer& SkyTree)
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
		else if (nEqlLattice < nPruned)
		{
			// Map non-skyline points with equal values.
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (nEqlLattice == 0)
				SkyTree.NodePointList.push_back(PointList[nPnt]);
			else
				EqlPointList[nEqlLattice].push_back(PointList[nPnt]);
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}
}


void PartialDominanceEM(vector<int>& AttList, int nBase, vector<Point>& PointList, vector<Point>* EqlPointList, SNode& SkyTree)
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
					if (!FilterPointEM(PointList[nPnt], AttList, EqlPointList, SkyTree.ChildNode[nChild]))
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

void PartialDominanceEM_pointer(vector<int>& AttList, int nBase, vector<int>& PointList, vector<int>* EqlPointList, SNode_pointer& SkyTree)
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
					if (!FilterPointEM_pointer(PointList[nPnt], AttList, EqlPointList, SkyTree.ChildNode[nChild]))
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



bool FilterPointEM(Point& CPoint, vector<int>& AttList, vector<Point>* EqlPointList, SNode& SkyNode)
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
	{
		// Map non-skyline points with equal values.
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (nEqlLattice < nPruned)
			EqlPointList[nEqlLattice].push_back(CPoint);
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		return false;
	}
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
						if (!FilterPointEM(CPoint, AttList, EqlPointList, SkyNode.ChildNode[nChild]))
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

bool FilterPointEM_pointer(int& CPoint, vector<int>& AttList, vector<int>* EqlPointList, SNode_pointer& SkyNode)
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
	{
		// Map non-skyline points with equal values.
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (nEqlLattice < nPruned)
			EqlPointList[nEqlLattice].push_back(CPoint);
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		return false;
	}
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
						if (!FilterPointEM_pointer(CPoint, AttList, EqlPointList, SkyNode.ChildNode[nChild]))
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
