#include "PivotSelection.h"


void SelectPivotPoint(vector<int>& AttList, vector<Point>& PointList)
{
	int nHead = 0, nTail = (int)PointList.size() - 1, nCPos = 1;
	int nCurAtt, nNumAtt = (int)AttList.size();

	Point CPoint, HPoint = PointList[nHead], Temp;
	float dCurDist, dMinDist = ComputeDistance(AttList, HPoint);

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
			PointList[nCPos] = PointList[nTail];
			PointList.pop_back();
			nTail--;
		}
		else if (!bDominate && bDominated)
		{
			PointList[nHead] = PointList[nCPos];
			PointList[nCPos] = PointList[nTail];
			PointList.pop_back();
			nTail--;

			HPoint = PointList[nHead];
			dMinDist = ComputeDistance(AttList, HPoint);
			nCPos = nHead + 1; // RESTART (BUG-FIX!)
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
					PointList[nCPos] = PointList[nTail];
					PointList.pop_back();
					nTail--;
				}
			}
			else
				nCPos++;
		}
	}
}

void SelectPivotPoint_pointer(vector<int>& AttList, vector<int>& PointList)
{
	int nHead = 0, nTail = (int)PointList.size() - 1, nCPos = 1;
	int nCurAtt, nNumAtt = (int)AttList.size();

	int CPoint, HPoint = PointList[nHead], Temp;
	float dCurDist, dMinDist = ComputeDistance_pointer(AttList, HPoint);

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
			PointList[nCPos] = PointList[nTail];
			PointList.pop_back();
			nTail--;
		}
		else if (!bDominate && bDominated)
		{
			PointList[nHead] = PointList[nCPos];
			PointList[nCPos] = PointList[nTail];
			PointList.pop_back();
			nTail--;

			HPoint = PointList[nHead];
			dMinDist = ComputeDistance_pointer(AttList, HPoint);
			nCPos = nHead + 1; // RESTART (BUG-FIX!)
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
					PointList[nCPos] = PointList[nTail];
					PointList.pop_back();
					nTail--;
				}
			}
			else
				nCPos++;
		}
	}
}


float ComputeDistance(vector<int>& AttList, Point CPoint)
{
	int nCurAtt, nNumAtt = (int)AttList.size();
	float dMax, dMin;

	dMax = dMin = CPoint[AttList[0]];
	for (int nAttID = 1; nAttID < nNumAtt; nAttID++)
	{
		nCurAtt = AttList[nAttID];
		if (dMax < CPoint[nCurAtt])
			dMax = CPoint[nCurAtt];
		else if (dMin > CPoint[nCurAtt])
			dMin = CPoint[nCurAtt];
	}

	return dMax - dMin;
}

float ComputeDistance_pointer(vector<int>& AttList, int CPoint)
{
	int nCurAtt, nNumAtt = (int)AttList.size();
	float dMax, dMin;

	dMax = dMin = GlobalData[CPoint][AttList[0]];
	for (int nAttID = 1; nAttID < nNumAtt; nAttID++)
	{
		nCurAtt = AttList[nAttID];
		if (dMax < GlobalData[CPoint][nCurAtt])
			dMax = GlobalData[CPoint][nCurAtt];
		else if (dMin > GlobalData[CPoint][nCurAtt])
			dMin = GlobalData[CPoint][nCurAtt];
	}

	return dMax - dMin;
}

bool EvaluatePoint_pointer(vector<int>& AttList, int nCPos, vector<int>& PointList)
{
	bool bSkyline = true, bDominate/*, eq*/;
	int nCurAtt, nNumAtt = (int)AttList.size();

	int CPoint = PointList[nCPos], SPoint;
	for (int nPnt = 1; nPnt < nCPos; nPnt++)
	{
#ifdef MEASURE
		nGMeasure++;
#endif

		SPoint = PointList[nPnt];

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		bDominate = true;/* eq = true;*/
		for (int nAttID = 0; nAttID < nNumAtt; nAttID++)
		{
			nCurAtt = AttList[nAttID];
			if (GlobalData[CPoint][nCurAtt] < GlobalData[SPoint][nCurAtt])
			{
				bDominate = false;
				/*eq = false;*/
				break;
			} /*else if (CPoint[nCurAtt] > SPoint[nCurAtt]) {
			    eq = false;
			}*/
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (bDominate /*&& !eq*/)
		{
			bSkyline = false;
			break;
		}
	}

	return bSkyline;
}


bool EvaluatePoint(vector<int>& AttList, int nCPos, vector<Point>& PointList)
{
	bool bSkyline = true, bDominate/*, eq*/;
	int nCurAtt, nNumAtt = (int)AttList.size();

	Point CPoint = PointList[nCPos], SPoint;
	for (int nPnt = 1; nPnt < nCPos; nPnt++)
	{
#ifdef MEASURE
		nGMeasure++;
#endif

		SPoint = PointList[nPnt];

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		bDominate = true;/* eq = true;*/
		for (int nAttID = 0; nAttID < nNumAtt; nAttID++)
		{
			nCurAtt = AttList[nAttID];
			if (CPoint[nCurAtt] < SPoint[nCurAtt])
			{
				bDominate = false;
				/*eq = false;*/
				break;
			} /*else if (CPoint[nCurAtt] > SPoint[nCurAtt]) {
			    eq = false;
			}*/
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (bDominate /*&& !eq*/)
		{
			bSkyline = false;
			break;
		}
	}

	return bSkyline;
}

