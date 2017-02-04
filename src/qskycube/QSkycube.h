#ifndef  _SKYCUBE_H
#define _SKYCUBE_H

#include "Common.h"
#include "BSkyTreeP.h"
#include "SkyTreeUtil.h"


void SBSkyTree_pointer(int nBase, int nPruned, vector<int>& PointList, SNode_pointer& SkyNode)
{
	bool bBase = false;
	int nNumChild = (int)SkyNode.nNumChild;
	for (int nChild = 0; nChild < nNumChild; nChild++)
	{
		if (SkyNode.ChildNode[nChild].nLatticeID == nBase)
			bBase = true;

		if (SkyNode.ChildNode[nChild].nLatticeID != nPruned)
			SBSkyTree_pointer(nBase, nPruned, PointList, SkyNode.ChildNode[nChild]);
	}

	if (!bBase)
	{
		for (auto it = SkyNode.NodePointList.begin(); it != SkyNode.NodePointList.end(); it++)
			PointList.push_back(*it);
	}
}


void SBSkyTree(int nBase, int nPruned, vector<Point>& PointList, SNode& SkyNode)
{
	bool bBase = false;
	int nNumChild = (int)SkyNode.nNumChild;
	for (int nChild = 0; nChild < nNumChild; nChild++)
	{
		if (SkyNode.ChildNode[nChild].nLatticeID == nBase)
			bBase = true;

		if (SkyNode.ChildNode[nChild].nLatticeID != nPruned)
			SBSkyTree(nBase, nPruned, PointList, SkyNode.ChildNode[nChild]);
	}

	if (!bBase)
	{
		for (vector<Point>::iterator it = SkyNode.NodePointList.begin(); it != SkyNode.NodePointList.end(); it++)
			PointList.push_back(*it);
	}
}

template<int BITSIZE> vector<Point> useSingleTree(int nNumAtt, int nCuboid, vector<int>* SubspaceList, SNode STreeCube)
{
	vector<Point> PointList;
	bitset<BITSIZE> CuboidBit = nCuboid;

	for (int nAttPos = 0; nAttPos < nNumAtt; nAttPos++)
	{
		if (CuboidBit.test(nAttPos))
		{
			int nPCuboid = 0;
			vector<int> Subspace = SubspaceList[nPCuboid];

			int nBase = 0, nNumPAtt = (int)Subspace.size();
			for (int nAttID = 0; nAttID < nNumPAtt; nAttID++)
			{
				if (Subspace[nAttID] == nAttPos + 1)
				{
					nBase = 1 << nAttID;
					break;
				}
			}
			int nPruned = (1 << nNumPAtt) - nBase - 1;

			SBSkyTree(nBase, nPruned, PointList, STreeCube);


			break;
		}
	}

	return PointList;
}

template<int BITSIZE> vector<Point> FindSingleParent(int nNumAtt, int nCuboid, vector<int>* SubspaceList, SNode* STreeCube)
{
	vector<Point> PointList;
	bitset<BITSIZE> CuboidBit = nCuboid;

	for (int nAttPos = 0; nAttPos < nNumAtt; nAttPos++)
	{
		if (CuboidBit.test(nAttPos))
		{
			int nPCuboid = nCuboid - (1 << nAttPos);
			vector<int> Subspace = SubspaceList[nPCuboid];

			int nBase = 0, nNumPAtt = (int)Subspace.size();
			for (int nAttID = 0; nAttID < nNumPAtt; nAttID++)
			{
				if (Subspace[nAttID] == nAttPos + 1)
				{
					nBase = 1 << nAttID;
					break;
				}
			}
			int nPruned = (1 << nNumPAtt) - nBase - 1;

			SBSkyTree(nBase, nPruned, PointList, STreeCube[nPCuboid]);
			break;
		}
	}

	return PointList;
}


template<int BITSIZE> vector<Point> FindSingleParentList(int nNumAtt, int nCuboid, vector<int>* SubspaceList, vector<Point>* STreeCube)
{
	vector<Point> PointList;
	bitset<BITSIZE> CuboidBit = nCuboid;

	for (int nAttPos = 0; nAttPos < nNumAtt; nAttPos++)
	{
		if (CuboidBit.test(nAttPos))
		{
			int nPCuboid = nCuboid - (1 << nAttPos);
			vector<int> Subspace = SubspaceList[nPCuboid];
/*
			int nBase = 0, nNumPAtt = (int)Subspace.size();
			for (int nAttID = 0; nAttID < nNumPAtt; nAttID++)
			{
				if (Subspace[nAttID] == nAttPos + 1)
				{
					nBase = 1 << nAttID;
					break;
				}
			}
			int nPruned = (1 << nNumPAtt) - nBase - 1;
*/
			//SBSkyTree(nBase, nPruned, PointList, STreeCube[nPCuboid]);
			PointList = STreeCube[nPCuboid];
			break;
		}
	}

	return PointList;
}

template<int BITSIZE> vector<int> FindSingleParent_pointer(int nNumAtt, int nCuboid, vector<int>* SubspaceList, SNode_pointer* STreeCube)
{
	vector<int> PointList;
	bitset<BITSIZE> CuboidBit = nCuboid;

	for (int nAttPos = 0; nAttPos < nNumAtt; nAttPos++)
	{
		if (CuboidBit.test(nAttPos))
		{
			int nPCuboid = nCuboid - (1 << nAttPos);
			vector<int> Subspace = SubspaceList[nPCuboid];

			int nBase = 0, nNumPAtt = (int)Subspace.size();
			for (int nAttID = 0; nAttID < nNumPAtt; nAttID++)
			{
				if (Subspace[nAttID] == nAttPos + 1)
				{
					nBase = 1 << nAttID;
					break;
				}
			}
			int nPruned = (1 << nNumPAtt) - nBase - 1;

			SBSkyTree_pointer(nBase, nPruned, PointList, STreeCube[nPCuboid]);
			break;
		}
	}

	return PointList;
}

template<int BITSIZE> vector<Point> FindMultipleParents(int nNumAtt, int nNumPnt, int nCuboid, vector<int>* SubspaceList, SNode* STreeCube)
{
	bitset<BITSIZE> CuboidBit = nCuboid;
	if (CuboidBit.count() == 1)
		return FindSingleParent<BITSIZE>(nNumAtt, nCuboid, SubspaceList, STreeCube);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int nParent = 0, nNumParent = CuboidBit.count();
	int * CountList = new int[nNumPnt];
	for (int nPnt = 0; nPnt < nNumPnt; nPnt++)
		CountList[nPnt] = 0;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	vector<Point> PointList;
	for (int nAttPos = 0; nAttPos < nNumAtt; nAttPos++)
	{
		if (CuboidBit.test(nAttPos))
		{
			int nPCuboid = nCuboid - (1 << nAttPos);
			vector<int> Subspace = SubspaceList[nPCuboid];

			int nBase = 0, nNumPAtt = (int)Subspace.size();
			for (int nAttID = 0; nAttID < nNumPAtt; nAttID++)
			{
				if (Subspace[nAttID] == nAttPos + 1)
				{
					nBase = 1 << nAttID;
					break;
				}
			}
			int nPruned = (1 << nNumPAtt) - nBase - 1;

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			vector<Point> CPointList;
			SBSkyTree(nBase, nPruned, CPointList, STreeCube[nPCuboid]);

			if (nParent == nNumParent - 1)
			{
				for (vector<Point>::iterator it = CPointList.begin(); it != CPointList.end(); it++)
					if (CountList[(int)(*it)[0]] == nParent) PointList.push_back(*it);
				break;
			}
			else
			{
				nParent++;
				for (vector<Point>::iterator it = CPointList.begin(); it != CPointList.end(); it++)
					CountList[(int)(*it)[0]]++;
			}
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}

	delete[] CountList;

	return PointList;
}

template<int NUM_DIMS> void ExecuteQSkycube(int nNumAtt, vector<Point>& PointList, vector<Point>* SkyCube)
{
	int nTotalCuboid = (1 << nNumAtt) - 1;

	vector<int>* SubspaceList = new vector<int>[nTotalCuboid];
	SNode* STreeCube = new SNode[nTotalCuboid];
	SetSubspaceList<NUM_DIMS>(nNumAtt, SubspaceList);		// Set the attribute ids for each subspace.

	ExecuteBSkyTree(SubspaceList[0], PointList, STreeCube[0]);
	InsertSkyline(SkyCube[0], STreeCube[0]);

	int nNumPnt = PointList.size();
	unsigned long long nMeasure = 0;
	for (int nCuboid = 1; nCuboid < nTotalCuboid; nCuboid++)
	{
		vector<Point> CPointList = FindMultipleParents<NUM_DIMS>(nNumAtt, nNumPnt, nCuboid, SubspaceList, STreeCube);

		nMeasure += (unsigned long long)CPointList.size();
		ExecuteBSkyTree(SubspaceList[nCuboid], CPointList, STreeCube[nCuboid]);
		InsertSkyline(SkyCube[nCuboid], STreeCube[nCuboid]);
	}

	nGMeasure = nMeasure;
	for (int nCuboid = 0; nCuboid < nTotalCuboid; nCuboid++)
		ClearSkyTree(STreeCube[nCuboid]);

	delete[] SubspaceList;
	delete[] STreeCube;
}





#endif
