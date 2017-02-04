#include <unordered_map>
#include <stdint.h>
#include <vector>
#include "QSkycube.h"
#include "BSkyTreeP.h"
#include "BSkyTreeEM.h"

#pragma once
#ifndef  _SKYCUBEG_H
#define _SKYCUBEG_H


bool inverse_snoob(int a, int b) {
	return __builtin_popcount(a) < __builtin_popcount(b);
}

void AddMissingSkyPoint(int nSAttID, std::vector<Point>& PointList, std::vector<Point>& Skyline)
{
	nGlobalAtt = nSAttID;
	sort(Skyline.begin(), Skyline.end(), CompareAtt);

	std::vector<Point> CPointList;
	int nPnt = 0, nSky = 0, nNumPnt = (int)PointList.size(), nNumSky = (int)Skyline.size();
	while (nPnt < nNumPnt && nSky < nNumSky)
	{
		if (PointList[nPnt][nSAttID] < Skyline[nSky][nSAttID])
			nPnt++;
		else if (PointList[nPnt][nSAttID] > Skyline[nSky][nSAttID])
			nSky++;
		else
		{
			// Add missing points that have the same values in nSAttID attribute.
			while (PointList[nPnt][nSAttID] == Skyline[nSky][nSAttID])
			{
				CPointList.push_back(PointList[nPnt++]);
				if (nPnt == nNumPnt) break;
			}
			nSky++;
		}
	}

	Skyline = CPointList;
}

void AddMissingSkyPointP(int nSAttID, std::vector<Point>& PointList, std::vector<Point>& Skyline)
{
	//nGlobalAtt = nSAttID;
	sort(Skyline.begin(), Skyline.end(), compareAttStruct(nSAttID));

	std::vector<Point> CPointList;
	int nPnt = 0, nSky = 0, nNumPnt = (int)PointList.size(), nNumSky = (int)Skyline.size();
	while (nPnt < nNumPnt && nSky < nNumSky)
	{
		if (PointList[nPnt][nSAttID] < Skyline[nSky][nSAttID])
			nPnt++;
		else if (PointList[nPnt][nSAttID] > Skyline[nSky][nSAttID])
			nSky++;
		else
		{
			// Add missing points that have the same values in nSAttID attribute.
			while (PointList[nPnt][nSAttID] == Skyline[nSky][nSAttID])
			{
				CPointList.push_back(PointList[nPnt++]);
				if (nPnt == nNumPnt) break;
			}
			nSky++;
		}
	}

	Skyline = CPointList;
}

void AddMissingSkyPointP_pointer(int nSAttID, std::vector<int>& PointList, std::vector<int>& Skyline)
{
	//nGlobalAtt = nSAttID;

	sort(Skyline.begin(), Skyline.end(), compareAttStruct_pointer(nSAttID));

	std::vector<int> CPointList;
	int nPnt = 0, nSky = 0, nNumPnt = (int)PointList.size(), nNumSky = (int)Skyline.size();
	while (nPnt < nNumPnt && nSky < nNumSky)
	{
		if (GlobalData[PointList[nPnt]][nSAttID] < GlobalData[Skyline[nSky]][nSAttID])
			nPnt++;
		else if (GlobalData[PointList[nPnt]][nSAttID] > GlobalData[Skyline[nSky]][nSAttID])
			nSky++;
		else
		{
			// Add missing points that have the same values in nSAttID attribute.
			while (GlobalData[PointList[nPnt]][nSAttID] == GlobalData[Skyline[nSky]][nSAttID])
			{
				CPointList.push_back(PointList[nPnt++]);
				if (nPnt == nNumPnt) break;
			}
			nSky++;
		}
	}

	Skyline = CPointList;
}

void AddMissingNonSkyPoint(int nCuboid, std::vector<Point>* EqlPointList, std::vector<Point>& Skyline)
{
	// Append missing points from non-skyline points.
	for (int nPCuboid = 1; nPCuboid <= nCuboid; nPCuboid++)
	{
		if ((nPCuboid & nCuboid) == nPCuboid)
		{
			if (!EqlPointList[nPCuboid].empty())
			{
				std::vector<Point>::iterator it = Skyline.end();
				Skyline.insert(it, EqlPointList[nPCuboid].begin(), EqlPointList[nPCuboid].end());
			}
		}
	}
}


void AddMissingNonSkyPoint_pointer(int nCuboid, std::vector<int>* EqlPointList, std::vector<int>& Skyline)
{
	// Append missing points from non-skyline points.
	for (int nPCuboid = 1; nPCuboid <= nCuboid; nPCuboid++)
	{
		if ((nPCuboid & nCuboid) == nPCuboid)
		{
			if (!EqlPointList[nPCuboid].empty())
			{
				std::vector<int>::iterator it = Skyline.end();
				Skyline.insert(it, EqlPointList[nPCuboid].begin(), EqlPointList[nPCuboid].end());
			}
		}
	}
}


// Extend the QSkycube algorithm using sorted lists (also used in TDS).
template<int NUM_DIMS> void ExecuteQSkycubeGS(int nNumAtt, std::vector<Point>& PointList, std::vector<Point>* SkyCube)
{
	int nTotalCuboid = (1 << nNumAtt) - 1;

	std::vector<int>* SubspaceList = new std::vector<int>[nTotalCuboid];
	SNode* STreeCube = new SNode[nTotalCuboid];
	SetSubspaceList<NUM_DIMS>(nNumAtt, SubspaceList);		// Set the attribute ids for each subspace.

	ExecuteBSkyTree(SubspaceList[0], PointList, STreeCube[0]);
	InsertSkyline(SkyCube[0], STreeCube[0]);

	std::vector<vector<Point> > SPointList(nNumAtt + 1);
	SortPointList(nNumAtt, PointList, SPointList);		// Sort PointList on every attribute.

	int nSAttID;
	for (int nCuboid = 1; nCuboid < nTotalCuboid; nCuboid++)
	{
		std::vector<Point> CPointList = FindSingleParent<NUM_DIMS>(nNumAtt, nCuboid, SubspaceList, STreeCube);

		// Add missing points from the sorted list (entire dataset).
		nSAttID = SubspaceList[nCuboid][0];
		AddMissingSkyPoint(nSAttID, SPointList[nSAttID], CPointList);

		ExecuteBSkyTree(SubspaceList[nCuboid], CPointList, STreeCube[nCuboid]);
		InsertSkyline(SkyCube[nCuboid], STreeCube[nCuboid]);
	}

	for (int nCuboid = 0; nCuboid < nTotalCuboid; nCuboid++)
		ClearSkyTree(STreeCube[nCuboid]);

	delete[] SubspaceList;
	delete[] STreeCube;
}



// Extend the QSkycube algorithm using the equivalence lattice).
template<int NUM_DIMS> void ExecuteQSkycubeGL(int nNumAtt, std::vector<Point>& PointList, std::vector<Point>* SkyCube)
{
	int nTotalCuboid = (1 << nNumAtt) - 1;

	std::vector<int>* SubspaceList = new std::vector<int>[nTotalCuboid];
	SNode* STreeCube = new SNode[nTotalCuboid];
	SetSubspaceList<NUM_DIMS>(nNumAtt, SubspaceList);		// Set the attribute ids for each subspace.


	// Exploit a lattice to cover points with equal values.
	std::vector<Point>* EqlPointList = new std::vector<Point>[nTotalCuboid];
	ExecuteBSkyTreeEM(SubspaceList[0], PointList, EqlPointList, STreeCube[0]);
	InsertSkyline(SkyCube[0], STreeCube[0]);
	std::vector<int> iteration_order;


	std::vector<vector<Point> > SPointList(nNumAtt + 1);
	SortPointList(nNumAtt, SkyCube[0], SPointList);	// Sort PointList on every single attribute.

	int nSAttID;
	unsigned long long nNumTotalPnt = 0;
	for (int nCuboid = 1; nCuboid < nTotalCuboid; nCuboid++){
		std::vector<Point> CPointList = FindSingleParent<NUM_DIMS>(nNumAtt, nCuboid, SubspaceList, STreeCube);

		// Add missing points from the sorted list (full space skyline).
		nSAttID = SubspaceList[nCuboid][0];
		AddMissingSkyPoint(nSAttID, SPointList[nSAttID], CPointList);

		// Add missing points from the lattice.
		AddMissingNonSkyPoint(nCuboid, EqlPointList, CPointList);

		nNumTotalPnt += (unsigned long long)CPointList.size();
		ExecuteBSkyTree(SubspaceList[nCuboid], CPointList, STreeCube[nCuboid]);

		InsertSkyline(SkyCube[nCuboid], STreeCube[nCuboid]);
	}

	nGMeasure = 0;
	nGMeasure = nNumTotalPnt;

	for (int nCuboid = 0; nCuboid < nTotalCuboid; nCuboid++)
		ClearSkyTree(STreeCube[nCuboid]);

	delete[] EqlPointList;
	delete[] SubspaceList;
	delete[] STreeCube;
}

//with memory usage reduction.
template<int NUM_DIMS> void ExecuteQSkycubeGLPMR(int nNumAtt, std::vector<Point>& PointList, std::vector<Point>* SkyCube, int t, int max_d)
{
	int nTotalCuboid = (1 << nNumAtt) - 1;

	std::vector<int>* SubspaceList = new std::vector<int>[nTotalCuboid];
	SNode* STreeCube = new SNode[nTotalCuboid];
	SetSubspaceList<NUM_DIMS>(nNumAtt, SubspaceList);		// Set the attribute ids for each subspace.


	// Exploit a lattice to cover points with equal values.
	std::vector<Point>* EqlPointList = new std::vector<Point>[nTotalCuboid];
	ExecuteBSkyTreeEM(SubspaceList[0], PointList, EqlPointList, STreeCube[0]);
	InsertSkyline(SkyCube[0], STreeCube[0]);
	std::vector<int> iteration_order;
	for(int i = 1; i < nTotalCuboid; i++){
		iteration_order.push_back(i);
	}
	std::sort(iteration_order.begin(),iteration_order.end(),inverse_snoob);

	std::vector<vector<int> > lattice_level_order;
	int last_popc = 0;
	for (auto it = iteration_order.begin(); it != iteration_order.end(); ++it) {
		int next_popc = __builtin_popcount(*it);

		if (next_popc != last_popc) {
			std::vector<int> next_list;
			next_list.push_back(*it);
			lattice_level_order.push_back(next_list);
			last_popc = next_popc;
		} else {
			lattice_level_order[lattice_level_order.size() - 1].push_back(*it);
		}
	}

	std::vector<vector<Point> > SPointList(nNumAtt + 1);
	SortPointList(nNumAtt, SkyCube[0], SPointList);	// Sort PointList on every single attribute.

	unsigned long long nNumTotalPnt = 0;

	int start_level = 0;
	if(max_d < nNumAtt) {
		start_level = nNumAtt - 1 - max_d;
	}

	for (unsigned int i = start_level; i < lattice_level_order.size(); i++) {

		if(i > 2+start_level){

			//clean up trees no longer used to minimize memory footprint
			for(auto it = lattice_level_order[i-2].begin(); it != lattice_level_order[i-2].end(); ++it){
				ClearSkyTree(STreeCube[*it]);
			}

		}

		if((i == 2+start_level)){
			ClearSkyTree(STreeCube[0]);
		}


#pragma omp parallel for schedule(dynamic, 1) num_threads( t )
		for (unsigned int j = 0; j < lattice_level_order[i].size(); j++) {
			int nCuboid = lattice_level_order[i][j];
			std::vector<Point> CPointList;
			//the first level we compute for
			if(i > start_level) {
				CPointList = FindSingleParent<NUM_DIMS>(nNumAtt, nCuboid, SubspaceList, STreeCube);
			} else {
				//use full skyline
				for(auto it = SkyCube[0].begin(); it != SkyCube[0].end(); it++) {
					CPointList.push_back(*it);
				}
			}
			// Add missing points from the sorted list (full space skyline).
			int nSAttID = SubspaceList[nCuboid][0];
			AddMissingSkyPointP(nSAttID, SPointList[nSAttID], CPointList);

			// Add missing points from the lattice.
			AddMissingNonSkyPoint(nCuboid, EqlPointList, CPointList);
			ExecuteBSkyTree(SubspaceList[nCuboid], CPointList, STreeCube[nCuboid]);

			InsertSkyline(SkyCube[nCuboid], STreeCube[nCuboid]);
		}
	}

	nGMeasure = 0;
	nGMeasure = nNumTotalPnt;

	delete[] EqlPointList;
	delete[] SubspaceList;
}


#endif
