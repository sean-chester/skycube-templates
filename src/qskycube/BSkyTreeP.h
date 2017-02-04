#ifndef  _BSKYTREEP_H
#define _BSKYTREEP_H

#include "Common.h"
#include "SkyTreeUtil.h"
#include "PivotSelection.h"

void ExecuteBSkyTree(vector<int>& AttList, vector<Point>& PointList, SNode& SkyTree);
void ExecuteBSkyTree_pointer(vector<int>& AttList, vector<int>& PointList, SNode_pointer& SkyTree);
void ExecuteBSkyTree(vector<int>& AttList, vector<Point>& PointList, vector<Point>& Skyline);

void ComputeSubBSkyTree(vector<int>& AttList, vector<Point>& PointList, SNode& SkyTree);
void ComputeSubBSkyTree_pointer(vector<int>& AttList, vector<int>& PointList, SNode_pointer& SkyTree);

void MapPointToRegion(vector<int>& AttList, vector<Point>& PointList, map<int, vector<Point> >& PointMap, SNode& SkyTree);
void MapPointToRegion_pointer(vector<int>& AttList, vector<int>& PointList, map<int, vector<int> >& PointMap, SNode_pointer& SkyTree);
void PartialDominance(vector<int>& AttList, int nBase, vector<Point>& PointList, SNode& SkyTree);
void PartialDominance_pointer(vector<int>& AttList, int nBase, vector<int>& PointList, SNode_pointer& SkyTree);
bool FilterPoint(Point& CPoint, vector<int>& AttList, SNode& SkyNode);
bool FilterPoint_pointer(int& CPoint, vector<int>& AttList, SNode_pointer& SkyNode);

#endif
