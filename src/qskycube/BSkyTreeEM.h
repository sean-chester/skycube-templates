#pragma once
#ifndef  _BSKYTREE_EM_H
#define _BSKYTREE_EM_H

#include "Common.h"
#include "SkyTreeUtil.h"
#include "PivotSelection.h"

void ExecuteBSkyTreeEM(vector<int>& AttList, vector<Point>& PointList, vector<Point>* EqlPointList, SNode& SkyTree);
void ExecuteBSkyTreeEM(vector<int>& AttList, vector<Point>& PointList, vector<Point>* EqlPointList, vector<Point>& Skyline);
void ExecuteBSkyTreeEM_pointer(vector<int>& AttList, vector<int>& PointList, vector<int>* EqlPointList, SNode_pointer& SkyTree);
void ComputeSubBSkyTreeEM(vector<int>& AttList, vector<Point>& PointList, vector<Point>* EqlPointList, SNode& SkyTree);
void ComputeSubBSkyTreeEM_Pointer(vector<int>& AttList, vector<int>& PointList, vector<int>* EqlPointList, SNode_pointer& SkyTree);

void SelectPivotPointEM(vector<int>& AttList, vector<Point>& PointList);
void SelectPivotPointEM_pointer(vector<int>& AttList, vector<int>& PointList);
void MapPointToRegionEM(vector<int>& AttList, vector<Point>& PointList, map<int, vector<Point> >& PointMap, vector<Point>* EqlPointList, SNode& SkyTree);
void MapPointToRegionEM_pointer(vector<int>& AttList, vector<int>& PointList, map<int, vector<int> >& PointMap, vector<int>* EqlPointList, SNode_pointer& SkyTree);
void PartialDominanceEM(vector<int>& AttList, int nBase, vector<Point>& PointList, vector<Point>* EqlPointList, SNode& SkyTree);
void PartialDominanceEM_pointer(vector<int>& AttList, int nBase, vector<int>& PointList, vector<int>* EqlPointList, SNode_pointer& SkyTree);
bool FilterPointEM(Point& CPoint, vector<int>& AttList, vector<Point>* EqlPointList, SNode& SkyNode);
bool FilterPointEM_pointer(int& CPoint, vector<int>& AttList, vector<int>* EqlPointList, SNode_pointer& SkyNode);

#endif
