#ifndef PIVOT_SELECTION_H
#define PIVOT_SELECTION_H

#include "Common.h"
#include "SkyTreeUtil.h"

void SelectPivotPoint(vector<int>& AttList, vector<Point>& PointList);
void SelectPivotPoint_pointer(vector<int>& AttList, vector<int>& PointList);
float ComputeDistance(vector<int>& AttList, Point CPoint);
float ComputeDistance_pointer(vector<int>& AttList, int CPoint);
bool EvaluatePoint(vector<int>& AttList, int nCPos, vector<Point>& PointList);
bool EvaluatePoint_pointer(vector<int>& AttList, int nCPos, vector<int>& PointList);

#endif
