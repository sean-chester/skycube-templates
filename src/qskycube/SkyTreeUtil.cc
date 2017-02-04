#include "SkyTreeUtil.h"


void InsertSkyline(vector<Point>& SkylineList, SNode& SkyNode)
{
	int nNumChild = (int)SkyNode.nNumChild;
	if (nNumChild > 0)
	{
		int nNumPnt = (int)SkyNode.NodePointList.size();
		for (int nPnt = 0; nPnt < nNumPnt; nPnt++)
			SkylineList.push_back(SkyNode.NodePointList[nPnt]);

		for (int nChild = 0; nChild < nNumChild; nChild++)
			InsertSkyline(SkylineList, SkyNode.ChildNode[nChild]);
	}
	else
	{
		int nNumPnt = (int)SkyNode.NodePointList.size();
		for (int nPnt = 0; nPnt < nNumPnt; nPnt++)
			SkylineList.push_back(SkyNode.NodePointList[nPnt]);
	}
}

void InsertSkyline_pointer(vector<int>& SkylineList, SNode_pointer& SkyNode)
{
	int nNumChild = (int)SkyNode.nNumChild;
	if (nNumChild > 0)
	{
		int nNumPnt = (int)SkyNode.NodePointList.size();
		for (int nPnt = 0; nPnt < nNumPnt; nPnt++)
			SkylineList.push_back(SkyNode.NodePointList[nPnt]);

		for (int nChild = 0; nChild < nNumChild; nChild++)
			InsertSkyline_pointer(SkylineList, SkyNode.ChildNode[nChild]);
	}
	else
	{
		int nNumPnt = (int)SkyNode.NodePointList.size();
		for (int nPnt = 0; nPnt < nNumPnt; nPnt++)
			SkylineList.push_back(SkyNode.NodePointList[nPnt]);
	}
}

void ClearSkyTree(SNode& SkyTree)
{
	stack<SNode> Stack;
	PushStack(Stack, SkyTree);

	while (!Stack.empty())
	{
		delete[] Stack.top().ChildNode;
		Stack.pop();
	}
}

void ClearSkyList(vector<Point>& SkyTree)
{
	vector<Point>().swap(SkyTree);
}

void ClearSkyTree_pointer(SNode_pointer& SkyTree)
{
	stack<SNode_pointer> Stack;
	PushStack_pointer(Stack, SkyTree);

	while (!Stack.empty())
	{
		delete[] Stack.top().ChildNode;
		Stack.pop();
	}
}


void PushStack(stack<SNode>& Stack, SNode& SkyNode)
{
	int nNumChild = (int)SkyNode.nNumChild;
	if (nNumChild > 0)
	{
		Stack.push(SkyNode);
		for (int nChild=0; nChild<nNumChild; nChild++)
			PushStack(Stack, SkyNode.ChildNode[nChild]);
	}
}

void PushStack_pointer(stack<SNode_pointer>& Stack, SNode_pointer& SkyNode)
{
	int nNumChild = (int)SkyNode.nNumChild;
	if (nNumChild > 0)
	{
		Stack.push(SkyNode);
		for (int nChild=0; nChild<nNumChild; nChild++)
			PushStack_pointer(Stack, SkyNode.ChildNode[nChild]);
	}
}

