#ifndef  _SKYTREEUTIL_H
#define _SKYTREEUTIL_H

#include "Common.h"
#include <unordered_map>
#include "../skycube/NodeAuxilliary.h"

/** The node used to build SkyTrees. */
struct SNode {
	/** The bitmask of this node relative to its parent. */
	int nLatticeID;

	/** A vector of the actual children data points. */
	vector<Point> NodePointList;

	/** The number of children that this node has. */
	int nNumChild;

	/** An array of meta data relating to each of the children nodes
	 * of this node, in one-to-one correspondence with the decoupled
	 * actual data points in NodePointList.
	 */
	SNode* ChildNode;

	/** Constructs an empty SNode. */
	SNode() { }

	/** Constructs an SNode with a given bitmask, NLatticeID. */
	SNode(int nLatticeID) {
		this->nLatticeID = nLatticeID;
	}
};

struct SNode_pointer {
  int nLatticeID;
  vector<int> NodePointList;

  int nNumChild;
  SNode_pointer* ChildNode;

  SNode_pointer() {
  }
  SNode_pointer(int nLatticeID) {
    this->nLatticeID = nLatticeID;
  }
};

void InsertSkyline(vector<Point>& SkylineList, SNode& SkyNode);
void InsertSkyline_pointer(vector<int>& SkylineList, SNode_pointer& SkyNode);
void ClearSkyTree(SNode& SkyTree);
void ClearSkyList(vector<Point>& SkyTree);
void ClearSkyTree_pointer(SNode_pointer& SkyTree);
void PushStack(stack<SNode>& Stack, SNode& SkyNode);
void PushStack_pointer(stack<SNode_pointer>& Stack, SNode_pointer& SkyNode);

template<int NUM_DIMS> void InsertSkylineHashcube(std::unordered_map<int, dom_vector<NUM_DIMS> >* domvectors, SNode& SkyNode, const int cuboid)
{

	/*
	int nNumChild = (int)SkyNode.nNumChild;
	if (nNumChild > 0)
	{
		int nNumPnt = (int)SkyNode.NodePointList.size();
		for (int nPnt = 0; nPnt < nNumPnt; nPnt++){
			//SkylineList.push_back(SkyNode.NodePointList[nPnt]);
			domvectors->operator[]((int)(SkyNode.NodePointList[nPnt][0])).clear(cuboid);
		}


		//for (int nChild = 0; nChild < nNumChild; nChild++)
			//InsertSkylineHashcube(domvectors, SkyNode.ChildNode[nChild],cuboid);
	}
	else
	{
		int nNumPnt = (int)SkyNode.NodePointList.size();
		for (int nPnt = 0; nPnt < nNumPnt; nPnt++) {
			//SkylineList.push_back(SkyNode.NodePointList[nPnt]);
			domvectors->operator[]((int)(SkyNode.NodePointList[nPnt][0])).clear(cuboid);
		}
	}*/
}

template<int NUM_DIMS> void InsertSkylineHashcube(std::unordered_map<int, dom_vector<NUM_DIMS> >* domvectors, vector<Point>& Skyline, const int cuboid)
{
	/*
	for(auto it = Skyline.begin(); it != Skyline.end(); ++it){
			domvectors->operator[]((int)((*it)[0])).clear(cuboid);
	}*/
}

template<int NUM_DIMS> void InsertSkylineHashcube_pointer(std::unordered_map<int, dom_vector<NUM_DIMS> >* domvectors, SNode_pointer& SkyNode, const int cuboid)
{
	/*
	int nNumChild = (int)SkyNode.nNumChild;
	if (nNumChild > 0)
	{
		int nNumPnt = (int)SkyNode.NodePointList.size();
		for (int nPnt = 0; nPnt < nNumPnt; nPnt++){
			//SkylineList.push_back(SkyNode.NodePointList[nPnt]);
			domvectors->operator[](SkyNode.NodePointList[nPnt]).clear(cuboid);
		}


		for (int nChild = 0; nChild < nNumChild; nChild++)
			InsertSkylineHashcube_pointer(domvectors, SkyNode.ChildNode[nChild],cuboid);
	}
	else
	{
		int nNumPnt = (int)SkyNode.NodePointList.size();
		for (int nPnt = 0; nPnt < nNumPnt; nPnt++) {
			//SkylineList.push_back(SkyNode.NodePointList[nPnt]);
			domvectors->operator[](SkyNode.NodePointList[nPnt]).clear(cuboid);
		}
	}*/
}

#endif
