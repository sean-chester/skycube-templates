#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <algorithm>
#include <unordered_set>
#include <vector>
#include "lattice_utilities.h"
#include "common2.h"
#include "../qskycube/Common.h"

#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#define DOM_P -1
#define DOM_Q 1
#define DOM_INCOMPARABLE 0
using namespace std;

//converts a multi vector to a matrix of floats
void redistribute_data(vector<vector<float> > datasetv, float** dataset) {
  for (unsigned int i = 0; i < datasetv.size(); i++) {
    float* x = dataset[i];
    vector<float> next = datasetv.at( i );
    for (unsigned int j = 0; j < datasetv.front().size(); j++) {
      x[j] = next[j];
    }
  }
}

void FreeDoubleArray(const unsigned row, float** matrix) {
  for (unsigned i = 0; i < row; i++)
    delete[] matrix[i];
  delete[] matrix;
}

float** AllocateDoubleArray(const unsigned row, const unsigned col) {
  float** matrix = new float*[row];
  for (unsigned i = 0; i < row; i++)
    matrix[i] = new float[col];

  return matrix;
}

vector<string> &my_split(const string &s, char delim, vector<string> &elems) {
  stringstream ss( s );
  string item;
  while ( std::getline( ss, item, delim ) ) {
    elems.push_back( item );
  }
  return elems;
}

void int_split(const string &s, char delim, vector<int> *elems) {
  stringstream ss( s );
  string item;
  while ( std::getline( ss, item, delim ) ) {
    elems->push_back( std::stoi(item) );
  }

}

vector<string> my_split(const string &s, char delim) {
  vector<string> elems;
  my_split( s, delim, elems );
  return elems;
}


//reads a number of type T from a string.
template<class T> bool from_string(T& t, const std::string& s,
		std::ios_base& (*f)(std::ios_base&)) {
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}


//used to split each line of the input into strings and convert them to floats
vector<float> split(string ins, bool line_numbers) {
	stringstream ss(ins);
	string s;
	float f;
	vector<float> returnvals;
	bool ignore_first = false;
	if (line_numbers) {
		ignore_first = true;
	}
	while (getline(ss, s, ',')) {
		if (ignore_first) {
			ignore_first = false;
		} else {
			if (from_string<float> (f, std::string(s), std::dec)) {
				returnvals.push_back(f);
			} else {
				cout << "PARSE FAILED" << endl;
			}
		}
	}
	return returnvals;
}

vector<float> split_int(string ins, bool line_numbers) {
	stringstream ss(ins);
	string s;
	vector<float> returnvals;
	bool ignore_first = false;
	if (line_numbers) {
		ignore_first = true;
	}
	while (getline(ss, s, ',')) {
		if (ignore_first) {
			ignore_first = false;
		} else {

			int numb;
			istringstream ( s ) >> numb;
			returnvals.push_back((float)numb);

		}
	}
	return returnvals;
}

vector<vector<float> > split_data(vector<vector<float> > input,int d){
	vector<vector<float> > ret;
	for(int i = 0; i < input.size(); i++){
		vector<float> tempVec;
		for(int j = 0; j < d; j++){
			tempVec.push_back(input.at(i).at(j));
		}
		ret.push_back(tempVec);
	}
	return ret;
}

//reads the file at filename and optionally removes line numbers and normalizes
//the input to a range of [0,1]
vector<vector<float> > read_data(const char *filename, bool line_number,
		bool normalize) {
	ifstream file(filename);
	string line;
	vector < vector<float> > lines;
	vector<float> tempVec;
	getline(file, line);
	tempVec = split(line, line_number);
	//how many entries do we have?
	int size = tempVec.size();
	float max[size];
	float min[size];
	for (unsigned int i = 0; i < tempVec.size(); i++) {
		max[i] = tempVec.at(i);
		min[i] = tempVec.at(i);
	}
	lines.push_back(tempVec);
	while (getline(file, line)) {
		tempVec = split(line, line_number);
		for (unsigned int i = 0; i < tempVec.size(); i++) {
			if (tempVec.at(i) > max[i])
				max[i] = tempVec.at(i);
			if (tempVec.at(i) < min[i])
				min[i] = tempVec.at(i);
		}
		lines.push_back(tempVec);
	}
	if (normalize) {
		for (unsigned int i = 0; i < lines.size(); i++) {
			tempVec = lines.at(i);
			for (unsigned int j = 0; j < tempVec.size(); j++) {
				tempVec.at(j) = (tempVec.at(j) - min[j]) / (max[j] - min[j]);
			}
			lines.at(i) = tempVec;
			//printf("%f",tempVec.at(0));
		}
	}
	return lines;
}

vector<Point> read_data2( string fname, bool contains_pids ) {
  ifstream ifs(fname.c_str());
  vector<Point> res;

  if ( !ifs.is_open() ) {
    cout << "..No such file: " << fname << endl;
    return res;
  }

  vector<string> tokens = my_split(fname, '-');
  string str_dims = tokens[2];
  string str_card = tokens[3];
  str_card = str_card.substr(0, str_card.size() - 4 );
  uint32_t dims = atoi(str_dims.c_str());
  uint32_t num_points = atoi(str_card.c_str());

  res.reserve(num_points);

  for (uint32_t i = 0; i < num_points; ++i) {
    Point p = new float[dims+1]; // p[0] is for ID
    p[0] = i;
    string line;
    getline(ifs, line);
    vector<string> dvalues = my_split(line, ',');
    for (uint32_t d = 0; d < dims; ++d) {
      p[d+1] = atof(dvalues[d].c_str() );
    }
    res.push_back(p);
  }
  ifs.close();
	return res;
}

int read_dims(const char *filename, bool line_number) {
	ifstream file(filename);
	vector < vector<float> > lines;
	if ( !file.is_open() ) {
	  cout << "..No such file: " << filename << endl;
	  return -1;
	}
	string line;
	vector<float> tempVec;
	getline(file, line);
	tempVec = split(line, line_number);
	//how many entries do we have?
	int size = tempVec.size();
	return size;
}

void ClearPointList(vector<Point>& PointList)
{
	int nNumPnt = (int)PointList.size();
	for (int nPntID = 0; nPntID < nNumPnt; nPntID++)
		delete[] PointList[nPntID];

	PointList.clear();
}

//converts a multi vector to a single vector by putting the floats from
//the multi vector into one big vector
vector<float> to_single_vector(vector<vector<float> > dataset) {
	vector<float> tempVec;
	vector<float> result;
	int k = dataset.front().size();
	for (unsigned int i = 0; i < dataset.size(); i++) {
		tempVec = dataset.at(i);
		for (int j = 0; j < k; j++) {
			result.push_back(tempVec.at(j));
		}
	}
	return result;
}

std::vector<int> getsubspace(unsigned int subspace, std::map<unsigned int,std::vector<int> > *hashcube){
	std::vector<int> ret;
	unsigned int index = subspace / 32;
	unsigned int mask  = (1 << (subspace % 32));

	for (auto it=hashcube[index].begin(); it!=hashcube[index].end(); ++it){
		//   std::cout << it->first << " => " << it->second << '\n';
		//all zero means that the mask bit was not set
		if(!(it->first & mask)){
			ret.insert(ret.end(),it->second.begin(),it->second.end());
		}
	}

	return ret;
}

template<int NUM_DIMS> vector<uint32_t>* ReorderAndConvertSkycube(uint32_t d, vector<Point>* skycube) {
  const int num_cuboids = (1 << d) - 1;
  vector<uint32_t>* reordered = new vector<uint32_t> [num_cuboids+1];


  vector<int>* subspace_list = new vector<int> [num_cuboids];
  SetSubspaceList<NUM_DIMS>(d, subspace_list);
  for (int cuboid = 0; cuboid < num_cuboids; ++cuboid) {
    bitset<NUM_DIMS> subspace(cuboid + 1);
    vector<int> dims;
    for (uint32_t i = 0; i < subspace.size(); ++i) {
      if (subspace.test(i))
        dims.push_back(i + 1);
    }

    int old_subspace = -1;
    for (int s = 0; s < num_cuboids; ++s) {
      vector<int>* dims2 = subspace_list + s;
      if (dims == (*dims2)) {
        old_subspace = s;
        break;
      }
    }
    assert(old_subspace != -1);
    for (uint32_t i = 0; i < skycube[old_subspace].size(); ++i) {
      reordered[cuboid+1].push_back(skycube[old_subspace][i][0]);
    }
  }

  delete[] subspace_list;
  return reordered;
}

#endif /* _UTILITIES_H_ */
