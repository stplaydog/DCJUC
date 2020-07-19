#pragma once
#ifndef _H_DISMASTER
#define _H_DISMASTER

#include <stdio.h>
#include <stdlib.h>

using namespace std;

class DisMaster {
public:
	int v_size;
	int e_size[2];
	int **v_idx;
	int **e_idx;

	char *file;

public: 
	DisMaster(char *f);
	~DisMaster();
	void init_graph();

private:
	void scan_graph();
	void read_graph();
};

#endif
