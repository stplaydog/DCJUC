#pragma once
#ifndef _H_BTITER_
#define _H_BTITER_

#include "bt.h"

template <class T> 
class BinTreeIterator {
private:
	BinTree<T> *tree;
public:
	BinTreeIterator(BinTree<T> *);
	~BinTreeIterator();
	bool has_next();
	T next();
};

#endif
