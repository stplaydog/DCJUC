#pragma once
#ifndef _H_HTITER_
#define _H_HTITER_

#include "ht.h"

template <class T> 
class HashTableIterator {
private:
	HashTable<T> *table;
	size_t iter_pos;
public:
	HashTableIterator(HashTable<T> *);
	~HashTableIterator();
	bool has_next();
	T next();
};

#endif
