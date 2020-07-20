#include "ht_iter.h"

template <class T>
HashTableIterator::HashTableIterator(HashTable<T> *in){
	table = in;
	iter_pos =0;	
}


template <class T>
HashTableIterator::~HashTableIterator(){

}


template <class T>
bool 
HashTableIterator::has_next(){
	size_t next=iter_pos;
	for(size_t i=iter_pos; i<table->num_elem; i++){
		if(idx[i] != -1){
			next = i;
			break;
		}
	}
	if(next != iter_pos)
		return true;
	else
		return false;
}


template <class T>
T 
HashTableIterator::next(){
	size_t next=iter_pos;
	for(size_t i=iter_pos; i<table->num_elem; i++){
		if(idx[i] != -1){
			next = i;
			break;
		}
	}
	T result = table->value[iter_pos];
	iter_pos = next;
	return result;
}


template class HashTableIterator<BinTree<int> *>;
