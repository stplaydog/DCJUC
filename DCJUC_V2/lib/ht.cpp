#include "ht.h"
#include "bt.h"

/**
 * allocate memory for the hash table
 * **/
template <class T, class T1>
HashTable<T,T1>::HashTable(){
	idx = (long long int*)malloc(sizeof(long long int)*num_elem);
	pos = (long long int*)malloc(sizeof(long long int)*num_elem);
	key = (T*)malloc(sizeof(T)*num_elem);
	value = (T1*)malloc(sizeof(T1)*num_elem);
	pos_idx = 0;
	for(size_t i=0; i<num_elem; i++)
		idx[i] = -1;
}

/**
 * destructor
 * **/
template <class T, class T1>
HashTable<T,T1>::~HashTable(){
	free(idx);
	free(pos);
	free(key);
	free(value);
}

/**
 * insert someting using hash table
 * 1) method of indexing:     pos = value % prime
 * 2) method of confliction:  pos = value % prime + i (i=0 to num_elem)
 * **/
template <class T, class T1>
void 
HashTable<T,T1>::insert(T const & k, T1 const & v){
	int position = k % num_elem; //get the position
	if(idx[position] == -1){ //no problem
		key[position] = k; 
		value[position] = v;
	}
	else{ //confliction
		for(size_t i=position; i<num_elem; i++){
			if(idx[i] == -1){
				key[position] = k; 
				value[position] = v;
				break;
			}
		}
	}
}

/**
 * remove element from hash table by value
 * hash method is the same as mentioned above, 
 * the only difference is that value needs tobe checked
 * **/
template <class T, class T1>
void 
HashTable<T,T1>::remove(T const & k){
	int position = k % num_elem; //get the position
	if(key[position] == k){ //no problem
		idx[position] = -1;
	}
	else{ //confliction
		for(size_t i=position; i<num_elem; i++){
			if(key[i] == k){
				idx[i] = -1;
				break;
			}
		}
	}
}

/**
 * check if a value is in the hash table
 * **/
template <class T, class T1>
T1 
HashTable<T,T1>::look_up(T const & k){
	int position = k % num_elem; //get the position
	if(key[position] == k){ //no problem
		return value[position];
	}
	else{ //confliction
		for(size_t i=position; i<num_elem; i++){
			if(key[i] == k){
				return value[i];
			}
		}
	}
}


template <class T, class T1>
vector<T1> 
HashTable<T,T1>::enumerate(){
	vector<T1> result;
	return result;
}

/**
 * instantiate something
 * **/
template class HashTable<int, BinTree<int> *>;
