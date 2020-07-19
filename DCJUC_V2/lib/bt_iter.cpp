#include "bt_iter.h"

template <class T>
BinTreeIterator::BinTreeIterator(BinTree<T> *in){
	tree = in;
}


template <class T>
BinTreeIterator::~BinTreeIterator(){

}


template <class T>
bool 
BinTreeIterator::has_next(){
	BinTree<T> *next = tree->get_inorder_successor(tree);
	if(next != NULL)
		return true;
	else
		return false;
}


template <class T>
T 
BinTreeIterator::next(){
	BinTree<T> *next = tree->get_inorder_successor(tree);
	T result = tree->value;
	tree = next;
	return result;
}


template class BinTreeIterator<int>;
