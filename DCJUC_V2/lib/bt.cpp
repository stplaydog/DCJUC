#include "bt.h"


/**
 * this is the constructor
 * **/
template <class T>
BinTree<T>::BinTree(T v){
	value = v; 
}

/**
 * this is the destructor
 * **/
template <class T>
BinTree<T>::~BinTree(){
}

/**
 * this is to insert an element into the binary search tree
 * **/
template <class T>
void
BinTree<T>::insert(T const & v){
	BinTree *tl = left;
	BinTree *tr = right;
	BinTree *tp;
	T vp = value;
	while(true){
		if(v < vp && tl != NULL){
			tp = tl;
			vp = tp->value;
			tl = tp->left;
			tr = tp->right;
		}
		else if(v > vp && tr != NULL){
			tp = tr;
			vp = tp->value;
			tl = tp->left;
			tr = tp->right;
		}
		else if(v < vp && tl == NULL){
			BinTree *node = new BinTree(v);
			tl = node;
			break;
		}
		else if(v > vp && tr == NULL){
			BinTree *node = new BinTree(v);
			tr = node;
			break;
		}
	}
}

/**
 * this is to delete an element from the binary search tree
 * [cite from wikipedia]
 * **/
template <class T>
void 
BinTree<T>::remove(T const & v){
	//traversal then divided into three cases
	bool found = false;
	BinTree *tl = left;
	BinTree *tr = right;
	BinTree *tp;
	T vp = value;
	while(found == false){
		if(v == vp){
			//divide into three cases
			//case 1, Deleting a leaf (node with no children)
			if(tl==NULL && tr == NULL){
				delete tp;
				break;
			}
			//case 2, Deleting a node with one child
			else if(tl!=NULL && tr == NULL){
				tp = tl;
				break;
			}
			else if(tl==NULL && tr != NULL){
				tp = tr;
				break;
			}
			//case 3, Deleting a node with two children
			//Call the node to be deleted N. Do not delete N. 
			//Instead, 
			// 1) choose either its in-order successor node 
			// 2) or its in-order predecessor node, R. 
			// 3) Copy the value of R to N, 
			// 4) then recursively call delete on R 
			// 5) until reaching one of the first two cases.
			else{
				//check inorder successor first
				BinTree *succ = get_inorder_successor(tp);
				if(succ != NULL){
					value = succ->value;
					tp = succ;
					vp = v;
					tl = succ->left;
					tr = succ->right;
					continue;
				}
				//check inorder predecessor then 
				else{
					BinTree *pred = get_inorder_predecessor(tp);
					value = pred->value;
					tp = pred;
					vp = v;
					tl = pred->left;
					tr = pred->right;
					continue;
				}
			}
		}
		else if(v < vp){
			if(tl == NULL){
				printf("can't find the result!\n");
				return;
			}
			tp = tl;
			vp = tp->value;
			tl = tp->left;
			tr = tp->right;
		}
		else if(v > vp){
			if(tr == NULL){
				printf("can't find the result!\n");
				return;
			}
			tp = tr;
			vp = tp->value;
			tl = tp->left;
			tr = tp->right;
		}
	}
}	

/**
 * this is to check if an element is in the binary search tree
 * **/
template <class T>
BinTree<T> * 
BinTree<T>::look_up(T const &){
	return this;
}

/**
 * this is to return a list of the leaf in the binary search tree
 * **/
template <class T>
vector<T> BinTree<T>::enumerate(){
	vector<T> result;
	vector<BinTree*> stack;
	stack.push_back(this);
	BinTree *current = this->left;
	while(stack.size()>0){
		if(current != NULL){
			stack.push_back(current);
			current = current->left;
		}
		else{
			current = stack[stack.size()-1];
			stack.pop_back();
			result.push_back(current->value);
			current = current->right;
		}
	}
	return result;
}

/**
 * this is to get predecessor of inorder traversal
 * 1) If u has a left child, l, then pred(u) is the rightmost descendent of l 
 * 2) Otherwise, pred(u) is the closest ancestor, v, of u (if any) 
 *    such that u is descended from the right child of v.
 * 3) If there is no such ancestor, then pred(u) is undefined.
 * **/
template <class T>
BinTree<T> * 
BinTree<T>::get_inorder_predecessor(BinTree *p){
	//if p has left child
	if(p->left != NULL){
		BinTree *l = p->left;
		BinTree *r = l->right;
		while(r != NULL){
			l = l->right;
			r = r->right;
		}
		return l;
	}
	//else
	else{
		BinTree *f = p->parent;
		BinTree *ff = f->parent;
		if(f == ff->right)
			return ff;
		else
			return NULL; //undefined
	}
}

/**
 * this is to get successor of inorder traversal
 * [cite from Thomas A. Anastasio's note]:
 * 1) If u has a right child, r, then succ(u) is the leftmost descendent of r,
 * 2) Otherwise, succ(u) is the closest ancestor, 
 *    v, of u (if any) such that u is descended from the left child of v. 
 * 3) If there is no such ancestor, then succ(u) is undefined.
 * **/
template <class T>
BinTree<T> * 
BinTree<T>::get_inorder_successor(BinTree *p){
	//if p has a right child	
	if(p->right != NULL){
		BinTree *r = p->right;
		BinTree *l = r->left;
		while (l != NULL){
			r = r->left;
			l = l->left;
		}
		return r; 
	}
	//else, 
	else{
		BinTree *f = p->parent;
		BinTree *ff = f->parent;
		if(f == ff->left)
			return ff;
		else
			return NULL; //undefined
	}
}

template class BinTree<int>;
