#include <cstdio>
#include <cstdlib>
#include <vector>


using namespace std;

template <class T>
class BinTree {
private:
	T value;	
	size_t num_elem;
	BinTree *left;
	BinTree *right;
	BinTree *parent;
public:
	BinTree(T);
	~BinTree();
	void insert(T const & v);
	void remove(T const & v);
	BinTree * look_up(T const & v);
	vector<T> enumerate();
	BinTree * get_inorder_predecessor(BinTree *p);
	BinTree * get_inorder_successor(BinTree *p);
};
