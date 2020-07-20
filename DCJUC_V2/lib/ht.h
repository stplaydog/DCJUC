#include <cstdio>
#include <cstdlib>
#include <vector>


using namespace std;

template <class T, class T1>
class HashTable {
private:
	long long int *idx;
	long long int *pos;
	long long int pos_idx;
	T *key;	
	T1 *value;
	//using constant as value
	static const size_t num_elem = 10007; // a prime number
public:
	HashTable();
	~HashTable();
	void insert(T const & , T1 const & );
	void remove(T const & );
	T1 look_up(T const & );
	vector<T1> enumerate();
};
