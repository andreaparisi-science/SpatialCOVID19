#ifndef  H_SPATIALCOVID_POLICY
#define  H_SPATIALCOVID_POLICY

#include <vector>
#include <algorithm>

//enum {POLICY_INSTANT = 1, POLICY_LINEAR, POLICY_BOOLEAN};
//enum  {POLICY_TYPE_RESERVED = 0, POLICY_TYPE_GLOBAL = 1, POLICY_TYPE_LOCAL, POLICY_TYPE_SIZE};
class Intervention;

class  PolicyQueue  {
public:
	PolicyQueue();
	~PolicyQueue();
	void  setStart( double v_time );
	void  addPolicy( Intervention *v_pol );
	void  applyPolicies( int v_ext, double time );
	void  updatePolicies();
	unsigned  size();
	void  clear();
private:
	double  actionTime;
	std::vector<Intervention*>  queue;
};

#endif



