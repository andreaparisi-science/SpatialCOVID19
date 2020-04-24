#include <vector>
#include <algorithm>

enum {POLICY_INSTANT = 1, POLICY_LINEAR, POLICY_BOOLEAN};
enum  {POLICY_TYPE_RESERVED = 0, POLICY_TYPE_GLOBAL = 1, POLICY_TYPE_LOCAL, POLICY_TYPE_SIZE};

class  Policy  {
public:
	Policy();
	Policy( const Policy &pol );
	Policy( int v_extent, double *v_valptr, double *v_action, double v_factor_at_end );
	Policy( int v_extent, double *v_valptr, double *v_action, double v_duration, double v_factor_at_end );
	Policy( int v_extent, bool *v_boolptr, double *v_action );
	~Policy();

	bool  applyPolicy(double actionTime, double ctime);
	double getEndTime() const;

	char   type;  // 1 - InstantPolicy; 2 - LinearPolicy;
	bool   started = false;
	double init_val;
	union  {
		double *valptr;
		bool   *boolptr;
	} ptr;
	double *action;
	double duration;
	double factor;
	int    extent; // Global, local, sublocal;
};





class  PolicyQueue  {
public:
	PolicyQueue();
	~PolicyQueue();
	void  setStart( double v_time );
	void  addPolicy( const Policy &v_pol );
	void  applyPolicies( int v_ext, double time );
	unsigned  size();
private:
	double  actionTime;
	std::vector<Policy>  queue;
};




