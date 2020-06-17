#ifndef  H_SPATIALCOVID_INTERVENTION
#define  H_SPATIALCOVID_INTERVENTION

#include "Policy.h"
#include <vector>


enum {POLICY_INSTANT = 1, POLICY_LINEAR, POLICY_BOOLEAN};

enum  {	POLICY_TRACING_PROB = 0,		// Contact tracing probability
		POLICY_SOCIALDIST_PROB, 		// Generalized reduction of social interactions
		POLICY_TRAVELREDUCTION, 		// Reduction of local travel intensity
		POLICY_TRAVELRED_ADMIN, 		// Reduction of inter-admin travel intensity
		POLICY_STAYATHOME_AGE, 			// Compliance of stay at home for eldest (non-working)
		POLICY_STAYATHOME_OTH, 			// Compliance of stay at home for working individuals
		POLICY_STAYATHOME_SCH, 			// Compliance of stay at home for school-aged individuals
		POLICY_FAMILY_TRANSMIT, 		// Increase in family transmission
		POLICY_STAYATHOME_FULL, 		// Whether stay-at-home for younger and older means avoiding all social contacts (ex. no shopping at all)
		POLICY_SCHOOL_CLOSURE, 			// Fraction of schools closed (generalized)
		POLICY_REDUCE_INFLIGHT,			// Stops external imports
		POLICY_TYPE_LAST
};


class  Intervention  {
public:
	Intervention();
	Intervention( const Intervention & );
	Intervention( int kind, int v_extent, double v_time, double v_duration, double v_targetValue );
	~Intervention();
	void  activate();
	void  setActivationTime( double );
	void  setDuration( double );
	double  getActivationTime();
	double  getDuration();
	bool  applyPolicy(double actionTime, double ctime);
	double getEndTime() const;
	int   getExtent();
	double activationTime;
	PolicyQueue * refQueue = nullptr;
private:
	char   type;  // 1 - InstantPolicy; 2 - LinearPolicy;
	bool   started = false;
	double init_val;
	union  {
		double *valptr;
		bool   *boolptr;
	} ptr;
	double duration;
	double factor;
	int    extent; // Global, local, sublocal;
	int kind;
};

#endif


