#include "Policy.h"


Policy::Policy() {}


Policy::Policy( const Policy &pol )  {
	type = pol.type;
	if (type == POLICY_BOOLEAN)  {
		ptr.boolptr = pol.ptr.boolptr;
	} else {
		ptr.valptr  = pol.ptr.valptr;
	}
	action   = pol.action;
	duration = pol.duration;
	factor   = pol.factor;
	init_val = pol.init_val;
	started  = pol.started;
}



Policy::Policy( int v_extent, double *v_valptr, double *v_action, double v_factor_at_end )  {
	type       = POLICY_INSTANT;
	ptr.valptr = v_valptr;
	action     = v_action;
	duration   = 0;
	factor     = v_factor_at_end;
	extent     = v_extent;
}



Policy::Policy( int v_extent, double *v_valptr, double *v_action, double v_duration, double v_factor_at_end )  {
	duration   = v_duration;
	type       = (duration == 0.0) ? POLICY_INSTANT : POLICY_LINEAR;
	ptr.valptr = v_valptr;
	action     = v_action;
	factor     = v_factor_at_end;
	extent     = v_extent;
}



Policy::Policy( int v_extent, bool *v_boolptr, double *v_action )  {
	type    = POLICY_BOOLEAN;
	ptr.boolptr = v_boolptr;
	action  = v_action;
	duration = 0;
	extent = v_extent;
}



Policy::~Policy() {}



bool  Policy::applyPolicy(double actionTime, double ctime)  {
	double  atime = *action + actionTime;
	if (type == POLICY_INSTANT)  {
		if (ctime >= atime)  {
			*ptr.valptr = factor;
			return true;
		}
		return false;
	} else if (type == POLICY_LINEAR)  {
		if (ctime >= atime && ctime <= atime+duration)  {
			if (!started)  {
				init_val = *ptr.valptr;
				started = true;
			}
//			*ptr.valptr = init_val*((ctime-atime)/duration)*factor;
			*ptr.valptr = init_val + ((ctime-atime)/duration)*(factor-init_val);
			return false;
		} else if (ctime > atime + duration)  {
			*ptr.valptr = factor;
			return true;
		} else {
			return false;
		}
	} else if (type == POLICY_BOOLEAN)  {
		if (ctime >= atime)  {
			*ptr.boolptr = !(*ptr.boolptr);
			return true;
		}
	}
	return false;
}



double Policy::getEndTime() const {
	return *action + duration;
}




PolicyQueue::PolicyQueue()  {
	actionTime = 1.e300;
}


PolicyQueue::~PolicyQueue()  {}


void  PolicyQueue::addPolicy( const Policy &v_pol )  {
	queue.push_back( v_pol );
	std::sort( queue.begin(), queue.end(), [](const Policy &a, const Policy &b) {return a.action < b.action;} );
}


void  PolicyQueue::applyPolicies( int v_extent, double time )  {
	for (auto el = queue.begin(); el != queue.end();)  {
		if (v_extent != -1 && (*el).extent != v_extent)  {
			el++;
			continue;  // Only policies of the given extent are examined, unless -1 was specified
		}

		if (*(*el).action + actionTime > time)  {
			break;
		} else {
			if (el->applyPolicy( actionTime, time ))  {
				queue.erase(el);
			} else {
				el++;
			}
		}
	}
}



void  PolicyQueue::setStart( double v_time ) {
	actionTime = v_time;
}



unsigned  PolicyQueue::size()  {
	return queue.size();
}


void  PolicyQueue::clear()  {
	queue.clear();
}


