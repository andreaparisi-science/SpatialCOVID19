#include "Policy.h"
#include "Intervention.h"



PolicyQueue::PolicyQueue()  {
	actionTime = 1.e300;
}


PolicyQueue::~PolicyQueue()  {}


void  PolicyQueue::addPolicy( Intervention *v_pol )  {
	v_pol->activate();
	queue.push_back( v_pol );
	v_pol->refQueue = this;
	updatePolicies();
}


void  PolicyQueue::updatePolicies()  {
	std::sort( queue.begin(), queue.end(), [](const Intervention *a, const Intervention *b) {return a->activationTime < b->activationTime;} );
}


void  PolicyQueue::applyPolicies( int v_extent, double time )  {
	for (auto el = queue.begin(); el != queue.end();)  {
		if (v_extent != -1 && (*el)->getExtent() != v_extent)  {
			el++;
			continue;  // Only policies of the given extent are examined, unless -1 was specified
		}

		if ((*el)->activationTime + actionTime > time)  {
			break;
		} else {
			if ((*el)->applyPolicy( actionTime, time ))  {
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


