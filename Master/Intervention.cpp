#include "Intervention.h"




// Vector storing policies (to be filled in country specific c++ setup files)
extern std::vector<double>  TRACING_PROB;
extern std::vector<double>  SOCIALDIST_PROB;
extern std::vector<double>  TRAVELREDUCTION;
extern std::vector<double>  TRAVELRED_ADMIN;
extern std::vector<double>  STAYATHOME_AGE;
extern std::vector<double>  STAYATHOME_OTH;
extern std::vector<double>  STAYATHOME_SCH;
extern std::vector<double>  FAMILY_TRANSMIT;
extern std::vector<double>  STAYATHOME_FULL;
extern std::vector<double>  SCHOOL_CLOSURE;


Intervention::Intervention() {}



Intervention::Intervention( const Intervention &intv )  {
	extent = intv.extent;
	activationTime = intv.activationTime;
	duration = intv.duration;
	type   = intv.type;
	factor = intv.factor;
	ptr.valptr = intv.ptr.valptr;
	kind   = intv.kind;
}



Intervention::~Intervention() {}


Intervention::Intervention( int v_kind, int v_extent, double v_time, double v_duration, double v_targetValue )  {
	extent = v_extent;
	activationTime = v_time;
	duration = v_duration;
	if (duration == 0.0)  {
		type = POLICY_INSTANT;
	} else {
		type = POLICY_LINEAR;
	}
	factor = v_targetValue;
	kind   = v_kind;
}



// This must be activated after vecors have been properly sized
void  Intervention::activate()  {
	switch( kind )  {
		case  POLICY_TRACING_PROB:
			ptr.valptr = &TRACING_PROB[extent];
			break;
		case  POLICY_SOCIALDIST_PROB:
			ptr.valptr = &SOCIALDIST_PROB[extent];
			break;
		case  POLICY_TRAVELREDUCTION:
			ptr.valptr = &TRAVELREDUCTION[extent];
			break;
		case  POLICY_TRAVELRED_ADMIN:
			ptr.valptr = &TRAVELRED_ADMIN[extent];
			break;
		case  POLICY_STAYATHOME_AGE:
			ptr.valptr = &STAYATHOME_AGE[extent];
			break;
		case  POLICY_STAYATHOME_OTH:
			ptr.valptr = &STAYATHOME_OTH[extent];
			break;
		case  POLICY_STAYATHOME_SCH:
			ptr.valptr = &STAYATHOME_SCH[extent];
			break;
		case  POLICY_FAMILY_TRANSMIT:
			ptr.valptr = &FAMILY_TRANSMIT[extent];
			break;
		case  POLICY_STAYATHOME_FULL:
			ptr.valptr = &STAYATHOME_FULL[extent];
			break;
		case  POLICY_SCHOOL_CLOSURE:
			ptr.valptr = &SCHOOL_CLOSURE[extent];
			break;
		default:
			throw  "Unknown Policy type";
	}
}



void  Intervention::setActivationTime( double v_time )  {
	activationTime = v_time;
}



double  Intervention::getActivationTime()  {
	return  activationTime;
}



bool  Intervention::applyPolicy(double actionTime, double ctime)  {
	double  atime = activationTime + actionTime;
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



double Intervention::getEndTime() const {
	return activationTime + duration;
}



int    Intervention::getExtent()  {
	return extent;
}




