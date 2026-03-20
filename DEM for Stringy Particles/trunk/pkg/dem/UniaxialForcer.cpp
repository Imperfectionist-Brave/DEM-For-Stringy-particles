// 2008 © Václav Šmilauer <eudoxos@arcig.cz>
#include "UniaxialForcer.hpp"

#include <core/Aabb.hpp>
#include <core/InteractionContainer.hpp>
#include <core/Scene.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((UniaxialForcer));
/************************ UniaxialStrainer **********************/
CREATE_LOGGER(UniaxialForcer);

void UniaxialForcer::init()
{
	needsInit = false;

	assert(posId >=0);
	assert(negId >=0);
	

	const shared_ptr<Body>& bp = Body::byId(posId, scene);
	const shared_ptr<Body>& bn = Body::byId(negId, scene);
	if (blockDisplacements) {
		bn->state->blockedDOFs = State::DOF_ALL;
		bp->state->blockedDOFs = (State::DOF_X | State::DOF_Y | State::DOF_Z) ^ (1 << axis);
	}
	else {
		bn->state->blockedDOFs = State::DOF_ALL;
	}
	posCoords= bp->state->pos[axis];
	negCoords =bn->state->pos[axis];
	
	currentLength=originalLength = axisCoord(posId) - axisCoord(negId);
	LOG_DEBUG(
	        "Reference node: positive #" << posId << " at " << axisCoord(posId) << "; negative #" << negId << " at "
	                                          << axisCoord(negId));
	LOG_INFO("Setting initial length to " << originalLength << " (between #" << negId << " and #" << posId << ")");
	if (originalLength <= 0)
		throw runtime_error(
		        ("UniaxialForcer: Initial length is negative or zero (swapped reference particles?)! " + boost::lexical_cast<string>(originalLength))
		                .c_str());
	/* this happens is nan propagates from e.g. brefcom consitutive law in case 2 bodies have _exactly_ the same position
	 * (the the normal strain is 0./0.=nan). That is an user's error, however and should not happen. */
	if (math::isnan(originalLength)) throw logic_error("UniaxialForcer: Initial length is NaN!");
	
	assert(originalLength > 0 && !math::isnan(originalLength));
	
	assert(loadTime_s > 0 && loadForce > 0);
	if (loadTime_s< scene->dt)throw runtime_error("loadTime_s<scene->dt!");
	loadStepNum = long(loadTime_s / scene->dt) + 1;
	currentStrain = 0;
	for(int i=0;i<3;i++){
		unitForce[i]=(i==axis?1:0);
	}
	std::cout<<"unitforce is "<<unitForce<<std::endl;
}

void UniaxialForcer::action()
{
	if (needsInit) init();
	// postconditions for initParams
	assert(posId >=0 && negId >= 0 && originalLength > 0);
	//nothing to do
	if (posId <0 || negId <0) return;

	currentLength= axisCoord(posId) - axisCoord(negId);
	currentStrain = (currentLength - originalLength) / originalLength;
	Vector3r stretchForce = Vector3r(0, 0, 0);
	if (currentIter < loadStepNum) {
		stretchForce = (currentIter + 1) * loadForce/ loadStepNum *unitForce;
		currentIter += 1;
	}
	else {
		solve_num += 1;
		stretchForce = loadForce*unitForce;
		if (solve_num >= force_num_max) {
			scene->stopAtIter = scene->iter + 1;
			std::cerr << "thread stretch completed!" << std::endl;
			std::cerr << "thread originalLength =" << originalLength*1000 <<" mm"<< std::endl;
			std::cerr << "thread currentLength =" << currentLength*1000 <<" mm"<< std::endl;
			std::cerr << "thread stretched Length =" << currentLength*1000-originalLength*1000<<" mm"<< std::endl;
			std::cerr << "thread currentStrain =" << currentStrain<< std::endl;
		}

	}
	addForce(scene, posId, stretchForce);
}

} // namespace yade
