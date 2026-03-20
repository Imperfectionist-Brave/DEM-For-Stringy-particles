// 2010 © Václav Šmilauer <eudoxos@arcig.cz>

#include "BodyContainer.hpp"
#include "Body.hpp"
#include "Clump.hpp"
#include "Scene.hpp"
#include <lib/base/LoggingUtils.hpp>
#ifdef YADE_OPENMP
#include <omp.h>
#endif

namespace yade { // Cannot have #include directive inside.


CREATE_LOGGER(BodyContainer);

Body::id_t BodyContainer::insert(shared_ptr<Body> b)
{
	const shared_ptr<Scene>& scene = Omega::instance().getScene();
	b->iterBorn                    = scene->iter;
	b->timeBorn                    = scene->time;
	b->id                          = body.size();
	scene->doSort                  = true;

	if (enableRedirection) {
		insertedBodies.push_back(b->id);
		dirty             = true;
		checkedByCollider = false;
	}
	body.push_back(b);
	// Notify ForceContainer about new id
	scene->forces.addMaxId(b->id);
	return b->id;
}

Body::id_t BodyContainer::insertAtId(shared_ptr<Body> b, Body::id_t candidate)
{
	if (not b) LOG_ERROR("Inserting null body");
	const shared_ptr<Scene>& scene = Omega::instance().getScene();
	if (enableRedirection) {
		dirty             = true;
		checkedByCollider = false;
		useRedirection    = true;
		insertedBodies.push_back(candidate); /*realBodies.push_back(candidate); */
	}                                            // if that special insertion is used, switch to algorithm optimized for non-full body container
	if (unsigned(candidate) >= size()) {
		body.resize(candidate + 1, nullptr);
		scene->forces.addMaxId(candidate);
	} else if (body[candidate]) {
		LOG_ERROR("invalid candidate id: " << candidate);
		return -1;
	}

	b->iterBorn   = scene->iter;
	b->timeBorn   = scene->time;
	b->id         = candidate;
	body[b->id]   = b;
	scene->doSort = true;
	return b->id;
}

void BodyContainer::clear()
{
	body.clear();
	dirty             = true;
	checkedByCollider = false;
}

BodyContainer::iterator BodyContainer::begin() { return iterator(body.begin(), body.end()); }

BodyContainer::iterator BodyContainer::end() { return iterator(body.end(), body.end()); }

size_t BodyContainer::size() const { return body.size(); }

shared_ptr<Body>& BodyContainer::operator[](unsigned int id) { return body[id]; }
const shared_ptr<Body>& BodyContainer::operator[](unsigned int id) const { return body[id]; }

bool BodyContainer::exists(Body::id_t id) const { return ((id >= 0) && ((size_t)id < body.size()) && ((bool)body[id])); }

bool BodyContainer::eraseAlreadyLocked(Body::id_t id, bool eraseClumpMembers)
{ //default is false (as before)
	if (!body[id]) return false;
	const shared_ptr<Scene>& scene = Omega::instance().getScene();
	if (enableRedirection) {
		useRedirection    = true;
		dirty             = true;
		checkedByCollider = false;
		erasedBodies.push_back(id);
	} // as soon as a body is erased we switch to algorithm optimized for non-full body container
	const shared_ptr<Body>& b = Body::byId(id);
	if ((b) and (b->isClumpMember())) {
		const shared_ptr<Body>  clumpBody = Body::byId(b->clumpId);
		const shared_ptr<Clump> clump     = YADE_PTR_CAST<Clump>(clumpBody->shape);
		Clump::del(clumpBody, b);
		if (clump->members.size() == 0) this->eraseAlreadyLocked(clumpBody->id, false); //Clump has no members any more. Remove it
	}

	if ((b) and (b->isClump())) {
		//erase all members if eraseClumpMembers is true:
		const shared_ptr<Clump>&    clump   = YADE_PTR_CAST<Clump>(b->shape);
		std::map<Body::id_t, Se3r>& members = clump->members;
		std::vector<Body::id_t>     idsToRemove;
		for (const auto& mm : members)
			idsToRemove.push_back(mm.first); // Prepare an array of ids, which need to be removed.
		for (Body::id_t memberId : idsToRemove) {
			if (eraseClumpMembers) {
				this->eraseAlreadyLocked(memberId, false); // erase members
			} else {
				//when the last members is erased, the clump will be erased automatically, see above
				Body::byId(memberId)->clumpId = Body::ID_NONE; // make members standalones
			}
		}
		body[id].reset();
		return true;
	}

	for (auto it = b->intrs.begin(), end = b->intrs.end(); it != end;) { //Iterate over all body's interactions
		Body::MapId2IntrT::iterator willBeInvalid = it;
		++it;
		scene->interactions->erase(willBeInvalid->second->getId1(), willBeInvalid->second->getId2(), willBeInvalid->second->linIx);
	}
	b->id = -1; //else it sits in the python scope without a chance to be inserted again
	body[id].reset();
	return true;
}

bool BodyContainer::erase(Body::id_t id, bool eraseClumpMembers)
{
	// Fix https://gitlab.com/yade-dev/trunk/-/issues/235 - wait with updating until some body removal is finished.
	// NOTE: This is called very rarely, because it is erasing bodies. Not typical simulation occurrence.
	//       But we still have to make sure to avoid crashes. Even those which are rare :)
	const std::lock_guard<std::mutex> lock(drawloopmutex);

	return this->eraseAlreadyLocked(id, eraseClumpMembers);
}

void BodyContainer::updateRealBodies()
{
	if (not enableRedirection) {
		LOG_ONCE_WARN("updateRealBodies returns because enableRedirection is false - please report bug");
		return;
	}
	if (not dirty) return; //already ok
	unsigned long size1 = realBodies.size();
	realBodies.clear();
	realBodies.reserve((long unsigned)(size1 * 1.3));
#ifdef YADE_MPI
	unsigned long size2 = subdomainBodies.size();
	subdomainBodies.clear();
	subdomainBodies.reserve((long unsigned)(size2 * 1.3));
	const int& subdomain = Omega::instance().getScene()->subdomain;
#endif
	for (const auto& b : body) {
		if (not b) continue;
		realBodies.push_back(b->getId());
#ifdef YADE_MPI
		if (b->subdomain == subdomain and not b->getIsSubdomain()) subdomainBodies.push_back(b->id);
#endif
	}
	dirty = false;
}


} // namespace yade
namespace yade {
YADE_PLUGIN((NodeNeighbor));
	
void NodeNeighbor::initial() {
				num = bodyids.size();
				/*L.clear();
				Areas.clear();
				Norm.clear();*/
				if (num == 6) {
					Vector3r*L1=new Vector3r[7];
					Real *Areas1=new Real[6];
					Vector3r *Norm1=new Vector3r[6];
					L=L1;
					Areas=Areas1;
					Norm=Norm1;
					/*L.reserve(7);
					Areas.reserve(6);
					Norm.reserve(6);*/
				}
				else {
					Vector3r*L1=new Vector3r[4];
					Real *Areas1=new Real[3];
					Vector3r *Norm1=new Vector3r[3];
					L=L1;
					Areas=Areas1;
					Norm=Norm1;
					/*L.reserve(4);
					Areas.reserve(3);
					Norm.reserve(3);*/
				}
				UnitForce = Vector3r(0,0,0);
			}
			//YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(NodeNeighbor, Serializable, "Storage node information.",
			//	((int, id, -1, , "Current timestep for integration."))
			//	((int, mask, -1, , "Current timestep for integration."))//0 逆时针计算 1 顺时针计算法向量
			//	,
			//	/*ctor*/
			//	,
			//	/* py */
			//	);
void NodeNeighbor::computeL() {
				assert(bodyStates.size() >= 4);
				std::vector<shared_ptr<State>>::iterator it = bodyStates.begin();
				if (num == 6){
					for(int i=0;i<6;++i,++it){
						L[i]=(*it)->pos - selfstate->pos;
					}
					L[6]=L[0];
				}
				else {//num==4
					for(int i=0;i<4;++i,++it){
						L[i]=(*it)->pos - selfstate->pos;
					}
				}
				
				/*L.clear();
				std::vector<shared_ptr<State>>::iterator it = bodyStates.begin();
				if (num == 6) {//中间点  
					for (; it != bodyStates.end(); ++it) {
						L.push_back((*it)->pos - selfstate->pos);
					}
					L.push_back(L[0]);//7组l
				}
				else {//num==4
					for (; it != bodyStates.end(); ++it) {
						L.push_back((*it)->pos - selfstate->pos);
					}
				}
				return L.size();//6或4*/
			}

void NodeNeighbor::computeAreaAndNorm() {
				if(num==6){
					if (!mask) {
					  for (int i=0; i<6;++i) {
						Vector3r vertical = (L[i]).cross(L[i+1]);
						Real scanl = vertical.norm();
						Areas[i]=0.5*scanl;
						Norm[i]=-vertical / scanl;
						}
					}

				        else  {
					  for (int i=0; i<6;++i) {
						Vector3r vertical = (L[i]).cross(L[i+1]);
						Real scanl = vertical.norm();
						Areas[i]=0.5*scanl;
						Norm[i]=vertical / scanl;
						}
					}
				}else{
				
					if (!mask) {
					  for (int i=0; i<3;++i) {
						Vector3r vertical = (L[i]).cross(L[i+1]);
						Real scanl = vertical.norm();
						Areas[i]=0.5*scanl;
						Norm[i]=-vertical / scanl;
						}
					}

				        else  {
					  for (int i=0; i<3;++i) {
						Vector3r vertical = (L[i]).cross(L[i+1]);
						Real scanl = vertical.norm();
						Areas[i]=0.5*scanl;
						Norm[i]=vertical / scanl;
						}
					}
				
				}
				
				//*******************
				/*assert(L.size() >= 4 && Norm.size() >= 4);
				Areas.clear();
				Norm.clear();
				std::vector<Vector3r>::iterator it1 = L.begin();
				std::vector<Vector3r>::iterator it2 = std::next(it1);
				if (!mask) {
					for (; it2 != L.end(); ++it1, ++it2) {
						Vector3r vertical = (*it1).cross(*it2);
						Real scanl = vertical.norm();
						Areas.push_back(0.5*scanl);
						Norm.push_back(-vertical / scanl);
					}
				}

				else  {
					for (; it2 != L.end(); ++it1, ++it2) {
						Vector3r vertical = (*it1).cross(*it2);
						Real scanl = vertical.norm();
						Areas.push_back(0.5*scanl);
						Norm.push_back(-vertical / scanl);
					}
				}
				return Areas.size();//6或3*/
			}


Vector3r NodeNeighbor::computeUnitForce()
			{
				Vector3r Force = Vector3r(0,0,0);
				computeL();//7 4
				computeAreaAndNorm();
				//int a = computeL();//7 4
				//int b = computeAreaAndNorm();
				for (int i = 0; i < num; i++) Force += 1.0 / 3 * Areas[i] * Norm[i];
				UnitForce = Force;
				return Force;
				/*if (a == (b + 1)) {
					for (int i = 0; i < (a - 1); i++) {
						Force += 1.0 / 3 * Areas[i] * Norm[i];
					}
					UnitForce = Force;
					return Force;
				}
				
				else throw std::runtime_error("Neighbor L shoud greater than Area ");*/
			}
}
namespace yade {

YADE_PLUGIN((NodeNeighborContainer));

std::vector<Body::id_t> NodeNeighborContainer::insert([[maybe_unused]]int _cows,[[maybe_unused]]int _cols,std::vector<std::vector<Body::id_t>> Matrix) {  //传入膜节点二维矩阵

				Ids = Matrix;//直接赋值
				std::vector<std::vector<Body::id_t>>::iterator first = Ids.begin();
				/*std::vector<std::vector<Body::id_t>>::iterator last = Ids.end();
				std::vector<std::vector<Body::id_t>>::iterator Real = std::prev(last);*/
				cows = Ids.size();
				cols = (*first).size();
				assert((cows == _cows) && (cols = _cols));
				const shared_ptr<Scene>& scene = Omega::instance().getScene();
				this->iterBorn = scene->iter;
				this->timeBorn = scene->time;
				/*std::vector<Body::id_t>::iterator colsPointerTF = first.begin();
				std::vector<Body::id_t>::iterator colsPointerTE = first.end();
				for (; colsPointerTF != colsPointerTE; colsPointerTF++) {
					topNode.clear();
					topNode.push_back(*colsPointerTF);
				}*/
				Self.clear();
				topNode.clear();
				for (int i = 0; i < cols; i++) {
					
					topNode.push_back(Ids[0][i]);
					Self.push_back(Ids[0][i]);
				}
				freeNode.clear();
				for (int i = 1; i < cows-1; i++) {
					for (int j = 0; j < cols; j++) {
						freeNode.push_back(Ids[i][j]);
						Self.push_back(Ids[i][j]);
					}
				}
				bottomNode.clear();
				for (int i = 0; i < cols; i++) {

					bottomNode.push_back(Ids[cows - 1][i]);
					Self.push_back(Ids[cows - 1][i]);
				}
				Neighbor.reserve(cows*cols);//预存空间ids
				extend();//扩充
				allocat();//存入NodeNeighbor

				return Self;
			}

void NodeNeighborContainer::extend() {
				for (int i = 0; i < cows; i++) {
					Body::id_t first = Ids[i][0];
					Body::id_t end = Ids[i][cols-1];
					Ids[i].insert(Ids[i].begin(), end);//前插入
					Ids[i].push_back(first);//后插入
				}
				
			}

			
void NodeNeighborContainer::allocat() {
				const shared_ptr<Scene>& scene = Omega::instance().getScene();
				int Ecows = cows;
				int Ecols = cols+2;
				std::vector<Body::id_t> outid;//容器临时创建传值
				std::vector<shared_ptr<State>> outStates;
				outid.clear();
				outStates.clear();
				for (int i = 0; i < Ecows; i++) {
					for (int j = 1; j < Ecols-1;j++){//cols 列个数
						
						if (i == 0) {//
							outid.push_back(Ids[i][j + 1]);
							outid.push_back(Ids[i+1][j]);
							outid.push_back(Ids[i+1][j-1]);
							outid.push_back(Ids[i][j - 1]);
							outStates.push_back(((*(scene->bodies))[Ids[i][j + 1]])->state);//*(scene->bodies)[Ids[i][j + 1]] 
							outStates.push_back(((*(scene->bodies))[Ids[i + 1][j]])->state);
							outStates.push_back(((*(scene->bodies))[Ids[i + 1][j - 1]])->state);
							outStates.push_back(((*(scene->bodies))[Ids[i][j - 1]])->state);
							Body::id_t selfid = Ids[i][j];
							shared_ptr<State>selfstate=(*(scene->bodies))[selfid]->state;
							bool mask = true;
							shared_ptr<NodeNeighbor> node(new NodeNeighbor(selfid, selfstate, mask, outid, outStates));
							Neighbor.push_back(node);//使用值传递
								//con.push_back(shared_ptr<person>(new person(45, "we1")));
						}
						else if (i == Ecows - 1) {
							outid.push_back(Ids[i][j + 1]);
							outid.push_back(Ids[i - 1][(i%2==0 ?j:j+1)]);//偶数行[0] 奇数行时[1] Ids[i - 1][j+1]
							outid.push_back(Ids[i - 1][(i%2==0 ?j - 1:j)]);//偶数行 奇数行时 Ids[i - 1][j]
							outid.push_back(Ids[i][j - 1]);
							outStates.push_back(((*(scene->bodies))[Ids[i][j + 1]])->state);
							outStates.push_back(((*(scene->bodies))[Ids[i - 1][(i % 2 == 0 ? j : j + 1)]])->state);
							outStates.push_back(((*(scene->bodies))[Ids[i - 1][(i % 2 == 0 ? j - 1 : j)]])->state);
							outStates.push_back(((*(scene->bodies))[Ids[i][j - 1]])->state);
							Body::id_t selfid = Ids[i][j];
							shared_ptr<State>selfstate = (*(scene->bodies))[selfid]->state;
							bool mask = false;
							shared_ptr<NodeNeighbor> node(new NodeNeighbor(selfid, selfstate, mask, outid, outStates));
							Neighbor.push_back(node);
						}
						else {
							outid.push_back(Ids[i][j + 1]);
							outid.push_back(Ids[i - 1][(i % 2 == 0 ? j : j + 1)]);
							outid.push_back(Ids[i - 1][(i % 2 == 0 ? j-1 : j )]);
							outid.push_back(Ids[i][j - 1]);
							outid.push_back(Ids[i+1][(i % 2 == 0 ? j - 1:j)]);
							outid.push_back(Ids[i+1][(i % 2 == 0 ? j  : j+1)]);
							outStates.push_back(((*(scene->bodies))[Ids[i][j + 1]])->state);
							outStates.push_back(((*(scene->bodies))[Ids[i - 1][(i % 2 == 0 ? j : j + 1)]])->state);
							outStates.push_back(((*(scene->bodies))[Ids[i - 1][(i % 2 == 0 ? j - 1 : j)]])->state);
							outStates.push_back(((*(scene->bodies))[Ids[i][j - 1]])->state);
							outStates.push_back(((*(scene->bodies))[Ids[i+1][(i % 2 == 0 ? j - 1 : j)]])->state);
							outStates.push_back(((*(scene->bodies))[Ids[i+1][(i % 2 == 0 ? j : j + 1)]])->state);
							Body::id_t selfid = Ids[i][j];
							shared_ptr<State>selfstate = (*(scene->bodies))[selfid]->state;
							bool mask = false;
							shared_ptr<NodeNeighbor> node(new NodeNeighbor(selfid, selfstate, mask, outid, outStates));
							Neighbor.push_back(node);

						}
						outid.clear();
						outStates.clear();
					}
				}
				
			}

			
void NodeNeighborContainer::clear() {
				Neighbor.clear();//清空容器
				Self.clear();
				Ids.clear();
			}
NodeNeighborContainer::iterator NodeNeighborContainer::begin() { return Neighbor.begin(); }
NodeNeighborContainer::iterator NodeNeighborContainer::end() { return Neighbor.end(); }
int NodeNeighborContainer::size()const { return Neighbor.size(); }
bool NodeNeighborContainer::exists(Body::id_t id)const{
				auto it = std::find(Self.begin(), Self.end(), id);
				return ((id>=0)&&(it != Self.end()));
			}
}
