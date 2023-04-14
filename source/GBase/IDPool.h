/**
 *  @file  	IDPool.h
 *  @brief 	This class is using for Id allocation.
 *  @author Gue Chen<guechen@buaa.edu.cn>
 *  @date 	Apr 7th, 2023
 **/
#ifndef __G_ID_POOL_H
#define __G_ID_POOL_H

#include <vector>

namespace GComponent {

class IDPoolBase {
public:
	IDPoolBase() = default;
	void Free(size_t id) {
		if (id == (current_id_ - 1)) --current_id_;
		else frees_.push_back(id);
	}

	void Clear() {
		current_id_ = 0;
		frees_.clear();
	}

	size_t Get() {
		const size_t size = frees_.size();
		
		if (size) {
			size_t ret = frees_.back();
			frees_.pop_back();
			return ret;
		}
		else {
			return current_id_++;
		}
	}

	inline size_t GetUsedID() const  { return current_id_ - frees_.size(); }

	inline size_t GetMaxUsed() const { return current_id_; }

protected:
	size_t				current_id_ = 0;
	std::vector<size_t> frees_		= {};
};

} // !namespace GComponent

#endif // !__G_ID_POOL_H
