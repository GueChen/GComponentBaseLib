/**
 *  @file  	gjk_data_structures.h
 *  @brief 	this file define mutiple GJK/EPA relation data structures for computing.
 *  @author Gue Chen<guechen@buaa.edu.cn>
 *  @date 	Apr 7th, 2023
 **/
#ifndef __GJK_DATA_STRUCTURE_H
#define __GJK_DATA_STRUCTURE_H

#include <GComponent/types.h>

#include <cinttypes>
#include <vector>
#include <array>

namespace GComponent {

struct TriangleFace {
/// Type alias
	template<class T>
	using AdjType = std::array<T, 3>;
/*___________________________________________________________________*/
/// Methods
	TriangleFace() = default;

	TriangleFace(const uint32_t i0, const uint32_t i1, const uint32_t i2){
		m_indices = { i0, i1, i2 };
		m_adj_faces.fill(nullptr);
		m_adj_edges.fill(-1);
	}

	bool Link(const int e0, TriangleFace* face, const int e1);

	auto operator<=>(const TriangleFace& other) const {
		return m_plane_dist <=> other.m_plane_dist;
	}

/*___________________________________________________________________*/
/// Filds
	// self properties
	Vec3f				   m_norm;
	float				   m_plane_dist   = std::numeric_limits<float>::lowest();
	uint32_t			   m_idx;
	
	// adj properties
	AdjType<TriangleFace*> m_adj_faces;
	AdjType<int>		   m_adj_edges;
	AdjType<uint32_t>	   m_indices;		// index of vertices

	// control flags
	bool				   is_obsolete	= false;
	bool				   is_in_queue  = false;

};

struct Edge {
/*___________________________________________________________________*/
/// Methods
	inline uint32_t	  GetSourcePointIdx() const { return m_face->m_indices[m_idx]; }
	inline uint32_t	  GetEndPointIdx() const	{ return m_face->m_indices[(m_idx + 1) % 3]; }

/*___________________________________________________________________*/
/// Filds
	TriangleFace* m_face;
	uint32_t	  m_idx;
};

/*______________________Implementation ________________________________*/
bool TriangleFace::Link(const int e0, TriangleFace* face, const int e1) {
	this->m_adj_faces[e0] = face;
	this->m_adj_edges[e0] = e1;
	face->m_adj_faces[e1] = this;
	face->m_adj_edges[e1] = e0;
	return (m_indices[e0] == face->m_indices[(e1 + 1) % 3]) &&
		   (m_indices[(e0 + 1) % 3] == face->m_indices[e1]);
}

} // !namespace GComponent

#endif // !__GJK_DATA_STRUCTURE_H
