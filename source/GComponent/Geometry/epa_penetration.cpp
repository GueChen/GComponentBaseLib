#include <GComponent/Geometry/gcollision_detection.h>
#include <GComponent/Geometry/gjk_data_structures.h>
#include <GBase/IDPool.h>

#include <queue>
#include <stack>
#include <algorithm>

// reference form NVIDIA PhysX GUEPA.cpp
namespace GComponent {

constexpr static const float kEpsilon = 1e-5;

struct TriangleFaceDistComparator {
	bool operator()(const TriangleFace* l, const TriangleFace* r) const {
		return *l > *r;
	}
};

constexpr const size_t MaxFaceNum	= 64;
constexpr const size_t MaxSupporNum = 64;
constexpr const size_t MaxEdgeNum   = 32;
class EPAImpl {
public:
	// triangle idx
	using Tri	  = std::array<uint32_t, 3>;
	using CRefTri = const Tri&;
public:
	EPAImpl(const GJKConvex&   convex_a, const GJKConvex& convex_b, 
			std::vector<Vec3f> polytope,
			std::vector<Vec3f> polytope_a,
			std::vector<Vec3f> polytope_b):
		a_(convex_a),
		b_(convex_b),
		polytope_(std::move(polytope)),
		poly_a_(std::move(polytope_a)),
		poly_b_(std::move(polytope_b))
	{}

	GJKStatus ComputePenetration(GJKOutput& output);

private:
	void GetContactInformation(GJKOutput& output, TriangleFace& tri);

	/// <summary>
	/// Found the closest point [closest] to origin on triangle face {tri}, expressed as [alpha, beta, gamma] 
	/// which knowned as barycentric coordinates. With these coordinates mapping to get the closest points pair
	/// on A, B geoemtries.
	/// <para>
	/// 找到三角面片 {tri} 上距离原点上的最近点 [closest], 获取它的三角重心表达式 [alpha, beta, gamma]. 使用该系数从 A，
	/// B 的支持点集中获取构筑最近点的一个点对，该点对即两个碰撞体间的最大穿透点.
	/// </para>
	/// </summary>
	/// <param name="tri"></param>
	/// <param name="closest_a"></param>
	/// <param name="closest_b"></param>
	void GetClosestPoint(TriangleFace& tri, Vec3f& pa, Vec3f& pb);

	/// <summary>
	/// 对 A, B 两个多面体分别做 support 并将结果保存在 minkowski 的多面体中
	/// </summary>
	/// <param name="dir"></param>
	void DoSupport		(const Vec3f& dir);

	/// <summary>
	/// In case the initialize polytope possess only one point, try to expanding it with a unit_x direction with support mapping
	/// <para>
	/// 在初始多面体近含一个点时，通过对多面体使用单位 x 方向的支撑映射来扩展多面体
	/// </para>
	/// </summary>
	/// <param name="num_vert"></param>
	/// <param name="upper_bound"></param>
	/// <returns></returns>
	bool ExpandPoint	(uint32_t& num_vert, float upper_bound);

	/// <summary>
	/// using current minkowski difference shape, find the most close to orthogonal axis direction
	/// <para>
	/// 使用当前线段的 Minkowski 差，找到与正交状态最接近的单位轴向作为支撑映射的搜索方向
	/// </para>
	/// </summary>
	/// <param name="num_vert"></param>
	/// <param name="upper_bound"></param>
	/// <returns></returns>
	bool ExpandSegment	(uint32_t& num_vert, float upper_bound);

	/// <summary>
	/// Creating a overlap shape as the initialize polytope
	/// <para>
	/// 使用重叠的正反面三角形作为初始的多面体图形
	/// </para>
	/// </summary>
	/// <param name="num_vert"></param>
	/// <param name="upper_bound"></param>
	/// <returns></returns>
	bool ExpandTriangle	(uint32_t& num_vert, float upper_bound);

	/// <summary>
	/// 通过顶点索引添加三角面至优先队列中
	/// </summary>
	/// <param name="tri"></param>
	/// <param name="upper"></param>
	/// <returns></returns>
	TriangleFace* AddFace(CRefTri tri, float upper);

	/// <summary>
	/// validate added face is normal and can add to queue structure
	/// </summary>
	/// <param name="tri"></param>
	/// <param name="tri_idx"></param>
	/// <param name="upper"></param>
	/// <returns></returns>
	bool TriangleValid(TriangleFace& tri, CRefTri tri_idx, float upper);
	
	void Silhouette(const Vec3f& p, TriangleFace& tri, const uint32_t idx);

private:
	const GJKConvex&   a_, 
				   &   b_;
	std::vector<Vec3f> poly_a_;
	std::vector<Vec3f> poly_b_;
	std::vector<Vec3f> polytope_;
	std::priority_queue<TriangleFace*, std::vector<TriangleFace*>, TriangleFaceDistComparator>
					   queue_;
	IDPoolBase		   face_manager_;
	TriangleFace	   faces_buffer_[MaxFaceNum];
	std::vector<Edge>  edges_buffer_;
};

GJKStatus EPAImpl::ComputePenetration(GJKOutput& output) {
	enum SimplexCase {
		Point = 1, Line, Triangle, Tetrhetron
	} initial_type = static_cast<SimplexCase>(poly_a_.size());
	
	float	 upper_bound = std::numeric_limits<float>::max();
	uint32_t num_verts	= 0;

	switch (initial_type) {
	case Point:{
		if (not ExpandPoint(num_verts,   upper_bound)) return EPA_FAIL;
		break;
	}
	case Line:{
		if (not ExpandSegment(num_verts, upper_bound)) return EPA_FAIL;
		break;
	}
	case Triangle: {
		if (not ExpandTriangle(num_verts, upper_bound))return EPA_FAIL;
		break;
	}
	case Tetrhetron: {
		// checking the index order, make the face normal outwards, shuffle if neccessary		
		const Vec3f l_01 = polytope_[1] - polytope_[0],
					l_02 = polytope_[2] - polytope_[0],
					l_03 = polytope_[3] - polytope_[0];
		const Vec3f norm = l_01.cross(l_02).normalized();
		if (norm.dot(l_03) > 0) {	
			std::swap(poly_a_[1], poly_a_[2]);
			std::swap(poly_b_[1], poly_b_[2]);
			std::swap(polytope_[1], polytope_[2]);
		}

		TriangleFace* f0 = AddFace({0, 1, 2}, upper_bound),
					* f1 = AddFace({0, 3, 1}, upper_bound),
					* f2 = AddFace({0, 2, 3}, upper_bound),
					* f3 = AddFace({1, 3, 2}, upper_bound);
		
		if (queue_.empty()) return EPA_FAIL;

		f0->Link(0, f1, 2); f0->Link(1, f3, 2); f0->Link(2, f2, 0);
		f1->Link(0, f2, 2); f1->Link(1, f3, 0);
		f2->Link(1, f3, 1);

		num_verts = 4;
	}
	}
	TriangleFace* face = nullptr;
	do {
		// closest face from queue
		face = queue_.top(); queue_.pop();
		face->is_in_queue = false;

		if (not face->is_obsolete) {
			const Vec3f face_norm = face->m_norm;
			const float face_dist = face->m_plane_dist;

			DoSupport(-face_norm);
			
			const float dist = polytope_.back().dot(face_norm);

			// found the minkowski difference border
			if (std::abs(dist - face_dist) <= kEpsilon) {
				GetContactInformation(output, *face);
				const Vec3f diff = output.closest_a - output.closest_b;				
				if (diff.norm() > std::abs(output.depth)) {
					return EPA_DEGENERATE;
				}				
				return EPA_CONTACT;
			}
			

			upper_bound = std::min(upper_bound, dist);

			// compute silhouette 
			edges_buffer_.clear();

			// flood fill silhouette
			face->is_obsolete = true;
			for (int i = 0; i < 3; ++i) {
				Silhouette(polytope_.back(), *face->m_adj_faces[i], face->m_adj_edges[i]);
			}

			// checking edge buffer is Valid 
			if (edges_buffer_.empty() || edges_buffer_.size() > MaxEdgeNum) {				
				GetContactInformation(output, *face);
				return EPA_DEGENERATE;
			}
			
			// checking edge buffer remianing is enough
			if (edges_buffer_.size() > MaxFaceNum - face_manager_.GetUsedID()) {
				GetContactInformation(output, *face);
				return EPA_DEGENERATE;
			}

			const uint32_t idx = num_verts++;
			Edge& edge		   = edges_buffer_.front();
			TriangleFace* first_face = AddFace({edge.GetEndPointIdx(),
												edge.GetSourcePointIdx(),
												idx},
												upper_bound);
			first_face->Link(0, edge.m_face, edge.m_idx);
			
			TriangleFace* last_face = first_face;
			for (int i = 1; i < edges_buffer_.size(); ++i) {
				Edge&		  edge	   = edges_buffer_[i];				
				TriangleFace* new_face = AddFace({ edge.GetEndPointIdx(), edge.GetSourcePointIdx(), idx},
												upper_bound);
				new_face->Link(0, edge.m_face, edge.m_idx);
				new_face->Link(2, last_face, 1);
				last_face = new_face;				
			}
			first_face->Link(2, last_face, 1);
			face_manager_.Free(face->m_idx);
		}
	}while(queue_.size() > 0 &&								// face not empty
		   upper_bound >= queue_.top()->m_plane_dist &&		// upper_bound always greater than rest plane
		   num_verts != MaxSupporNum);						// support point num not greter than max limitation
	
	GetContactInformation(output, *face);
	return EPA_DEGENERATE;
}

void EPAImpl::GetContactInformation(GJKOutput& output, TriangleFace& face){	
	GetClosestPoint(face, output.closest_a, output.closest_b);
	output.normal = -face.m_norm;
	output.depth  = -std::abs(face.m_plane_dist);
}


void EPAImpl::GetClosestPoint(TriangleFace& tri, Vec3f& closest_a, Vec3f& closest_b)
{
	const Vec3f closest = tri.m_plane_dist * tri.m_norm;
	
	const Vec3f p0 = polytope_[tri.m_indices[0]],
				p1 = polytope_[tri.m_indices[1]],
				p2 = polytope_[tri.m_indices[2]];
	const Vec3f v0 = p1 - p0, v1 = p2 - p0, v2 = closest - p0;

	// caculate barycentric coordinates
	float d00 = v0.dot(v0), d01 = v0.dot(v1),
		  d11 = v1.dot(v1),
		  d20 = v2.dot(v0), d21 = v2.dot(v1);

	float det	= d00 * d11 - d01 * d01;
	if (det >= kEpsilon) {		
		float beta   = (d11 * d20 - d01 * d21) / det;
		float gamma  = (d00 * d21 - d01 * d20) / det;
		float alpha  = 1 - beta - gamma;
		closest_a = alpha * poly_a_[tri.m_indices[0]] + beta * poly_a_[tri.m_indices[1]] + gamma * poly_a_[tri.m_indices[2]];
		closest_b = alpha * poly_b_[tri.m_indices[0]] + beta * poly_b_[tri.m_indices[1]] + gamma * poly_b_[tri.m_indices[2]];
	}
	else {		
		closest_a = poly_a_[tri.m_indices[0]];
		closest_b = poly_b_[tri.m_indices[0]];
	}
}

void EPAImpl::DoSupport(const Vec3f& dir)
{
	poly_a_  .push_back(a_.Support(-dir));
	poly_b_  .push_back(b_.Support( dir));
	polytope_.push_back(poly_a_.back() - poly_b_.back());
}

bool EPAImpl::ExpandPoint(uint32_t& num_vert, float upper_bound)
{
	const Vec3f dir = Vec3f::UnitX();
		
	DoSupport(dir);
	
	// guaranteeing not degenerating to point case after expanding
	if ((polytope_.back() - polytope_[0]).squaredNorm() < kEpsilon) {
		return false;
	}

	return ExpandSegment(num_vert, upper_bound);
}

bool EPAImpl::ExpandSegment(uint32_t& num_vert, float upper_bound)
{
	enum Dir { X = 0, Y, Z };
	const Vec3f v = polytope_[1] - polytope_[0];
	
	float min_abs = std::abs(v.x());
	int	  min_idx = 0;
	for (int i = 1; i < 3; ++i) {
		float cur_abs = std::abs(v(i));
		if (cur_abs < min_abs) {
			min_abs = cur_abs;
			min_idx = i;
		}
	}	
	
	switch (min_idx) {
	case X: DoSupport(v.cross(Vec3f::UnitX()).normalized()); break;
	case Y: DoSupport(v.cross(Vec3f::UnitY()).normalized()); break;
	case Z: DoSupport(v.cross(Vec3f::UnitZ()).normalized());
	}		
		
	return ExpandTriangle(num_vert, upper_bound);
}

bool EPAImpl::ExpandTriangle(uint32_t& num_vert, float upper_bound)
{
	num_vert = 3;
	
	TriangleFace* f0 = AddFace({0, 1, 2}, upper_bound);
	TriangleFace* f1 = AddFace({1, 0, 2}, upper_bound);
	
	if (queue_.empty()) {
		return false;
	}
	else {
		f0->Link(0, f1, 0); f0->Link(1, f1, 2); f0->Link(2, f1, 1);
		return true;
	}
}

TriangleFace* EPAImpl::AddFace(CRefTri tri, float upper)
{
	//assert(std::unique(tri.cbegin(), tri.cend()) != tri.cend() && "duplicate ele");
	assert(face_manager_.GetUsedID() < MaxFaceNum);
	
	const size_t  face_id = face_manager_.Get();
	TriangleFace* face	  = new(&faces_buffer_[face_id])TriangleFace(tri[0], tri[1], tri[2]);
	face->m_idx = face_id;
	
	if (TriangleValid(*face, tri, upper)) {
		queue_.push(face);
		face->is_in_queue = true;
	}
	else {
		face->is_in_queue = false;
	}

	return face;
}


bool EPAImpl::TriangleValid(TriangleFace& tri, CRefTri tri_idx, float upper)
{
	std::array<Vec3f, 3> mink_tri;
	for (int i = 0; i < tri_idx.size(); ++i) {
		mink_tri[i] = polytope_[tri_idx[i]];
	}
	
	const Vec3f l_01 = mink_tri[1] - mink_tri[0],
				l_02 = mink_tri[2] - mink_tri[0];
	const Vec3f norm = l_01.cross(l_02).normalized();

	if (norm.squaredNorm() < kEpsilon) return false;

	float sign_dist = norm.dot(mink_tri[0]);	
	if (sign_dist < 0) {
		tri.m_norm		 = -norm;
		tri.m_plane_dist = -sign_dist;
	}
	else {
		tri.m_norm		 = norm;
		tri.m_plane_dist = sign_dist;
	}
		
	return sign_dist <= upper;
}

void EPAImpl::Silhouette(const Vec3f& w, TriangleFace& tri, const uint32_t idx)
{
	std::stack<Edge> edge_stack;
	
	auto try_emplace_edge = [&stack = edge_stack](auto&& faces, auto&& edges, int idx) {		
		if (faces[idx]->is_obsolete) {			
			stack.emplace(faces[idx], edges[idx]);
		}
	};

	edge_stack.emplace(&tri, idx);

	while (!edge_stack.empty()) {		
		auto [face, index] = edge_stack.top(); edge_stack.pop();		
		assert(face->m_adj_faces[0] && face->m_adj_faces[1] && face->m_adj_faces[2]);

		if (!face->is_obsolete) {
			Vec3f mink_edge = w - polytope_[face->m_indices[0]];
			const float point_plane_dist = mink_edge.dot(face->m_norm);

			auto& adj_faces = face->m_adj_faces;
			auto& adj_edges = face->m_adj_edges;

			if (point_plane_dist <= 0) {				
				edges_buffer_.emplace_back(face, index);
			}
			else {
				face->is_obsolete = true;
			
				try_emplace_edge(adj_faces, adj_edges, (index + 1) % 3);
				try_emplace_edge(adj_faces, adj_edges, (index + 2) % 3);
				
				if (!face->is_in_queue) {
					face_manager_.Free(face->m_idx);
				}
			}
		}
	}
}

/*__________________________________OUTSIDE INTERFACE __________________________________________*/
GJKStatus GComponent::PenetrationEPA(GJKOutput& output, const GJKConvex& a, const GJKConvex& b, 
									 const std::vector<Vec3f>& simplex, 
									 const std::vector<Vec3f>& simplex_a, 
									 const std::vector<Vec3f>& simplex_b)
{	
	return EPAImpl(a, b, simplex, simplex_a, simplex_b).ComputePenetration(output);
}

}