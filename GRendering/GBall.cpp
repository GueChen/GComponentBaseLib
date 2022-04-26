#include "GBall.h"
using namespace PhsycalLagrangian;

GBall::GBall(vec3 o, float r):
	center(o), radius(r)
{
	setupMesh(36, 72);
	std::tie(ID, VBO) = GMeshObject::genVABO(&mesh.vertices[0], sizeof(Vertex) * mesh.vertices.size());
	GMeshObject::EnableVertexAttrbArrays(3, 3, 2);
	EBO = GMeshObject::genEBO(mesh.indicies);
	TextureID = GMeshObject::loadFromFile("Resources/world.jpg");
}

void GBall::setupMesh(int la, int lo) {
	const auto dLati = pi<float>() / la;
	const auto dLongi = 2.0f * pi<float>() / lo;
	auto deltaV = 1.0f / la;
	auto deltaU = 1.0f / ( lo - 1 );
	for (int i = 0; i <= la; ++i)
	{
		auto theta = -pi<float>() + dLati * i;
		auto sTh = sin(theta);
		auto cTh = cos(theta);
		auto V = 1.0f - deltaV * i;
		for (int j = 0; j < lo; ++j)
		{
			auto U = 1.0f - deltaU * j;
			auto phi = dLongi * j;
			auto sPh = sin(phi);
			auto cPh = cos(phi);
			vec3 dir = { sTh * cPh, cTh, sTh * sPh };
			auto vertex = Vertex{ radius * dir, dir, {U, V}};
			mesh.vertices.push_back(vertex);
		}
	}
	{
		// 填充球面
		for (int i = 0; i <= la; ++i)
		{
			auto lineAdd = i * lo;
			for (int j = 0; j < lo - 1; ++j)
			{
				auto pre = lineAdd + j;
				auto cross = pre + lo + 1;
				mesh.indicies.push_back({ pre ,pre + 1, cross - 1 });
				mesh.indicies.push_back({ pre + 1, cross , cross - 1 });
			}
		}
		// 填充球面的缝合处
		for (int i = 0; i <= la; ++i)
		{
			auto aft = lo * i;
			auto pre = aft + lo - 1;
			mesh.indicies.push_back({ pre, aft, pre + lo });
			mesh.indicies.push_back({ aft, aft + lo, pre + lo });
		}
	}
}

void GBall::Draw(MyShader& shader) const noexcept
{
	shader.setMat4("model", translate(identity<mat4>(), center));
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, TextureID);
	glBindVertexArray(ID);
	glDrawElements(GL_TRIANGLES, 3* mesh.indicies.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}