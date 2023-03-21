#include "GCube.h"

using namespace GCONST;

//TODO: 这里所有方盒子可通用一例
//		有待优化使用 { 单例模式 }
GCube::GCube()
{
	
	// 数据句柄生成
	std::tie(model_id_, VBO) = GMeshObject::genVABO((void *)VertexData, sizeof(VertexData));

	// 数据传递
	GMeshObject::EnableVertexAttrbArrays(3, 3, 2);

	//TODO: 这里感觉可以用RAII优化
	// Wrapper::HelperBinder可以RAII
	glBindVertexArray(0);

	//TODO: 可由此派生带Texture的盒子类
	// DiffuseTexureID = loadFromFile("xxx.jpg");
}

void GCube::Draw(MyShader& shader) const noexcept
{
	glBindVertexArray(model_id_);

	glDrawArrays(GL_TRIANGLES, 0, 36);
	
	glBindVertexArray(0);
}